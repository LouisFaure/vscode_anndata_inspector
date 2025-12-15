/**
 * H5AD File Reader using system hdf5-tools (h5ls, h5dump)
 * 
 * Reads AnnData h5ad files and extracts experimental design information.
 * Uses system tools to avoid memory limits and improve performance on large files.
 */

import * as cp from 'child_process';
import * as fs from 'fs';
import * as path from 'path';
import * as os from 'os';
import * as util from 'util';
import * as crypto from 'crypto';

const execFile = util.promisify(cp.execFile);

export interface FactorInfo {
    name: string;
    categories: string[];
    counts: number[];
    type: 'experimental' | 'replicate' | 'classification' | 'batch';
}

export interface MatrixInfo {
    name: string;
    type: 'integer' | 'float';
}

export interface ObsmInfo {
    name: string;
    columns: number;
}

export interface H5ADMetadata {
    totalCells: number;
    factors: Map<string, FactorInfo>;
    species: string;
    geneCount: number;
    matrices: MatrixInfo[];
    obsm: ObsmInfo[];
}

// Helper to manage temp files
const tempFiles: Set<string> = new Set();

function getTempFilePath(prefix: string): string {
    const tmpDir = os.tmpdir();
    const name = `${prefix}_${crypto.randomBytes(6).toString('hex')}`;
    const p = path.join(tmpDir, name);
    tempFiles.add(p);
    return p;
}

function cleanupTempFile(p: string) {
    try {
        if (fs.existsSync(p)) fs.unlinkSync(p);
        tempFiles.delete(p);
    } catch (e) {}
}

// Ensure cleanup on exit
process.on('exit', () => {
    for (const file of tempFiles) {
        try { if (fs.existsSync(file)) fs.unlinkSync(file); } catch (e) {}
    }
});

async function runCommand(command: string, args: string[]): Promise<string> {
    try {
        // Increase buffer for large outputs (e.g. h5ls -r on complex files)
        const { stdout } = await execFile(command, args, { maxBuffer: 1024 * 1024 * 50 });
        return stdout;
    } catch (error: any) {
        throw new Error(`Command failed: ${command} ${args.join(' ')}\n${error.message}`);
    }
}

async function ensureLocalFile(filePath: string): Promise<{ path: string, isTemp: boolean }> {
    if (filePath.startsWith('http://') || filePath.startsWith('https://')) {
        const tempPath = getTempFilePath('h5ad_download');
        const response = await fetch(filePath);
        if (!response.ok) throw new Error(`Failed to fetch ${filePath}`);
        const buffer = await response.arrayBuffer();
        fs.writeFileSync(tempPath, Buffer.from(buffer));
        return { path: tempPath, isTemp: true };
    }
    return { path: filePath, isTemp: false };
}

interface DatasetInfo {
    shape: number[];
    typeClass: 'string' | 'integer' | 'float' | 'other';
    typeSize: number; // in bytes
    isSigned: boolean;
}

async function getDatasetInfo(filePath: string, datasetPath: string): Promise<DatasetInfo> {
    // Remove -d to avoid dumping data and exceeding buffer
    const stdout = await runCommand('h5ls', ['-v', `${filePath}/${datasetPath}`]);
    
    // Parse Shape
    // Output example: "    Dataset {10/10, 5/5}" or "    Dataset {10}"
    const shapeMatch = stdout.match(/Dataset \{([0-9, \/]+)\}/);
    const shape = shapeMatch 
        ? shapeMatch[1].split(',').map(s => parseInt(s.split('/')[0].trim()))
        : [];

    // Parse Type
    // Output example: "    Type:      Integer, 64-bit, little-endian"
    // "    Type:      String, length = 10, null pad, ascii character set"
    // "    Type:      native double"
    // "    Type:      native int"
    const typeMatch = stdout.match(/Type:\s+(.+)$/m);
    const typeStr = (typeMatch ? typeMatch[1] : '').toLowerCase();

    let typeClass: DatasetInfo['typeClass'] = 'other';
    let typeSize = 1;
    let isSigned = true;

    if (typeStr.includes('integer') || typeStr.includes('int')) {
        typeClass = 'integer';
        if (typeStr.includes('8-bit')) typeSize = 1;
        else if (typeStr.includes('16-bit')) typeSize = 2;
        else if (typeStr.includes('32-bit')) typeSize = 4;
        else if (typeStr.includes('64-bit')) typeSize = 8;
        else typeSize = 4; // Default native integer

        if (typeStr.includes('unsigned')) isSigned = false;
    } else if (typeStr.includes('float') || typeStr.includes('double')) {
        typeClass = 'float';
        if (typeStr.includes('32-bit')) typeSize = 4;
        else if (typeStr.includes('64-bit')) typeSize = 8;
        else typeSize = 4;
    } else if (typeStr.includes('string')) {
        typeClass = 'string';
        const lenMatch = typeStr.match(/length = (\d+)/);
        if (lenMatch) typeSize = parseInt(lenMatch[1]);
        else typeSize = 0; // Variable length?
    }

    return { shape, typeClass, typeSize, isSigned };
}

async function readDatasetBinary(filePath: string, datasetPath: string, info: DatasetInfo): Promise<Buffer> {
    const tempBin = getTempFilePath('h5dump_bin');
    try {
        // -d dataset -b LE (Little Endian) -o output
        await runCommand('h5dump', ['-d', datasetPath, '-b', 'LE', '-o', tempBin, filePath]);
        if (fs.existsSync(tempBin)) {
            return fs.readFileSync(tempBin);
        }
        return Buffer.alloc(0);
    } finally {
        cleanupTempFile(tempBin);
    }
}

/**
 * List all categorical factors in the /obs section of an h5ad file
 */
export async function listFactors(filePath: string): Promise<string[]> {
    const { path: localPath, isTemp } = await ensureLocalFile(filePath);
    try {
        // Optimize: Only list /obs to avoid scanning the whole file
        // h5ls -r file/obs
        let stdout = '';
        let isObsScan = false;
        try {
            stdout = await runCommand('h5ls', ['-r', localPath + '/obs']);
            isObsScan = true;
        } catch (e) {
            // Fallback to full scan if /obs direct access fails
            stdout = await runCommand('h5ls', ['-r', localPath]);
        }

        const lines = stdout.split('\n');
        
        // Find groups in obs that have both categories and codes
        // Map: groupPath -> Set of children
        const groups = new Map<string, Set<string>>();
        
        for (const line of lines) {
            // Line format: /path/to/obj  Type
            const parts = line.trim().split(/\s+/);
            if (parts.length < 2) continue;
            
            let objPath = parts[0];
            
            if (isObsScan) {
                // If we listed /obs, paths are relative to it (e.g. /factor or factor)
                // We need to reconstruct full path /obs/factor
                if (objPath.startsWith('/')) {
                    // /factor -> /obs/factor
                    // But check if it already starts with /obs/ (unlikely based on observation)
                    if (!objPath.startsWith('/obs/')) {
                        objPath = '/obs' + objPath;
                    }
                } else {
                    // factor -> /obs/factor
                    objPath = '/obs/' + objPath;
                }
            }

            if (objPath.startsWith('/obs/')) {
                const parent = path.dirname(objPath);
                const name = path.basename(objPath);
                
                if (!groups.has(parent)) groups.set(parent, new Set());
                groups.get(parent)?.add(name);
            }
        }

        const factors: string[] = [];
        for (const [groupPath, children] of groups) {
            if (children.has('categories') && children.has('codes')) {
                // groupPath is like /obs/factorName
                const factorName = path.basename(groupPath);
                factors.push(factorName);
            }
        }
        return factors.sort();

    } finally {
        if (isTemp) cleanupTempFile(localPath);
    }
}

/**
 * Extract category names from a categorical factor
 */
export async function extractCategories(filePath: string, factorName: string): Promise<string[]> {
    const { path: localPath, isTemp } = await ensureLocalFile(filePath);
    try {
        const datasetPath = `/obs/${factorName}/categories`;
        const info = await getDatasetInfo(localPath, datasetPath);
        
        if (info.typeClass === 'string' && info.typeSize > 0) {
            // Fixed length strings
            const buffer = await readDatasetBinary(localPath, datasetPath, info);
            const categories: string[] = [];
            const count = info.shape[0] || 0;
            
            for (let i = 0; i < count; i++) {
                const start = i * info.typeSize;
                const end = start + info.typeSize;
                if (end <= buffer.length) {
                    // Read string and trim nulls
                    let str = buffer.toString('utf8', start, end);
                    // Remove null termination and padding
                    str = str.replace(/\0/g, '');
                    categories.push(str);
                }
            }
            return categories;
        } else {
            // Variable length or other: fallback to text dump
            // h5dump -d path -y -w 0
            const stdout = await runCommand('h5dump', ['-d', datasetPath, '-y', '-w', '0', localPath]);
            // Parse text output
            // Look for strings in quotes
            const matches = stdout.matchAll(/"(.*?)"/g);
            const categories: string[] = [];
            for (const m of matches) {
                if (m[1] !== 'HDF5') { // Skip header if matched accidentally
                     categories.push(m[1]);
                }
            }
            // Filter out header noise if any
            // Usually data block is inside DATA { ... }
            // Simple regex might be risky but h5dump output is structured.
            // Better: extract content inside DATA { ... } first.
            const dataBlock = stdout.match(/DATA \{([\s\S]*?)\}/);
            if (dataBlock) {
                const inner = dataBlock[1];
                const innerMatches = inner.matchAll(/"(.*?)"/g);
                return Array.from(innerMatches, m => m[1]);
            }
            return [];
        }
    } finally {
        if (isTemp) cleanupTempFile(localPath);
    }
}

/**
 * Count cells per category for a factor
 */
export async function countCells(filePath: string, factorName: string): Promise<Map<string, number>> {
    const { path: localPath, isTemp } = await ensureLocalFile(filePath);
    try {
        const categories = await extractCategories(localPath, factorName);
        const counts = new Map<string, number>();
        categories.forEach(c => counts.set(c, 0));

        const codesPath = `/obs/${factorName}/codes`;
        const info = await getDatasetInfo(localPath, codesPath);
        
        if (info.typeClass === 'integer') {
            const buffer = await readDatasetBinary(localPath, codesPath, info);
            
            // Create typed array view
            let view: ArrayLike<number>;
            if (info.typeSize === 1) view = info.isSigned ? new Int8Array(buffer.buffer) : new Uint8Array(buffer.buffer);
            else if (info.typeSize === 2) view = info.isSigned ? new Int16Array(buffer.buffer) : new Uint16Array(buffer.buffer);
            else if (info.typeSize === 4) view = info.isSigned ? new Int32Array(buffer.buffer) : new Uint32Array(buffer.buffer);
            else view = new Int32Array(buffer.buffer); // Fallback for 64-bit (might lose precision but codes are usually small)

            for (let i = 0; i < view.length; i++) {
                const code = view[i];
                if (code >= 0 && code < categories.length) {
                    const cat = categories[code];
                    counts.set(cat, (counts.get(cat) || 0) + 1);
                }
            }
        }
        return counts;
    } finally {
        if (isTemp) cleanupTempFile(localPath);
    }
}

/**
 * Get total cell count from the h5ad file
 */
export async function getTotalCells(filePath: string): Promise<number> {
    const { path: localPath, isTemp } = await ensureLocalFile(filePath);
    try {
        // Try obs/_index
        try {
            const info = await getDatasetInfo(localPath, '/obs/_index');
            if (info.shape.length > 0) return info.shape[0];
        } catch (e) {}

        // Try X
        try {
            const info = await getDatasetInfo(localPath, '/X');
            if (info.shape.length > 0) return info.shape[0];
        } catch (e) {}

        return 0;
    } finally {
        if (isTemp) cleanupTempFile(localPath);
    }
}

/**
 * Detect species from gene naming patterns in /var
 */
export async function detectSpecies(filePath: string): Promise<string> {
    const { path: localPath, isTemp } = await ensureLocalFile(filePath);
    try {
        const possiblePaths = ['/var/_index', '/var/index', '/var/gene_symbols', '/var/gene_names'];
        let geneNames: string[] = [];

        for (const varPath of possiblePaths) {
            try {
                // Read first 50 items
                // h5dump -d path -s 0 -c 50 -y -w 0
                const stdout = await runCommand('h5dump', ['-d', varPath, '-s', '0', '-c', '50', '-y', '-w', '0', localPath]);
                const dataBlock = stdout.match(/DATA \{([\s\S]*?)\}/);
                if (dataBlock) {
                    const inner = dataBlock[1];
                    // Match strings
                    const matches = Array.from(inner.matchAll(/"(.*?)"/g), m => m[1]);
                    if (matches.length > 0) {
                        geneNames = matches;
                        break;
                    }
                }
            } catch (e) {}
        }

        if (geneNames.length === 0) return 'unknown';

        // Analyze naming patterns
        let uppercase = 0;
        let titlecase = 0;
        let lowercase = 0;

        for (const gene of geneNames) {
            if (/^[A-Z][A-Z0-9]+$/.test(gene)) uppercase++;
            else if (/^[A-Z][a-z0-9]+$/.test(gene)) titlecase++;
            else if (/^[a-z][a-z0-9]+$/.test(gene)) lowercase++;
        }

        const total = uppercase + titlecase + lowercase;
        if (total === 0) return 'unknown';

        if (uppercase / total > 0.5) return 'human';
        else if (titlecase / total > 0.5) return 'mouse';
        else if (lowercase / total > 0.5) return 'zebrafish';

        return 'unknown';
    } finally {
        if (isTemp) cleanupTempFile(localPath);
    }
}

/**
 * Get gene count from the h5ad file
 */
export async function getGeneCount(filePath: string): Promise<number> {
    const { path: localPath, isTemp } = await ensureLocalFile(filePath);
    try {
        try {
            const info = await getDatasetInfo(localPath, '/var/_index');
            if (info.shape.length > 0) return info.shape[0];
        } catch (e) {}

        try {
            const info = await getDatasetInfo(localPath, '/X');
            if (info.shape.length > 1) return info.shape[1];
        } catch (e) {}

        return 0;
    } finally {
        if (isTemp) cleanupTempFile(localPath);
    }
}

/**
 * Read first 5 values from a dataset to check for integer/float content
 */
async function readFirstFiveValues(filePath: string, datasetPath: string): Promise<number[]> {
    try {
        // h5dump -d datasetPath -s 0 -c 5 -y -w 0 filePath
        const stdout = await runCommand('h5dump', ['-d', datasetPath, '-s', '0', '-c', '5', '-y', '-w', '0', filePath]);
        const dataBlock = stdout.match(/DATA \{([\s\S]*?)\}/);
        if (dataBlock) {
            const inner = dataBlock[1];
            // Split by comma or whitespace
            const items = inner.split(/,?\s+/).filter(s => s.trim() !== '');
            const values: number[] = [];
            for (const item of items) {
                // Remove quotes if string (shouldn't happen for matrix data)
                const cleanItem = item.replace(/"/g, '');
                const val = parseFloat(cleanItem);
                if (!isNaN(val)) values.push(val);
                if (values.length >= 5) break;
            }
            return values;
        }
    } catch (e) {}
    return [];
}

/**
 * Check matrix type (integer or float)
 * For sparse matrices (isDense=false), checks actual values if storage is float.
 * For dense matrices (isDense=true), relies on storage type to save time.
 */
async function checkMatrixType(filePath: string, datasetPath: string, isDense: boolean): Promise<'integer' | 'float'> {
    try {
        const info = await getDatasetInfo(filePath, datasetPath);
        if (info.typeClass === 'integer') return 'integer';
        
        if (isDense) {
            // If dense, just say it is dense (float) and skip value check
            return 'float';
        }

        // If sparse and float storage, check first 5 non-zero values
        const values = await readFirstFiveValues(filePath, datasetPath);
        if (values.length === 0) return 'float'; // Default
        
        const allIntegers = values.every(v => Number.isInteger(v));
        return allIntegers ? 'integer' : 'float';
    } catch (e) {
        return 'float';
    }
}

/**
 * Analyze matrices (.X and layers) for data types
 */
export async function analyzeMatrices(filePath: string): Promise<MatrixInfo[]> {
    const { path: localPath, isTemp } = await ensureLocalFile(filePath);
    try {
        const matrices: MatrixInfo[] = [];
        
        // Check X
        try {
            // Check if X is group (sparse) or dataset (dense)
            // h5ls localPath/X
            const stdout = await runCommand('h5ls', [localPath + '/X']);
            
            let isSparse = false;
            // Sparse matrix in AnnData is a Group containing data, indices, indptr
            if (stdout.includes('data') && stdout.includes('indices') && stdout.includes('indptr')) {
                isSparse = true;
            }
            
            if (isSparse) {
                const type = await checkMatrixType(localPath, '/X/data', false);
                matrices.push({ name: 'X', type });
            } else {
                // Dense
                const type = await checkMatrixType(localPath, '/X', true);
                matrices.push({ name: 'X', type });
            }
        } catch (e) {}

        // Check layers
        try {
            const stdout = await runCommand('h5ls', [localPath + '/layers']);
            const lines = stdout.split('\n');
            for (const line of lines) {
                const parts = line.trim().split(/\s+/);
                if (parts.length >= 1) {
                    const name = parts[0];
                    if (name) {
                        try {
                            // Check if layer is sparse or dense
                            const layerPath = `/layers/${name}`;
                            // We need to check if it is a group or dataset.
                            // h5ls localPath/layers/name
                            let isLayerSparse = false;
                            try {
                                const layerStdout = await runCommand('h5ls', [localPath + layerPath]);
                                if (layerStdout.includes('data') && layerStdout.includes('indices') && layerStdout.includes('indptr')) {
                                    isLayerSparse = true;
                                }
                            } catch (e) {}

                            if (isLayerSparse) {
                                const type = await checkMatrixType(localPath, `${layerPath}/data`, false);
                                matrices.push({ name: `layers/${name}`, type });
                            } else {
                                const type = await checkMatrixType(localPath, layerPath, true);
                                matrices.push({ name: `layers/${name}`, type });
                            }
                        } catch (e) {}
                    }
                }
            }
        } catch (e) {}

        return matrices;
    } finally {
        if (isTemp) cleanupTempFile(localPath);
    }
}

/**
 * Analyze obsm entries
 */
export async function analyzeObsm(filePath: string): Promise<ObsmInfo[]> {
    const { path: localPath, isTemp } = await ensureLocalFile(filePath);
    try {
        const obsmInfo: ObsmInfo[] = [];
        try {
            const stdout = await runCommand('h5ls', [localPath + '/obsm']);
            const lines = stdout.split('\n');
            for (const line of lines) {
                const parts = line.trim().split(/\s+/);
                if (parts.length >= 1) {
                    const name = parts[0];
                    if (name) {
                        try {
                            const info = await getDatasetInfo(localPath, `/obsm/${name}`);
                            const columns = info.shape.length > 1 ? info.shape[1] : 1;
                            obsmInfo.push({ name, columns });
                        } catch (e) {}
                    }
                }
            }
        } catch (e) {}
        return obsmInfo;
    } finally {
        if (isTemp) cleanupTempFile(localPath);
    }
}

/**
 * Read complete h5ad metadata for design analysis
 */
export async function readH5ADMetadata(filePath: string): Promise<H5ADMetadata> {
    // We can reuse the local file path to avoid multiple downloads/copies
    const { path: localPath, isTemp } = await ensureLocalFile(filePath);
    
    try {
        const totalCells = await getTotalCells(localPath);
        const species = await detectSpecies(localPath);
        const geneCount = await getGeneCount(localPath);
        const factorNames = await listFactors(localPath);
        const matrices = await analyzeMatrices(localPath);
        const obsm = await analyzeObsm(localPath);

        const factors = new Map<string, FactorInfo>();

        for (const name of factorNames) {
            const categories = await extractCategories(localPath, name);
            const cellCounts = await countCells(localPath, name);
            const counts = categories.map(cat => cellCounts.get(cat) || 0);

            let type: FactorInfo['type'] = 'experimental';
            const nameLower = name.toLowerCase();

            if (nameLower.includes('cell_type') || nameLower.includes('celltype') ||
                nameLower.includes('cluster') || nameLower.includes('annotation')) {
                type = 'classification';
            } else if (nameLower.includes('sample') || nameLower.includes('replicate') ||
                nameLower.includes('donor') || nameLower.includes('patient')) {
                type = 'replicate';
            } else if (nameLower.includes('batch') || nameLower.includes('library')) {
                type = 'batch';
            }

            factors.set(name, { name, categories, counts, type });
        }

        return { totalCells, factors, species, geneCount, matrices, obsm };
    } finally {
        if (isTemp) cleanupTempFile(localPath);
    }
}

