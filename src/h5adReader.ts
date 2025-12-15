/**
 * H5AD File Reader using h5wasm
 * 
 * Reads AnnData h5ad files and extracts experimental design information.
 */

import type * as h5wasmType from 'h5wasm';
import * as path from 'path';
import * as fs from 'fs';

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

let h5wasmModule: typeof h5wasmType | null = null;

async function getH5Wasm(): Promise<typeof h5wasmType> {
    if (!h5wasmModule) {
        // Dynamic import workaround for CommonJS environment
        // We use 'h5wasm/node' which is the correct entry point for Node.js
        // Using eval to prevent TypeScript from transpiling import() to require()
        const mod = await (eval('import("h5wasm/node")') as Promise<typeof h5wasmType>);
        h5wasmModule = mod;
    }
    await h5wasmModule.ready;
    return h5wasmModule;
}

/**
 * Helper to fetch file content from URL
 */
async function fetchFileBuffer(url: string): Promise<Buffer> {
    const response = await fetch(url);
    if (!response.ok) {
        throw new Error(`Failed to fetch ${url}: ${response.statusText}`);
    }
    const arrayBuffer = await response.arrayBuffer();
    return Buffer.from(arrayBuffer);
}

/**
 * Helper to open an H5AD file from disk or URL using h5wasm VFS
 */
async function openH5File(filePath: string): Promise<h5wasmType.File> {
    const h5 = await getH5Wasm();
    
    let fileBuffer: Buffer;
    let fileName: string;
    
    // Check if it's a URL or local file path
    if (filePath.startsWith('http://') || filePath.startsWith('https://')) {
        fileBuffer = await fetchFileBuffer(filePath);
        // Extract filename from URL
        fileName = filePath.split('/').pop() || 'h5ad_file.h5ad';
    } else {
        fileBuffer = fs.readFileSync(filePath);
        fileName = path.basename(filePath);
    }

    // Write to virtual file system
    // FS is available on the module after ready
    if (h5.FS) {
        // Use /tmp to avoid read-only file system errors
        try { h5.FS.mkdir('/tmp'); } catch (e) {}
        const vfsPath = `/tmp/${fileName}`;
        h5.FS.writeFile(vfsPath, new Uint8Array(fileBuffer));
        return new h5.File(vfsPath, 'r');
    }

    return new h5.File(fileName, 'r');
}

/**
 * List all categorical factors in the /obs section of an h5ad file
 */
export async function listFactors(filePath: string): Promise<string[]> {
    const file = await openH5File(filePath);
    const h5 = await getH5Wasm();
    const factors: string[] = [];

    try {
        const obs = file.get('obs');
        if (obs && obs instanceof h5.Group) {
            const keys = obs.keys();

            for (const key of keys) {
                const item = obs.get(key);
                // Check if it's a Group (categorical in h5ad format)
                if (item && item instanceof h5.Group) {
                    // Verify it has 'categories' and 'codes' datasets (categorical structure)
                    const groupKeys = item.keys();
                    if (groupKeys.includes('categories') && groupKeys.includes('codes')) {
                        factors.push(key);
                    }
                }
            }
        }
    } finally {
        file.close();
        // Clean up VFS to save memory
        try { if (h5.FS) h5.FS.unlink(`/tmp/${path.basename(filePath)}`); } catch (e) { }
    }

    return factors.sort();
}

/**
 * Extract category names from a categorical factor
 */
export async function extractCategories(filePath: string, factorName: string): Promise<string[]> {
    const file = await openH5File(filePath);
    const h5 = await getH5Wasm();
    const categories: string[] = [];

    try {
        const categoriesPath = `obs/${factorName}/categories`;
        const dataset = file.get(categoriesPath);

        if (dataset && dataset instanceof h5.Dataset) {
            const values = dataset.value;

            if (Array.isArray(values) || ArrayBuffer.isView(values)) {
                // Handle different array types
                const arr = values as any[];
                for (const val of arr) {
                    categories.push(String(val));
                }
            }
        }
    } finally {
        file.close();
        try { if (h5.FS) h5.FS.unlink(`/tmp/${path.basename(filePath)}`); } catch (e) { }
    }

    return categories;
}

/**
 * Count cells per category for a factor
 */
export async function countCells(filePath: string, factorName: string): Promise<Map<string, number>> {
    const file = await openH5File(filePath);
    const h5 = await getH5Wasm();
    const counts = new Map<string, number>();

    try {
        // We need categories to map codes
        // Getting them from file again to avoid re-opening file multiple times if called separately
        // But here we are already inside an open file context, so we read directly

        let categories: string[] = [];
        const catDataset = file.get(`obs/${factorName}/categories`);
        if (catDataset && catDataset instanceof h5.Dataset) {
            const values = catDataset.value;
            if (Array.isArray(values) || ArrayBuffer.isView(values)) {
                const arr = values as any[];
                categories = arr.map(v => String(v));
            }
        }

        // Get codes
        const codesPath = `obs/${factorName}/codes`;
        const codesDataset = file.get(codesPath);

        if (codesDataset && codesDataset instanceof h5.Dataset) {
            const codes = codesDataset.value as any; // Typed array

            // Initialize counts
            for (const cat of categories) {
                counts.set(cat, 0);
            }

            // Count occurrences
            // Iterate over typed array
            if (codes.length) {
                for (let i = 0; i < codes.length; i++) {
                    const code = codes[i];
                    if (code >= 0 && code < categories.length) {
                        const cat = categories[code];
                        counts.set(cat, (counts.get(cat) || 0) + 1);
                    }
                }
            }
        }
    } finally {
        file.close();
        try { if (h5.FS) h5.FS.unlink(`/tmp/${path.basename(filePath)}`); } catch (e) { }
    }

    return counts;
}

/**
 * Get total cell count from the h5ad file
 */
export async function getTotalCells(filePath: string): Promise<number> {
    const file = await openH5File(filePath);
    const h5 = await getH5Wasm();

    try {
        // Check obs/_index for cell count
        const indexDataset = file.get('obs/_index');
        if (indexDataset && indexDataset instanceof h5.Dataset) {
            return indexDataset.shape ? indexDataset.shape[0] : 0;
        }

        // Fallback: check X matrix shape
        const xDataset = file.get('X');
        if (xDataset && xDataset instanceof h5.Dataset) {
            return xDataset.shape ? xDataset.shape[0] : 0;
        }
    } finally {
        file.close();
        try { if (h5.FS) h5.FS.unlink(`/tmp/${path.basename(filePath)}`); } catch (e) { }
    }

    return 0;
}

/**
 * Detect species from gene naming patterns in /var
 */
export async function detectSpecies(filePath: string): Promise<string> {
    const file = await openH5File(filePath);
    const h5 = await getH5Wasm();

    try {
        // Try to get gene names from var/_index or var/gene_symbols etc.
        const possiblePaths = ['var/_index', 'var/index', 'var/gene_symbols', 'var/gene_names'];
        let geneNames: string[] = [];

        for (const varPath of possiblePaths) {
            const dataset = file.get(varPath);
            if (dataset && dataset instanceof h5.Dataset) {
                const values = dataset.value;
                if ((Array.isArray(values) || ArrayBuffer.isView(values)) && values.length > 0) {
                    // Take first 50
                    const arr = values as any[];
                    geneNames = arr.slice(0, 50).map(v => String(v));
                    break;
                }
            }
        }

        if (geneNames.length === 0) {
            return 'unknown';
        }

        // Analyze naming patterns
        let uppercase = 0;  // Human: ACTB, GAPDH
        let titlecase = 0;  // Mouse: Actb, Gapdh
        let lowercase = 0;  // Zebrafish: actb, gapdh

        for (const gene of geneNames) {
            if (/^[A-Z][A-Z0-9]+$/.test(gene)) {
                uppercase++;
            } else if (/^[A-Z][a-z0-9]+$/.test(gene)) {
                titlecase++;
            } else if (/^[a-z][a-z0-9]+$/.test(gene)) {
                lowercase++;
            }
        }

        const total = uppercase + titlecase + lowercase;
        if (total === 0) {
            return 'unknown';
        }

        if (uppercase / total > 0.5) {
            return 'human';
        } else if (titlecase / total > 0.5) {
            return 'mouse';
        } else if (lowercase / total > 0.5) {
            return 'zebrafish';
        }

        return 'unknown';
    } finally {
        file.close();
        try { if (h5.FS) h5.FS.unlink(`/tmp/${path.basename(filePath)}`); } catch (e) { }
    }
}

/**
 * Get gene count from the h5ad file
 */
export async function getGeneCount(filePath: string): Promise<number> {
    const file = await openH5File(filePath);
    const h5 = await getH5Wasm();

    try {
        // Check var/_index or var/index
        const possiblePaths = ['var/_index', 'var/index'];
        for (const p of possiblePaths) {
            const dataset = file.get(p);
            if (dataset && dataset instanceof h5.Dataset) {
                return dataset.shape ? dataset.shape[0] : 0;
            }
        }

        // Check any dataset in var
        const varGroup = file.get('var');
        if (varGroup && varGroup instanceof h5.Group) {
             const keys = varGroup.keys();
             if (keys.length > 0) {
                 const first = varGroup.get(keys[0]);
                 if (first && first instanceof h5.Dataset) {
                     return first.shape ? first.shape[0] : 0;
                 }
             }
        }

        // Fallback: check X matrix shape (n_obs, n_vars)
        const xDataset = file.get('X');
        if (xDataset && xDataset instanceof h5.Dataset) {
            return xDataset.shape && xDataset.shape.length > 1 ? xDataset.shape[1] : 0;
        }
    } finally {
        file.close();
        try { if (h5.FS) h5.FS.unlink(`/tmp/${path.basename(filePath)}`); } catch (e) { }
    }

    return 0;
}

/**
 * Analyze matrices (.X and layers) for data types
 */
export async function analyzeMatrices(filePath: string): Promise<MatrixInfo[]> {
    const file = await openH5File(filePath);
    const h5 = await getH5Wasm();
    const matrices: MatrixInfo[] = [];

    try {
        // Check .X
        const xType = await getMatrixType(file, 'X', h5);
        matrices.push({ name: 'X', type: xType });

        // Check layers
        const layers = file.get('layers');
        if (layers && layers instanceof h5.Group) {
            for (const key of layers.keys()) {
                const type = await getMatrixType(file, `layers/${key}`, h5);
                matrices.push({ name: `layers/${key}`, type });
            }
        }
    } finally {
        file.close();
        try { if (h5.FS) h5.FS.unlink(`/tmp/${path.basename(filePath)}`); } catch (e) { }
    }
    return matrices;
}

async function getMatrixType(file: h5wasmType.File, path: string, h5: typeof h5wasmType): Promise<'integer' | 'float'> {
    const item = file.get(path);
    if (!item) return 'float';

    let dataDataset: h5wasmType.Dataset | null = null;

    if (item instanceof h5.Dataset) {
        dataDataset = item;
    } else if (item instanceof h5.Group) {
        const dataItem = item.get('data');
        if (dataItem && dataItem instanceof h5.Dataset) {
            dataDataset = dataItem;
        }
    }

    if (dataDataset && dataDataset.shape) {
        // Check dtype
        // If it's explicitly integer type
        if (String(dataDataset.dtype).includes('int')) return 'integer';

        // Read first few values to check if they are effectively integers
        try {
            let selection: any[] = [];
            if (dataDataset.shape.length === 1) {
                selection = [[0, Math.min(5, dataDataset.shape[0])]];
            } else {
                selection = dataDataset.shape.map((dim, i) => {
                    if (i === 0) return [0, 1]; // 1 row
                    if (i === 1) return [0, Math.min(5, dim)]; // 5 cols
                    return [0, 1];
                });
            }
            
            const values = dataDataset.slice(selection);
            if (values) {
                let arr: any[] = [];
                if (Array.isArray(values)) {
                    arr = (values as any[]).flat(Infinity);
                } else {
                    arr = Array.from(values as any);
                }
                
                // Check if all are integers
                const allInt = arr.every((n: any) => typeof n === 'number' && n % 1 === 0);
                return allInt ? 'integer' : 'float';
            }
        } catch (e) {
            // Fallback if slice fails
            return 'float';
        }
    }

    return 'float';
}

/**
 * Analyze obsm entries
 */
export async function analyzeObsm(filePath: string): Promise<ObsmInfo[]> {
    const file = await openH5File(filePath);
    const h5 = await getH5Wasm();
    const obsmInfo: ObsmInfo[] = [];

    try {
        const obsm = file.get('obsm');
        if (obsm && obsm instanceof h5.Group) {
            for (const key of obsm.keys()) {
                const item = obsm.get(key);
                if (item && item instanceof h5.Dataset && item.shape) {
                    // shape is (n_obs, n_features)
                    const columns = item.shape.length > 1 ? item.shape[1] : 1;
                    obsmInfo.push({ name: key, columns });
                }
            }
        }
    } finally {
        file.close();
        try { if (h5.FS) h5.FS.unlink(`/tmp/${path.basename(filePath)}`); } catch (e) { }
    }
    return obsmInfo;
}

/**
 * Read complete h5ad metadata for design analysis
 */
export async function readH5ADMetadata(filePath: string): Promise<H5ADMetadata> {
    // We get individual pieces. Optimization: open file once and read everything.
    // However, to ensure consistency and speed, let's implement a single pass version
    // or just call them sequentially. calling sequentially is fine for now.

    const totalCells = await getTotalCells(filePath);
    const species = await detectSpecies(filePath);
    const geneCount = await getGeneCount(filePath);
    const factorNames = await listFactors(filePath);
    const matrices = await analyzeMatrices(filePath);
    const obsm = await analyzeObsm(filePath);

    const factors = new Map<string, FactorInfo>();

    for (const name of factorNames) {
        const categories = await extractCategories(filePath, name);
        const cellCounts = await countCells(filePath, name);
        const counts = categories.map(cat => cellCounts.get(cat) || 0);

        // Determine factor type based on name heuristics
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
}
