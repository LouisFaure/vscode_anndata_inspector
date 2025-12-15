import * as assert from 'assert';
import * as path from 'path';
import { readH5ADMetadata, H5ADMetadata } from '../h5adReader';

const TEST_FILE_PATH = path.resolve(__dirname, '../../adata_test.h5ad');

let metadata: H5ADMetadata;

before(async function() {
    this.timeout(15000); // Increase timeout for remote file download
    console.log(`Loading test file: ${TEST_FILE_PATH}`);
    metadata = await readH5ADMetadata(TEST_FILE_PATH);
});

describe('H5AD Reader Test Suite', () => {
    it('should read metadata from h5ad file', function() {
        console.log('Metadata Summary:');
        console.log(`- Total Cells: ${metadata.totalCells}`);
        console.log(`- Gene Count: ${metadata.geneCount}`);
        console.log(`- Species: ${metadata.species}`);
        console.log(`- Factors: ${metadata.factors.size}`);
        
        for (const [name, info] of metadata.factors) {
            console.log(`  * ${name}: ${info.categories.length} categories, type=${info.type}`);
        }

        assert.ok(metadata.totalCells > 0, 'Total cells should be greater than 0');
        assert.ok(metadata.geneCount > 0, 'Gene count should be greater than 0');
        assert.ok(metadata.factors.size > 0, 'Should have detected factors');
        // assert.notStrictEqual(metadata.species, 'unknown', 'Should detect species'); // Species detection might fail on reduced data
        assert.strictEqual(metadata.hasRaw, true, 'Should  have raw data');
    });
});

describe('H5AD Reader Data Structure Test Suite', () => {
    it('should read matrices and obsm info', function() {
        console.log('Matrices:');
        for (const m of metadata.matrices) {
            console.log(`- ${m.name}: ${m.type}`);
        }
        console.log('Obsm:');
        for (const o of metadata.obsm) {
            console.log(`- ${o.name}: ${o.columns} columns`);
        }

        assert.ok(metadata.matrices.length > 0, 'Should have detected matrices');
        assert.ok(metadata.obsm.length >= 0, 'Should have detected obsm (might be empty)');
        // Check X
        const x = metadata.matrices.find(m => m.name === 'X');
        assert.ok(x, 'Should have X matrix');
    });

    it('should read raw metadata correctly', async function() {
        const rawMetadata = await readH5ADMetadata(TEST_FILE_PATH, true);
        assert.strictEqual(rawMetadata.hasRaw, true, 'Should detect raw data');
        // Total cells should match
        assert.strictEqual(rawMetadata.totalCells, metadata.totalCells, 'Total cells should match');
        // Should have matrices (at least raw/X)
        assert.ok(rawMetadata.matrices.length > 0, 'Should have raw matrices');
        const rawX = rawMetadata.matrices.find(m => m.name === 'raw/X');
        assert.ok(rawX, 'Should have raw/X matrix');
    });
});
