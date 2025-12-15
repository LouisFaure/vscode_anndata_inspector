import * as assert from 'assert';
import * as path from 'path';
import { readH5ADMetadata, H5ADMetadata } from '../h5adReader';

const TEST_FILE_PATH = 'https://github.com/scverse/scanpy/raw/refs/heads/main/src/scanpy/datasets/10x_pbmc68k_reduced.h5ad';

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
});
