# AnnData Inspector (VS Code Extension)

Inspect and visualize single cell metadata from AnnData (`.h5ad`) files directly in VS Code.

## Features

- **Metadata Summary**: Detects how many cells and genes the h5ad data contains. Automatically identifies organism (Human, Mouse, Zebrafish, Drosophila) from gene naming patterns.
- **Data Structure**: Detects the count matrices in `.X` and `.layers` and checked the presence of integer counts. Detects `.obsm` entries to show the various dimensionality reductions saved.
- **Pure TypeScript**: No Python dependencies required - uses WebAssembly for fast HDF5 reading.

## Usage

### Inspect a File
1. Open any `.h5ad` file in VS Code (it will open in the custom viewer automatically).
2. OR use the command `AnnData: Inspect H5AD File` from the Command Palette.

## Requirements

- VS Code 1.85.0 or newer

## Installation via Source

1. Clone the repository
2. Run `npm install`
3. Run `npm run compile`
4. Press F5 to debug in VS Code
