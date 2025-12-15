/**
 * AnnData Design Inspector - VS Code Extension
 * 
 * Main extension entry point.
 */

import * as vscode from 'vscode';
import * as path from 'path';
import { H5ADEditorProvider } from './h5adEditorProvider';
import { readH5ADMetadata, FactorInfo } from './h5adReader';

export function activate(context: vscode.ExtensionContext) {
    console.log('AnnData Design Inspector is now active');

    // Register custom editor for h5ad files
    context.subscriptions.push(H5ADEditorProvider.register(context));

    // Register inspect file command
    context.subscriptions.push(
        vscode.commands.registerCommand('anndata.inspectFile', async () => {
            const fileUris = await vscode.window.showOpenDialog({
                canSelectMany: false,
                filters: {
                    'H5AD Files': ['h5ad'],
                    'All Files': ['*']
                },
                title: 'Select H5AD File to Inspect'
            });

            if (fileUris && fileUris.length > 0) {
                const uri = fileUris[0];
                await vscode.commands.executeCommand('vscode.openWith', uri, H5ADEditorProvider.viewType);
            }
        })
    );
}



export function deactivate() {
    console.log('AnnData Design Inspector is now deactivated');
}
