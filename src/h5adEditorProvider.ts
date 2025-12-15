/**
 * H5AD Custom Editor Provider
 * 
 * Provides a custom editor for viewing h5ad files in VS Code.
 * Displays experimental design information in a rich webview.
 */

import * as vscode from 'vscode';
import * as path from 'path';
import { readH5ADMetadata, FactorInfo } from './h5adReader';

export class H5ADEditorProvider implements vscode.CustomReadonlyEditorProvider {
    public static readonly viewType = 'anndata.h5adViewer';

    constructor(private readonly context: vscode.ExtensionContext) { }

    public static register(context: vscode.ExtensionContext): vscode.Disposable {
        const provider = new H5ADEditorProvider(context);
        return vscode.window.registerCustomEditorProvider(
            H5ADEditorProvider.viewType,
            provider,
            {
                webviewOptions: { retainContextWhenHidden: true },
                supportsMultipleEditorsPerDocument: false
            }
        );
    }

    async openCustomDocument(
        uri: vscode.Uri,
        _openContext: vscode.CustomDocumentOpenContext,
        _token: vscode.CancellationToken
    ): Promise<vscode.CustomDocument> {
        return { uri, dispose: () => { } };
    }

    async resolveCustomEditor(
        document: vscode.CustomDocument,
        webviewPanel: vscode.WebviewPanel,
        _token: vscode.CancellationToken
    ): Promise<void> {
        webviewPanel.webview.options = {
            enableScripts: true
        };

        // Show loading state
        webviewPanel.webview.html = this.getLoadingHtml();

        try {
            // Read and analyze the h5ad file
            const filePath = document.uri.fsPath;
            const metadata = await readH5ADMetadata(filePath);

            // Generate webview content
            webviewPanel.webview.html = this.getWebviewContent(
                webviewPanel.webview,
                document.uri.fsPath,
                metadata
            );

            // Handle messages from webview
            webviewPanel.webview.onDidReceiveMessage(
                async message => {
                    // No messages handled currently
                },
                undefined,
                this.context.subscriptions
            );
        } catch (error) {
            webviewPanel.webview.html = this.getErrorHtml(error);
        }
    }



    private getLoadingHtml(): string {
        return `<!DOCTYPE html>
<html>
<head>
    <style>
        body {
            display: flex;
            justify-content: center;
            align-items: center;
            height: 100vh;
            margin: 0;
            font-family: var(--vscode-font-family);
            color: var(--vscode-foreground);
            background: var(--vscode-editor-background);
        }
        .loader {
            text-align: center;
        }
        .spinner {
            width: 40px;
            height: 40px;
            border: 3px solid var(--vscode-foreground);
            border-top: 3px solid transparent;
            border-radius: 50%;
            animation: spin 1s linear infinite;
            margin: 0 auto 16px;
        }
        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }
    </style>
</head>
<body>
    <div class="loader">
        <div class="spinner"></div>
        <p>Loading H5AD file...</p>
    </div>
</body>
</html>`;
    }

    private getErrorHtml(error: any): string {
        return `<!DOCTYPE html>
<html>
<head>
    <style>
        body {
            padding: 24px;
            font-family: var(--vscode-font-family);
            color: var(--vscode-foreground);
            background: var(--vscode-editor-background);
        }
        .error {
            color: var(--vscode-errorForeground);
            background: var(--vscode-inputValidation-errorBackground);
            border: 1px solid var(--vscode-inputValidation-errorBorder);
            padding: 16px;
            border-radius: 4px;
        }
        pre {
            white-space: pre-wrap;
            word-wrap: break-word;
        }
    </style>
</head>
<body>
    <h1>Error Loading H5AD File</h1>
    <div class="error">
        <p><strong>Failed to read the H5AD file.</strong></p>
        <pre>${error instanceof Error ? error.message : String(error)}</pre>
    </div>
    <p>Please ensure:</p>
    <ul>
        <li>The file is a valid AnnData H5AD file</li>
        <li>The file is not corrupted</li>
        <li>You have read permissions for the file</li>
    </ul>
</body>
</html>`;
    }

    private getWebviewContent(
        webview: vscode.Webview,
        filePath: string,
        metadata: any
    ): string {
        const fileName = path.basename(filePath);
        const fileSize = this.formatFileSize(filePath);

        const factorsRecord = this.mapFactorsToRecord(metadata.factors);

        // Build matrices rows
        const matricesRows = (metadata.matrices || []).map((m: any) => {
            return `<tr>
                <td><code>${m.name}</code></td>
                <td><span class="badge ${m.type}">${this.capitalize(m.type)}</span></td>
            </tr>`;
        }).join('');

        // Build obsm rows
        const obsmRows = (metadata.obsm || []).map((o: any) => {
            return `<tr>
                <td><code>${o.name}</code></td>
                <td>${o.columns} columns</td>
            </tr>`;
        }).join('');

        // Build factors table rows
        const factorsRows = Object.entries(factorsRecord)
            .map(([name, info]: [string, any]) => {
                return `<tr>
                    <td><strong>${this.toTitleCase(name)}</strong></td>
                    <td>${info.categories.length}</td>
                    <td>${info.categories.slice(0, 5).join(', ')}${info.categories.length > 5 ? '...' : ''}</td>
                </tr>`;
            })
            .join('');


        return `<!DOCTYPE html>
<html>
<head>
    <style>
        :root {
            --bg-color: var(--vscode-editor-background);
            --fg-color: var(--vscode-foreground);
            --border-color: var(--vscode-panel-border);
            --accent-color: var(--vscode-textLink-foreground);
            --muted-color: var(--vscode-descriptionForeground);
        }
        body {
            padding: 24px;
            font-family: var(--vscode-font-family);
            color: var(--fg-color);
            background: var(--bg-color);
            line-height: 1.6;
            max-width: 1000px;
            margin: 0 auto;
        }
        h1 {
            border-bottom: 2px solid var(--border-color);
            padding-bottom: 12px;
            margin-bottom: 24px;
        }
        h2 {
            color: var(--accent-color);
            margin-top: 32px;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 8px;
        }
        h3 {
            font-size: 1.1em;
            margin-top: 0;
            margin-bottom: 16px;
        }
        .header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 16px;
        }
        .meta-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px;
            margin: 24px 0;
        }
        .meta-card {
            background: var(--vscode-input-background);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 16px;
        }
        .meta-label {
            font-size: 0.85em;
            color: var(--muted-color);
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        .meta-value {
            font-size: 1.5em;
            font-weight: bold;
            margin-top: 4px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 16px 0;
        }
        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }
        th {
            background: var(--vscode-input-background);
            font-weight: 600;
        }
        code {
            font-family: var(--vscode-editor-font-family);
            background: var(--vscode-textCodeBlock-background);
            padding: 2px 4px;
            border-radius: 4px;
        }
        .badge {
            display: inline-block;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 0.85em;
            font-weight: 500;
        }
        .badge.integer {
            background-color: var(--vscode-charts-blue);
            color: white;
        }
        .badge.float {
            background-color: var(--vscode-charts-orange);
            color: white;
        }
        .grammar-box {
            background: var(--vscode-textCodeBlock-background);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 16px 24px;
            font-family: var(--vscode-editor-font-family);
            font-size: 1.2em;
            margin: 16px 0;
        }
        .design-diagram {
            background: var(--vscode-input-background);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 24px;
            margin: 16px 0;
            text-align: center;
        }
        .diagram-node {
            display: inline-block;
            background: var(--vscode-button-background);
            color: var(--vscode-button-foreground);
            padding: 8px 16px;
            border-radius: 4px;
            margin: 4px;
        }
        .diagram-arrow {
            color: var(--muted-color);
            font-size: 1.5em;
            margin: 0 8px;
        }
        .diagram-row {
            margin: 8px 0;
        }
        .context-box {
            background: var(--vscode-input-background);
            border-left: 4px solid var(--accent-color);
            padding: 16px;
            margin: 16px 0;
        }
        button {
            background: var(--vscode-button-background);
            color: var(--vscode-button-foreground);
            border: none;
            padding: 10px 20px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
        }
        button:hover {
            background: var(--vscode-button-hoverBackground);
        }
        .section-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 24px;
            margin-top: 24px;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>📊 ${fileName}</h1>
    </div>
    
    <div class="meta-grid">
        <div class="meta-card">
            <div class="meta-label">Total Cells</div>
            <div class="meta-value">${metadata.totalCells.toLocaleString()}</div>
        </div>
        <div class="meta-card">
            <div class="meta-label">Genes</div>
            <div class="meta-value">${metadata.geneCount.toLocaleString()}</div>
        </div>
        <div class="meta-card">
            <div class="meta-label">Species</div>
            <div class="meta-value">${this.getSpeciesEmoji(metadata.species)} ${this.capitalize(metadata.species)}</div>
        </div>
    </div>

    <h2>Data Structure</h2>
    <div class="section-grid">
        <div>
            <h3>Data Matrices</h3>
            <table>
                <thead>
                    <tr>
                        <th>Name</th>
                        <th>Type</th>
                    </tr>
                </thead>
                <tbody>
                    ${matricesRows}
                </tbody>
            </table>
        </div>
        
        <div>
            <h3>Embeddings (obsm)</h3>
            <table>
                <thead>
                    <tr>
                        <th>Name</th>
                        <th>Dimensions</th>
                    </tr>
                </thead>
                <tbody>
                    ${obsmRows}
                </tbody>
            </table>
        </div>
    </div>
    
    <h2>Experimental Factors</h2>
    <table>
        <thead>
            <tr>
                <th>Factor</th>
                <th>Levels</th>
                <th>Categories</th>
            </tr>
        </thead>
        <tbody>
            ${factorsRows}
        </tbody>
    </table>
    
    <script>
        const vscode = acquireVsCodeApi();
    </script>
</body>
</html>`;
    }

    private formatFileSize(filePath: string): string {
        try {
            const fs = require('fs');
            const stats = fs.statSync(filePath);
            const bytes = stats.size;

            if (bytes < 1024) return bytes + ' B';
            if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + ' KB';
            if (bytes < 1024 * 1024 * 1024) return (bytes / (1024 * 1024)).toFixed(1) + ' MB';
            return (bytes / (1024 * 1024 * 1024)).toFixed(1) + ' GB';
        } catch {
            return 'Unknown';
        }
    }

    private mapFactorsToRecord(factors: Map<string, FactorInfo>): Record<string, FactorInfo> {
        const record: Record<string, FactorInfo> = {};
        for (const [name, info] of factors) {
            record[name] = info;
        }
        return record;
    }

    private toTitleCase(name: string): string {
        return name.replace(/_/g, ' ').split(' ').map(w => w.charAt(0).toUpperCase() + w.slice(1)).join(' ');
    }

    private capitalize(str: string): string {
        return str.charAt(0).toUpperCase() + str.slice(1);
    }

    private getSpeciesEmoji(species: string): string {
        const emojis: Record<string, string> = {
            human: '🧬',
            mouse: '🐭',
            zebrafish: '🐟',
            drosophila: '🪰',
            unknown: '❓'
        };
        return emojis[species] || '🧬';
    }
}
