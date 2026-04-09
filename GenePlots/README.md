# Gene Analysis Pipeline - Python Suite

This Python analysis suite provides comprehensive visualization and statistical analysis of gene data from multiple databases.

## 🚀 Quick Start

### Prerequisites
```bash
# Install required Python packages
pip install -r requirements.txt
```

### Running the Analysis
```bash
# Basic usage
python main_analysis.py /path/to/GeneDownloaderScriptsR/AllPackagesGenes/* phenotype_name

# Example for migraine
python main_analysis.py /data/ascher02/uqmmune1/ANNOVAR/GeneDownloaderScriptsR/AllPackagesGenes/* migraine
```

## 📊 Analysis Modules

### 1. Source Comparison (`Analysis1_SourceComparison.py`)
- **Purpose**: Compare gene counts and associations across databases
- **Outputs**: 
  - Multi-panel visualization with bar charts and heatmaps
  - Database comparison report (CSV)
  - Statistical summaries

### 2. Gene Frequency Analysis (`Analysis2_GeneFrequency.py`)
- **Purpose**: Analyze how often genes appear across databases
- **Outputs**:
  - Frequency distribution plots
  - Top genes identification
  - High-confidence gene lists

### 3. Database Overlap Analysis (`Analysis3_DatabaseOverlap.py`)
- **Purpose**: Analyze overlap between databases using Jaccard similarity
- **Outputs**:
  - Jaccard similarity matrix
  - Pairwise overlap visualizations
  - Clustering analysis

### 4. Gene Set Enrichment (`Analysis4_GeneSetEnrichment.py`)
- **Purpose**: Functional enrichment analysis with pathway information
- **Outputs**:
  - Pathway enrichment plots
  - Volcano plots for fold enrichment
  - Gene-pathway mapping

### 5. Network Analysis (`Analysis5_NetworkAnalysis.py`)
- **Purpose**: Protein-protein interaction network analysis
- **Outputs**:
  - Network topology visualization
  - Hub gene identification
  - Community detection analysis

### 6. Statistical Summary (`Analysis6_StatisticalSummary.py`)
- **Purpose**: Comprehensive statistical overview and final reports
- **Outputs**:
  - Multi-panel statistical dashboard
  - Executive summary report
  - Database quality scores
  - Gene confidence analysis

## 📁 Output Structure

After running the analysis, you'll find:

```
GenePlots/
├── main_analysis.py           # Main pipeline coordinator
├── Analysis1_SourceComparison.py
├── Analysis2_GeneFrequency.py
├── Analysis3_DatabaseOverlap.py
├── Analysis4_GeneSetEnrichment.py
├── Analysis5_NetworkAnalysis.py
├── Analysis6_StatisticalSummary.py
├── requirements.txt
├── README.md
├── plots/                     # All visualization outputs
│   ├── [phenotype]_Analysis1_SourceComparison.png
│   ├── [phenotype]_Analysis2_GeneFrequency.png
│   ├── [phenotype]_Analysis3_DatabaseOverlap.png
│   ├── [phenotype]_Analysis4_GeneSetEnrichment.png
│   ├── [phenotype]_Analysis5_NetworkAnalysis.png
│   ├── [phenotype]_Analysis6_StatisticalSummary.png
│   └── *.pdf versions
├── reports/                   # CSV reports and summaries
│   ├── [phenotype]_Analysis1_DatabaseComparison.csv
│   ├── [phenotype]_Analysis2_GeneFrequency.csv
│   ├── [phenotype]_Analysis3_JaccardMatrix.csv
│   ├── [phenotype]_Analysis4_EnrichmentResults.csv
│   ├── [phenotype]_Analysis5_NetworkMetrics.csv
│   ├── [phenotype]_Analysis6_DetailedStatistics.csv
│   ├── [phenotype]_Analysis6_GeneConfidenceAnalysis.csv
│   ├── [phenotype]_Analysis6_FinalSummary.csv
│   └── [phenotype]_ExecutiveSummary.txt
└── data/                      # Processed data files
    ├── [phenotype]_all_genes.csv
    └── [phenotype]_summary_stats.csv
```

## 🔧 Advanced Usage

### Running Individual Analyses
Each analysis module can be imported and run separately:

```python
from Analysis1_SourceComparison import run_analysis
from pathlib import Path
import pandas as pd

# Load your data
gene_data = {}
for file_path in Path("AllPackagesGenes").glob("*.csv"):
    gene_data[file_path.stem] = pd.read_csv(file_path)

# Run specific analysis
run_analysis(gene_data, "migraine")
```

### Customizing Analyses
- Modify plot parameters in each analysis module
- Adjust statistical thresholds
- Add custom databases or filtering

## 📋 Data Format Requirements

The analysis expects CSV files with gene information. Common gene column names supported:
- `Gene`
- `GeneSymbol` 
- `Gene_Symbol`
- `gene_symbol`
- `symbol`

## 🎯 Key Features

- **Comprehensive Analysis**: 6 different analytical perspectives
- **High-Quality Visualizations**: Publication-ready plots in PNG and PDF
- **Detailed Reports**: CSV files with all statistical results
- **Executive Summaries**: Human-readable findings and recommendations
- **Modular Design**: Each analysis can be run independently
- **Quality Assessment**: Database quality scoring and confidence levels
- **Network Analysis**: Protein interaction network insights

## 🔍 Interpreting Results

### Gene Confidence Levels
- **High Confidence**: Genes found in 3+ databases
- **Medium Confidence**: Genes found in 2 databases  
- **Low Confidence**: Genes found in 1 database only

### Database Quality Scores
Calculated based on:
- Gene count (40% weight)
- Coverage percentage (30% weight)
- Productivity (records per gene, 30% weight)

### Network Metrics
- **Degree Centrality**: Number of connections
- **Betweenness Centrality**: Bridge importance
- **Eigenvector Centrality**: Influence based on connected neighbors

## 🐛 Troubleshooting

### Common Issues
1. **Missing gene columns**: Ensure CSV files have gene identifier columns
2. **Memory issues**: For large datasets, consider processing subsets
3. **Plot rendering**: Install matplotlib backend for your system

### Getting Help
- Check the executive summary for key findings
- Review individual module outputs for detailed analysis
- Examine CSV reports for raw statistical data

## 📊 Complete Workflow Example

1. **Run R pipeline**: `Rscript download_genes.R migraine`
2. **Run Python analysis**: `python main_analysis.py AllPackagesGenes/* migraine`
3. **Review results**: Check `GenePlots/reports/migraine_ExecutiveSummary.txt`
4. **Examine plots**: Open PNG files in `GenePlots/plots/`
5. **Deep dive**: Analyze CSV reports in `GenePlots/reports/`

This provides a complete gene analysis pipeline from data collection to publication-ready visualizations!
