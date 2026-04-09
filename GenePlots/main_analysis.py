#!/usr/bin/env python3
"""
Main Gene Analysis Pipeline
===========================
Comprehensive analysis of genes from multiple databases for a given phenotype.

Usage:
    python main_analysis.py <data_directory> <phenotype>

Example:
    python main_analysis.py AllPackagesGenes migraine
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set up plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def setup_directories():
    """Create necessary directories for analysis outputs"""
    directories = [
        'GenePlots',
        'GenePlots/data',
        'GenePlots/plots',
        'GenePlots/reports'
    ]
    
    for dir_name in directories:
        Path(dir_name).mkdir(parents=True, exist_ok=True)
    
    print(f"📁 Created analysis directories: {', '.join(directories)}")

def load_gene_data(data_directory, phenotype):
    """
    Load all gene data files for the specified phenotype
    
    Args:
        data_directory (str): Path to directory containing gene data files
        phenotype (str): Phenotype name to analyze
    
    Returns:
        dict: Dictionary of DataFrames with database names as keys
    """
    data_path = Path(data_directory)
    if not data_path.exists():
        raise FileNotFoundError(f"Data directory not found: {data_directory}")
    
    # Define database mappings
    database_files = {
        'KEGG': f'{phenotype}_kegg_robust.csv',
        'Gene_Ontology': f'{phenotype}_go.csv',
        'ClinVar': f'{phenotype}_clinvar.csv',
        'HPO': f'{phenotype}_hpo.csv',
        'GWAS': f'{phenotype}_gwas.csv',
        'GTEx': f'{phenotype}_gtex.csv',
        'PubMed': f'{phenotype}_pubmed.csv',
        'UniProt': f'{phenotype}_uniprot.csv',
        'STRING': f'{phenotype}_string.csv',
        'Reactome': f'{phenotype}_reactome.csv',
        'OMIM': f'{phenotype}_omim.csv',
        'OpenTargets': f'{phenotype}_opentargets.csv',
        'DisGeNET': f'{phenotype}_disgenet.csv',
        'STRING_DB': f'{phenotype}_string_db.csv'
    }
    
    # Load summary files
    summary_file = data_path / f'{phenotype}_SOURCES_SUMMARY.csv'
    all_genes_file = data_path / f'{phenotype}_ALL_SOURCES_GENES.csv'
    
    gene_data = {}
    
    # Load individual database files
    for db_name, filename in database_files.items():
        file_path = data_path / filename
        if file_path.exists():
            try:
                df = pd.read_csv(file_path)
                gene_data[db_name] = df
                print(f"✅ Loaded {db_name}: {len(df)} records")
            except Exception as e:
                print(f"⚠️ Error loading {db_name}: {e}")
        else:
            print(f"❌ File not found: {filename}")
    
    # Load summary data
    if summary_file.exists():
        gene_data['SUMMARY'] = pd.read_csv(summary_file)
        print(f"✅ Loaded summary: {len(gene_data['SUMMARY'])} sources")
    
    if all_genes_file.exists():
        gene_data['ALL_GENES'] = pd.read_csv(all_genes_file)
        print(f"✅ Loaded all genes: {len(gene_data['ALL_GENES'])} unique genes")
    
    return gene_data

def run_analysis_modules(gene_data, phenotype):
    """
    Run all analysis modules
    
    Args:
        gene_data (dict): Dictionary of DataFrames
        phenotype (str): Phenotype name
    """
    print("\n🔬 RUNNING ANALYSIS MODULES")
    print("=" * 50)
    
    # Import and run each analysis module
    analysis_modules = [
        'Analysis1_SourceComparison',
        'Analysis2_GeneFrequency', 
        'Analysis3_DatabaseOverlap',
        'Analysis4_GeneSetEnrichment',
        'Analysis5_NetworkAnalysis',
        'Analysis6_StatisticalSummary'
    ]
    
    for module_name in analysis_modules:
        try:
            print(f"\n🔄 Running {module_name}...")
            
            # Import the module dynamically
            module = __import__(module_name)
            
            # Run the analysis
            if hasattr(module, 'run_analysis'):
                module.run_analysis(gene_data, phenotype)
                print(f"✅ {module_name} completed successfully")
            else:
                print(f"⚠️ {module_name} missing run_analysis function")
                
        except ImportError as e:
            print(f"❌ Could not import {module_name}: {e}")
        except Exception as e:
            print(f"❌ Error in {module_name}: {e}")

def main():
    """Main analysis pipeline"""
    parser = argparse.ArgumentParser(
        description='Comprehensive Gene Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python main_analysis.py AllPackagesGenes migraine
    python main_analysis.py /path/to/data diabetes
        """
    )
    
    parser.add_argument('data_directory', 
                       help='Directory containing gene data files')
    parser.add_argument('phenotype', 
                       help='Phenotype name to analyze')
    parser.add_argument('--output-dir', 
                       default='GenePlots',
                       help='Output directory for analysis results')
    
    args = parser.parse_args()
    
    print("🧬 GENE ANALYSIS PIPELINE")
    print("=" * 50)
    print(f"📁 Data directory: {args.data_directory}")
    print(f"🎯 Phenotype: {args.phenotype}")
    print(f"📊 Output directory: {args.output_dir}")
    
    # Setup directories
    setup_directories()
    
    # Load data
    print(f"\n📥 LOADING DATA FOR: {args.phenotype.upper()}")
    print("=" * 50)
    gene_data = load_gene_data(args.data_directory, args.phenotype)
    
    if not gene_data:
        print("❌ No data loaded. Please check the data directory and phenotype name.")
        sys.exit(1)
    
    # Run analyses
    run_analysis_modules(gene_data, args.phenotype)
    
    print(f"\n🎉 ANALYSIS COMPLETE!")
    print("=" * 50)
    print(f"📊 Results saved in: {args.output_dir}/")
    print("📈 Plots saved in: GenePlots/plots/")
    print("📄 Reports saved in: GenePlots/reports/")

if __name__ == "__main__":
    main()
