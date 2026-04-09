#!/usr/bin/env python3
"""
Analysis 2: Gene Frequency Analysis
===================================
Analyzes how frequently genes appear across different databases and identifies
the most commonly found genes for the phenotype.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from collections import Counter

def run_analysis(gene_data, phenotype):
    """
    Run gene frequency analysis
    
    Args:
        gene_data (dict): Dictionary of DataFrames with gene data
        phenotype (str): Phenotype name
    """
    print("🔬 Gene Frequency Analysis")
    
    # Collect all genes and their source databases
    gene_sources = {}
    all_genes = []
    
    # Define gene column names to look for
    gene_columns = ['Gene', 'GeneSymbol', 'Gene_Symbol', 'gene_symbol', 'symbol']
    
    for db_name, df in gene_data.items():
        if db_name in ['SUMMARY', 'ALL_GENES']:
            continue
            
        # Find the gene column
        gene_col = None
        for col in gene_columns:
            if col in df.columns:
                gene_col = col
                break
        
        if gene_col:
            genes_in_db = df[gene_col].dropna().unique()
            for gene in genes_in_db:
                if gene != '' and pd.notna(gene):
                    if gene not in gene_sources:
                        gene_sources[gene] = []
                    gene_sources[gene].append(db_name)
                    all_genes.append(gene)
    
    if not gene_sources:
        print("❌ No genes found in the data")
        return
    
    # Calculate gene frequencies
    gene_frequency = {gene: len(sources) for gene, sources in gene_sources.items()}
    frequency_counter = Counter(gene_frequency.values())
    
    # Create analysis plots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'Gene Frequency Analysis for {phenotype.title()}', fontsize=16, fontweight='bold')
    
    # Plot 1: Distribution of gene frequencies
    ax1 = axes[0, 0]
    frequencies = list(frequency_counter.keys())
    counts = list(frequency_counter.values())
    
    bars1 = ax1.bar(frequencies, counts, color=sns.color_palette("viridis", len(frequencies)))
    ax1.set_title('Distribution of Gene Frequencies', fontweight='bold')
    ax1.set_xlabel('Number of Databases')
    ax1.set_ylabel('Number of Genes')
    
    # Add value labels
    for i, bar in enumerate(bars1):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 2: Top 20 most frequent genes
    ax2 = axes[0, 1]
    top_genes = sorted(gene_frequency.items(), key=lambda x: x[1], reverse=True)[:20]
    
    if top_genes:
        gene_names = [gene for gene, freq in top_genes]
        frequencies = [freq for gene, freq in top_genes]
        
        bars2 = ax2.barh(range(len(gene_names)), frequencies, 
                         color=sns.color_palette("plasma", len(gene_names)))
        ax2.set_title('Top 20 Most Frequent Genes', fontweight='bold')
        ax2.set_xlabel('Number of Databases')
        ax2.set_ylabel('Genes')
        ax2.set_yticks(range(len(gene_names)))
        ax2.set_yticklabels(gene_names)
        ax2.invert_yaxis()
        
        # Add value labels
        for i, bar in enumerate(bars2):
            width = bar.get_width()
            ax2.text(width + 0.1, bar.get_y() + bar.get_height()/2.,
                    f'{int(width)}', ha='left', va='center', fontweight='bold')
    
    # Plot 3: Gene frequency vs rank (Zipf-like distribution)
    ax3 = axes[1, 0]
    sorted_frequencies = sorted(gene_frequency.values(), reverse=True)
    ranks = range(1, len(sorted_frequencies) + 1)
    
    ax3.loglog(ranks, sorted_frequencies, 'o-', alpha=0.7, markersize=3)
    ax3.set_title('Gene Frequency Distribution (Log-Log)', fontweight='bold')
    ax3.set_xlabel('Gene Rank (log scale)')
    ax3.set_ylabel('Frequency (log scale)')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Database coverage analysis
    ax4 = axes[1, 1]
    
    # Calculate what percentage of genes appear in each number of databases
    max_freq = max(gene_frequency.values())
    coverage_data = []
    coverage_labels = []
    
    for i in range(1, max_freq + 1):
        genes_in_i_dbs = sum(1 for freq in gene_frequency.values() if freq == i)
        coverage_data.append(genes_in_i_dbs)
        coverage_labels.append(f'{i} DB{"s" if i > 1 else ""}')
    
    colors = sns.color_palette("Set3", len(coverage_data))
    wedges, texts, autotexts = ax4.pie(coverage_data, labels=coverage_labels, autopct='%1.1f%%', 
                                       colors=colors, startangle=90)
    ax4.set_title('Gene Distribution by Database Coverage', fontweight='bold')
    
    plt.tight_layout()
    
    # Save the plot
    output_path = Path('GenePlots/plots')
    output_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path / f'{phenotype}_Analysis2_GeneFrequency.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(output_path / f'{phenotype}_Analysis2_GeneFrequency.pdf', 
                bbox_inches='tight')
    plt.close()
    
    # Create detailed gene frequency report
    gene_freq_report = []
    for gene, sources in gene_sources.items():
        gene_freq_report.append({
            'Gene': gene,
            'Frequency': len(sources),
            'Databases': ', '.join(sources),
            'Coverage_Percentage': (len(sources) / len([db for db in gene_data.keys() 
                                                      if db not in ['SUMMARY', 'ALL_GENES']]) * 100)
        })
    
    gene_freq_df = pd.DataFrame(gene_freq_report)
    gene_freq_df = gene_freq_df.sort_values('Frequency', ascending=False)
    
    # Save detailed report
    report_path = Path('GenePlots/reports')
    report_path.mkdir(parents=True, exist_ok=True)
    gene_freq_df.to_csv(report_path / f'{phenotype}_Analysis2_GeneFrequency_Report.csv', index=False)
    
    # Create summary statistics
    stats_summary = {
        'Total_Unique_Genes': len(gene_sources),
        'Average_Gene_Frequency': np.mean(list(gene_frequency.values())),
        'Median_Gene_Frequency': np.median(list(gene_frequency.values())),
        'Max_Gene_Frequency': max(gene_frequency.values()),
        'Most_Frequent_Gene': max(gene_frequency.items(), key=lambda x: x[1])[0],
        'Genes_in_Single_Database': sum(1 for freq in gene_frequency.values() if freq == 1),
        'Genes_in_Multiple_Databases': sum(1 for freq in gene_frequency.values() if freq > 1),
        'High_Confidence_Genes_3plus': sum(1 for freq in gene_frequency.values() if freq >= 3)
    }
    
    stats_df = pd.DataFrame([stats_summary])
    stats_df.to_csv(report_path / f'{phenotype}_Analysis2_GeneFrequency_Statistics.csv', index=False)
    
    # Create high-confidence gene list (genes in 3+ databases)
    high_confidence_genes = gene_freq_df[gene_freq_df['Frequency'] >= 3].copy()
    if len(high_confidence_genes) > 0:
        high_confidence_genes.to_csv(report_path / f'{phenotype}_Analysis2_HighConfidenceGenes.csv', index=False)
        print(f"📈 High-confidence genes (3+ databases): {len(high_confidence_genes)}")
    
    print(f"✅ Gene frequency analysis completed")
    print(f"📊 Plot saved: GenePlots/plots/{phenotype}_Analysis2_GeneFrequency.png")
    print(f"📄 Report saved: GenePlots/reports/{phenotype}_Analysis2_GeneFrequency_Report.csv")
    print(f"🧬 Total unique genes: {len(gene_sources)}")
    print(f"🏆 Most frequent gene: {stats_summary['Most_Frequent_Gene']} (in {stats_summary['Max_Gene_Frequency']} databases)")

if __name__ == "__main__":
    print("This module should be run through main_analysis.py")
