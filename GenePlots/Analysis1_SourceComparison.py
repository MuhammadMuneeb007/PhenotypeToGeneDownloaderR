#!/usr/bin/env python3
"""
Analysis 1: Source Comparison Analysis
=====================================
Compares gene counts and associations across different databases.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

def run_analysis(gene_data, phenotype):
    """
    Run source comparison analysis
    
    Args:
        gene_data (dict): Dictionary of DataFrames with gene data
        phenotype (str): Phenotype name
    """
    print("📊 Source Comparison Analysis")
    
    # Extract summary data
    if 'SUMMARY' not in gene_data:
        print("❌ Summary data not found")
        return
    
    summary_df = gene_data['SUMMARY'].copy()
    
    # Create comprehensive comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'Database Source Comparison for {phenotype.title()}', fontsize=16, fontweight='bold')
    
    # Plot 1: Gene counts by source (bar plot)
    ax1 = axes[0, 0]
    bars = ax1.bar(range(len(summary_df)), summary_df['Unique_Genes'], 
                   color=sns.color_palette("viridis", len(summary_df)))
    ax1.set_title('Unique Genes per Database', fontweight='bold')
    ax1.set_xlabel('Database')
    ax1.set_ylabel('Number of Unique Genes')
    ax1.set_xticks(range(len(summary_df)))
    ax1.set_xticklabels(summary_df['Source'], rotation=45, ha='right')
    
    # Add value labels on bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 2: Total associations by source
    ax2 = axes[0, 1]
    bars2 = ax2.bar(range(len(summary_df)), summary_df['Total_Associations'],
                    color=sns.color_palette("plasma", len(summary_df)))
    ax2.set_title('Total Associations per Database', fontweight='bold')
    ax2.set_xlabel('Database')
    ax2.set_ylabel('Number of Associations')
    ax2.set_xticks(range(len(summary_df)))
    ax2.set_xticklabels(summary_df['Source'], rotation=45, ha='right')
    
    # Add value labels
    for i, bar in enumerate(bars2):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + max(summary_df['Total_Associations'])*0.01,
                f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 3: Ratio of associations to unique genes
    ax3 = axes[1, 0]
    ratio = summary_df['Total_Associations'] / summary_df['Unique_Genes']
    bars3 = ax3.bar(range(len(summary_df)), ratio,
                    color=sns.color_palette("coolwarm", len(summary_df)))
    ax3.set_title('Association-to-Gene Ratio', fontweight='bold')
    ax3.set_xlabel('Database')
    ax3.set_ylabel('Associations per Gene')
    ax3.set_xticks(range(len(summary_df)))
    ax3.set_xticklabels(summary_df['Source'], rotation=45, ha='right')
    
    # Add value labels
    for i, bar in enumerate(bars3):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{height:.1f}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 4: Data completeness heatmap
    ax4 = axes[1, 1]
    
    # Create completeness matrix
    databases = summary_df['Source'].tolist()
    completeness_data = []
    
    for db in databases:
        has_data = []
        for db_check in databases:
            if db == db_check:
                has_data.append(1.0)  # Self comparison
            else:
                # Check if both databases have data
                db_genes = summary_df[summary_df['Source'] == db]['Unique_Genes'].iloc[0]
                check_genes = summary_df[summary_df['Source'] == db_check]['Unique_Genes'].iloc[0]
                has_data.append(min(db_genes, check_genes) / max(db_genes, check_genes) if max(db_genes, check_genes) > 0 else 0)
        completeness_data.append(has_data)
    
    sns.heatmap(completeness_data, 
                xticklabels=[s[:8] for s in databases],
                yticklabels=[s[:8] for s in databases],
                annot=True, fmt='.2f', cmap='RdYlBu_r', ax=ax4)
    ax4.set_title('Database Size Similarity Matrix', fontweight='bold')
    
    plt.tight_layout()
    
    # Save the plot
    output_path = Path('GenePlots/plots')
    output_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path / f'{phenotype}_Analysis1_SourceComparison.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(output_path / f'{phenotype}_Analysis1_SourceComparison.pdf', 
                bbox_inches='tight')
    plt.close()
    
    # Create detailed statistics report
    report_data = {
        'Database': summary_df['Source'],
        'Unique_Genes': summary_df['Unique_Genes'],
        'Total_Associations': summary_df['Total_Associations'],
        'Avg_Associations_per_Gene': summary_df['Total_Associations'] / summary_df['Unique_Genes'],
        'Percentage_of_Total_Genes': (summary_df['Unique_Genes'] / summary_df['Unique_Genes'].sum() * 100).round(2)
    }
    
    detailed_report = pd.DataFrame(report_data)
    detailed_report = detailed_report.sort_values('Unique_Genes', ascending=False)
    
    # Save detailed report
    report_path = Path('GenePlots/reports')
    report_path.mkdir(parents=True, exist_ok=True)
    detailed_report.to_csv(report_path / f'{phenotype}_Analysis1_SourceComparison_Report.csv', index=False)
    
    # Save summary statistics
    stats_summary = {
        'Total_Databases': len(summary_df),
        'Total_Unique_Genes_Across_All': summary_df['Unique_Genes'].sum(),
        'Average_Genes_per_Database': summary_df['Unique_Genes'].mean(),
        'Median_Genes_per_Database': summary_df['Unique_Genes'].median(),
        'Most_Productive_Database': summary_df.loc[summary_df['Unique_Genes'].idxmax(), 'Source'],
        'Least_Productive_Database': summary_df.loc[summary_df['Unique_Genes'].idxmin(), 'Source'],
        'Total_Associations_All_Databases': summary_df['Total_Associations'].sum()
    }
    
    stats_df = pd.DataFrame([stats_summary])
    stats_df.to_csv(report_path / f'{phenotype}_Analysis1_Summary_Statistics.csv', index=False)
    
    print(f"✅ Source comparison analysis completed")
    print(f"📊 Plot saved: GenePlots/plots/{phenotype}_Analysis1_SourceComparison.png")
    print(f"📄 Report saved: GenePlots/reports/{phenotype}_Analysis1_SourceComparison_Report.csv")
    print(f"📈 Top database: {stats_summary['Most_Productive_Database']} ({summary_df['Unique_Genes'].max()} genes)")

if __name__ == "__main__":
    # Test the module independently
    print("This module should be run through main_analysis.py")
