#!/usr/bin/env python3
"""
Analysis 6: Statistical Summary
===============================
Provides comprehensive statistical summary and generates final reports
with key findings and recommendations.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from collections import Counter
import scipy.stats as stats

def run_analysis(gene_data, phenotype):
    """
    Run comprehensive statistical summary analysis
    
    Args:
        gene_data (dict): Dictionary of DataFrames with gene data
        phenotype (str): Phenotype name
    """
    print("📊 Statistical Summary Analysis")
    
    # Initialize data collection
    database_stats = {}
    all_genes = set()
    gene_sources = {}
    gene_columns = ['Gene', 'GeneSymbol', 'Gene_Symbol', 'gene_symbol', 'symbol']
    
    # Collect data from all databases
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
            genes = set(df[gene_col].dropna().unique())
            genes = {gene for gene in genes if gene != '' and pd.notna(gene)}
            
            database_stats[db_name] = {
                'total_records': len(df),
                'unique_genes': len(genes),
                'avg_records_per_gene': len(df) / len(genes) if genes else 0,
                'gene_list': genes
            }
            
            all_genes.update(genes)
            
            # Track gene sources
            for gene in genes:
                if gene not in gene_sources:
                    gene_sources[gene] = []
                gene_sources[gene].append(db_name)
    
    if not database_stats:
        print("❌ No data available for statistical analysis")
        return
    
    # Create comprehensive statistical visualization
    fig = plt.figure(figsize=(20, 16))
    fig.suptitle(f'Comprehensive Statistical Summary for {phenotype.title()}', 
                fontsize=18, fontweight='bold')
    
    # Create a complex subplot layout
    gs = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)
    
    # Plot 1: Database comparison overview
    ax1 = fig.add_subplot(gs[0, :2])
    db_names = list(database_stats.keys())
    unique_genes = [database_stats[db]['unique_genes'] for db in db_names]
    total_records = [database_stats[db]['total_records'] for db in db_names]
    
    x = np.arange(len(db_names))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, unique_genes, width, label='Unique Genes', alpha=0.8, color='skyblue')
    bars2 = ax1.bar(x + width/2, total_records, width, label='Total Records', alpha=0.8, color='lightcoral')
    
    ax1.set_title('Database Overview', fontweight='bold', fontsize=14)
    ax1.set_xlabel('Databases')
    ax1.set_ylabel('Count')
    ax1.set_xticks(x)
    ax1.set_xticklabels([name[:8] for name in db_names], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + max(max(unique_genes), max(total_records))*0.01,
                    f'{int(height)}', ha='center', va='bottom', fontsize=8)
    
    # Plot 2: Gene frequency distribution
    ax2 = fig.add_subplot(gs[0, 2:])
    gene_frequencies = [len(sources) for sources in gene_sources.values()]
    frequency_counts = Counter(gene_frequencies)
    
    freqs = list(frequency_counts.keys())
    counts = list(frequency_counts.values())
    
    bars = ax2.bar(freqs, counts, color=sns.color_palette("viridis", len(freqs)), alpha=0.8)
    ax2.set_title('Gene Frequency Distribution', fontweight='bold', fontsize=14)
    ax2.set_xlabel('Number of Databases per Gene')
    ax2.set_ylabel('Number of Genes')
    ax2.grid(True, alpha=0.3)
    
    # Add percentage labels
    total_genes = sum(counts)
    for i, bar in enumerate(bars):
        height = bar.get_height()
        percentage = (height / total_genes) * 100
        ax2.text(bar.get_x() + bar.get_width()/2., height + max(counts)*0.01,
                f'{percentage:.1f}%', ha='center', va='bottom', fontsize=8)
    
    # Plot 3: Database productivity analysis
    ax3 = fig.add_subplot(gs[1, :2])
    productivity_ratios = [database_stats[db]['avg_records_per_gene'] for db in db_names]
    
    bars = ax3.barh(range(len(db_names)), productivity_ratios, 
                   color=sns.color_palette("plasma", len(db_names)), alpha=0.8)
    ax3.set_title('Database Productivity (Records per Gene)', fontweight='bold', fontsize=14)
    ax3.set_xlabel('Average Records per Gene')
    ax3.set_yticks(range(len(db_names)))
    ax3.set_yticklabels([name[:8] for name in db_names])
    ax3.grid(True, alpha=0.3)
    
    # Add value labels
    for i, bar in enumerate(bars):
        width = bar.get_width()
        ax3.text(width + max(productivity_ratios)*0.01, bar.get_y() + bar.get_height()/2.,
                f'{width:.1f}', ha='left', va='center', fontsize=8)
    
    # Plot 4: Coverage analysis
    ax4 = fig.add_subplot(gs[1, 2:])
    
    # Calculate what percentage of total genes each database covers
    total_unique_genes = len(all_genes)
    coverage_percentages = [(database_stats[db]['unique_genes'] / total_unique_genes * 100) 
                           for db in db_names]
    
    # Create pie chart for top databases
    top_n = min(8, len(db_names))
    sorted_coverage = sorted(zip(db_names, coverage_percentages), key=lambda x: x[1], reverse=True)
    top_dbs = [db for db, _ in sorted_coverage[:top_n]]
    top_percentages = [pct for _, pct in sorted_coverage[:top_n]]
    
    if len(sorted_coverage) > top_n:
        top_dbs.append('Others')
        top_percentages.append(sum([pct for _, pct in sorted_coverage[top_n:]]))
    
    colors = sns.color_palette("Set3", len(top_dbs))
    wedges, texts, autotexts = ax4.pie(top_percentages, labels=[db[:8] for db in top_dbs], 
                                      autopct='%1.1f%%', colors=colors, startangle=90)
    ax4.set_title('Gene Coverage by Database', fontweight='bold', fontsize=14)
    
    # Plot 5: Statistical distributions
    ax5 = fig.add_subplot(gs[2, :2])
    
    # Box plot of database sizes
    db_sizes = [database_stats[db]['unique_genes'] for db in db_names]
    db_records = [database_stats[db]['total_records'] for db in db_names]
    
    box_data = [db_sizes, db_records]
    box_labels = ['Unique Genes', 'Total Records']
    
    bp = ax5.boxplot(box_data, labels=box_labels, patch_artist=True)
    colors = ['lightblue', 'lightcoral']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax5.set_title('Database Size Distributions', fontweight='bold', fontsize=14)
    ax5.set_ylabel('Count')
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Correlation analysis
    ax6 = fig.add_subplot(gs[2, 2:])
    
    # Scatter plot: unique genes vs total records
    scatter = ax6.scatter(db_sizes, db_records, s=100, alpha=0.7, 
                         c=productivity_ratios, cmap='viridis')
    
    # Add database labels
    for i, db in enumerate(db_names):
        ax6.annotate(db[:6], (db_sizes[i], db_records[i]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    # Add correlation line
    if len(db_sizes) > 1:
        slope, intercept, r_value, p_value, std_err = stats.linregress(db_sizes, db_records)
        line_x = np.array([min(db_sizes), max(db_sizes)])
        line_y = slope * line_x + intercept
        ax6.plot(line_x, line_y, 'r--', alpha=0.8, 
                label=f'R² = {r_value**2:.3f}, p = {p_value:.3f}')
        ax6.legend()
    
    ax6.set_title('Database Size Correlation', fontweight='bold', fontsize=14)
    ax6.set_xlabel('Unique Genes')
    ax6.set_ylabel('Total Records')
    ax6.grid(True, alpha=0.3)
    
    plt.colorbar(scatter, ax=ax6, label='Records/Gene Ratio')
    
    # Plot 7: Quality metrics
    ax7 = fig.add_subplot(gs[3, :2])
    
    # Calculate quality scores for each database
    quality_metrics = []
    for db in db_names:
        stats_db = database_stats[db]
        
        # Quality factors
        gene_count_score = min(stats_db['unique_genes'] / 100, 1.0)  # Max score at 100 genes
        coverage_score = (stats_db['unique_genes'] / total_unique_genes)
        productivity_score = min(stats_db['avg_records_per_gene'] / 10, 1.0)  # Max score at 10 records/gene
        
        # Combined quality score (weighted average)
        quality_score = (gene_count_score * 0.4 + coverage_score * 0.3 + productivity_score * 0.3) * 100
        quality_metrics.append(quality_score)
    
    bars = ax7.barh(range(len(db_names)), quality_metrics,
                   color=sns.color_palette("RdYlGn", len(db_names)), alpha=0.8)
    ax7.set_title('Database Quality Score', fontweight='bold', fontsize=14)
    ax7.set_xlabel('Quality Score (0-100)')
    ax7.set_yticks(range(len(db_names)))
    ax7.set_yticklabels([name[:8] for name in db_names])
    ax7.grid(True, alpha=0.3)
    
    # Add score labels
    for i, bar in enumerate(bars):
        width = bar.get_width()
        ax7.text(width + 1, bar.get_y() + bar.get_height()/2.,
                f'{width:.1f}', ha='left', va='center', fontsize=8)
    
    # Plot 8: Summary statistics table
    ax8 = fig.add_subplot(gs[3, 2:])
    ax8.axis('tight')
    ax8.axis('off')
    
    # Create summary table
    summary_data = {
        'Metric': [
            'Total Databases',
            'Total Unique Genes',
            'Total Records',
            'Avg Genes per DB',
            'Avg Records per DB',
            'Most Productive DB',
            'Largest Gene Set',
            'Gene Overlap Rate',
            'Quality Leader'
        ],
        'Value': [
            len(db_names),
            len(all_genes),
            sum(database_stats[db]['total_records'] for db in db_names),
            f"{np.mean([database_stats[db]['unique_genes'] for db in db_names]):.1f}",
            f"{np.mean([database_stats[db]['total_records'] for db in db_names]):.1f}",
            db_names[np.argmax(productivity_ratios)][:10],
            db_names[np.argmax(db_sizes)][:10],
            f"{(1 - len(all_genes) / sum(len(database_stats[db]['gene_list']) for db in db_names)) * 100:.1f}%",
            db_names[np.argmax(quality_metrics)][:10]
        ]
    }
    
    table = ax8.table(cellText=[[summary_data['Metric'][i], summary_data['Value'][i]] 
                               for i in range(len(summary_data['Metric']))],
                     colLabels=['Metric', 'Value'],
                     cellLoc='left',
                     loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    
    # Style the table
    for i in range(len(summary_data['Metric']) + 1):
        table[(i, 0)].set_facecolor('#E6E6FA')
        table[(i, 1)].set_facecolor('#F0F8FF')
    
    ax8.set_title('Summary Statistics', fontweight='bold', fontsize=14)
    
    plt.tight_layout()
    
    # Save the plot
    output_path = Path('GenePlots/plots')
    output_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path / f'{phenotype}_Analysis6_StatisticalSummary.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(output_path / f'{phenotype}_Analysis6_StatisticalSummary.pdf', 
                bbox_inches='tight')
    plt.close()
    
    # Generate comprehensive final report
    report_path = Path('GenePlots/reports')
    report_path.mkdir(parents=True, exist_ok=True)
    
    # Detailed database statistics
    detailed_stats = []
    for db in db_names:
        stats_db = database_stats[db]
        detailed_stats.append({
            'Database': db,
            'Total_Records': stats_db['total_records'],
            'Unique_Genes': stats_db['unique_genes'],
            'Records_per_Gene': stats_db['avg_records_per_gene'],
            'Coverage_Percentage': (stats_db['unique_genes'] / len(all_genes)) * 100,
            'Quality_Score': quality_metrics[db_names.index(db)]
        })
    
    detailed_df = pd.DataFrame(detailed_stats)
    detailed_df = detailed_df.sort_values('Quality_Score', ascending=False)
    detailed_df.to_csv(report_path / f'{phenotype}_Analysis6_DetailedStatistics.csv', index=False)
    
    # Gene-level analysis
    gene_analysis = []
    for gene, sources in gene_sources.items():
        gene_analysis.append({
            'Gene': gene,
            'Database_Count': len(sources),
            'Databases': ';'.join(sources),
            'Confidence_Level': 'High' if len(sources) >= 3 else 'Medium' if len(sources) == 2 else 'Low'
        })
    
    gene_df = pd.DataFrame(gene_analysis)
    gene_df = gene_df.sort_values('Database_Count', ascending=False)
    gene_df.to_csv(report_path / f'{phenotype}_Analysis6_GeneConfidenceAnalysis.csv', index=False)
    
    # Final summary report
    final_summary = {
        'Analysis_Date': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'Phenotype': phenotype,
        'Total_Databases_Analyzed': len(db_names),
        'Total_Unique_Genes_Found': len(all_genes),
        'Total_Records_Processed': sum(database_stats[db]['total_records'] for db in db_names),
        'High_Confidence_Genes': len([g for g in gene_analysis if g['Confidence_Level'] == 'High']),
        'Medium_Confidence_Genes': len([g for g in gene_analysis if g['Confidence_Level'] == 'Medium']),
        'Low_Confidence_Genes': len([g for g in gene_analysis if g['Confidence_Level'] == 'Low']),
        'Most_Productive_Database': db_names[np.argmax(productivity_ratios)],
        'Largest_Database': db_names[np.argmax(db_sizes)],
        'Highest_Quality_Database': db_names[np.argmax(quality_metrics)],
        'Average_Database_Quality': np.mean(quality_metrics),
        'Gene_Overlap_Rate': (1 - len(all_genes) / sum(len(database_stats[db]['gene_list']) for db in db_names)) * 100,
        'Recommended_Focus_Databases': ';'.join([db for db, score in zip(db_names, quality_metrics) if score > 70])
    }
    
    final_df = pd.DataFrame([final_summary])
    final_df.to_csv(report_path / f'{phenotype}_Analysis6_FinalSummary.csv', index=False)
    
    # Create executive summary text report
    executive_summary = f"""
EXECUTIVE SUMMARY - GENE ANALYSIS REPORT
========================================
Phenotype: {phenotype.title()}
Analysis Date: {final_summary['Analysis_Date']}

KEY FINDINGS:
✓ Total unique genes identified: {final_summary['Total_Unique_Genes_Found']}
✓ High-confidence genes (3+ databases): {final_summary['High_Confidence_Genes']}
✓ Databases analyzed: {final_summary['Total_Databases_Analyzed']}
✓ Total records processed: {final_summary['Total_Records_Processed']}

DATABASE PERFORMANCE:
🏆 Most productive: {final_summary['Most_Productive_Database']}
📊 Largest dataset: {final_summary['Largest_Database']}
⭐ Highest quality: {final_summary['Highest_Quality_Database']}

CONFIDENCE LEVELS:
🟢 High confidence genes: {final_summary['High_Confidence_Genes']} ({final_summary['High_Confidence_Genes']/len(all_genes)*100:.1f}%)
🟡 Medium confidence genes: {final_summary['Medium_Confidence_Genes']} ({final_summary['Medium_Confidence_Genes']/len(all_genes)*100:.1f}%)
🔴 Low confidence genes: {final_summary['Low_Confidence_Genes']} ({final_summary['Low_Confidence_Genes']/len(all_genes)*100:.1f}%)

RECOMMENDATIONS:
1. Focus on high-confidence genes for further validation
2. Prioritize databases: {final_summary['Recommended_Focus_Databases'].replace(';', ', ')}
3. Consider additional validation for medium-confidence genes
4. Gene overlap rate of {final_summary['Gene_Overlap_Rate']:.1f}% indicates good database complementarity

For detailed analysis, see individual report files and visualizations.
"""
    
    with open(report_path / f'{phenotype}_ExecutiveSummary.txt', 'w') as f:
        f.write(executive_summary)
    
    print(f"✅ Statistical summary analysis completed")
    print(f"📊 Plot saved: GenePlots/plots/{phenotype}_Analysis6_StatisticalSummary.png")
    print(f"📄 Reports saved: GenePlots/reports/{phenotype}_Analysis6_*.csv")
    print(f"📋 Executive summary: GenePlots/reports/{phenotype}_ExecutiveSummary.txt")
    print(f"🎯 Key findings:")
    print(f"   • {len(all_genes)} unique genes found")
    print(f"   • {final_summary['High_Confidence_Genes']} high-confidence genes")
    print(f"   • Best database: {final_summary['Highest_Quality_Database']}")

if __name__ == "__main__":
    print("This module should be run through main_analysis.py")
