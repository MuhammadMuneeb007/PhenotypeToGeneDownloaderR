#!/usr/bin/env python3
"""
Analysis 3: Database Overlap Analysis
=====================================
Analyzes the overlap between different databases using Venn diagrams,
UpSet plots, and Jaccard similarity measures.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from itertools import combinations
import matplotlib.patches as patches
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles

def run_analysis(gene_data, phenotype):
    """
    Run database overlap analysis
    
    Args:
        gene_data (dict): Dictionary of DataFrames with gene data
        phenotype (str): Phenotype name
    """
    print("🔗 Database Overlap Analysis")
    
    # Extract gene sets from each database
    gene_sets = {}
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
            genes = set(df[gene_col].dropna().unique())
            genes = {gene for gene in genes if gene != '' and pd.notna(gene)}
            if genes:
                gene_sets[db_name] = genes
    
    if len(gene_sets) < 2:
        print("❌ Need at least 2 databases for overlap analysis")
        return
    
    # Calculate Jaccard similarities first
    db_names = list(gene_sets.keys())
    n_dbs = len(db_names)
    jaccard_matrix = np.zeros((n_dbs, n_dbs))
    
    for i, db1 in enumerate(db_names):
        for j, db2 in enumerate(db_names):
            if i == j:
                jaccard_matrix[i, j] = 1.0
            else:
                intersection = len(gene_sets[db1] & gene_sets[db2])
                union = len(gene_sets[db1] | gene_sets[db2])
                jaccard_matrix[i, j] = intersection / union if union > 0 else 0
    
    # Create publication-quality heatmap visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'Database Overlap Analysis for {phenotype.title()}', 
                fontsize=16, fontweight='bold')
    
    # Plot 1: Jaccard Similarity Heatmap
    mask = np.triu(np.ones_like(jaccard_matrix), k=1)
    
    # Create custom colormap for better publication quality
    cmap = sns.diverging_palette(250, 10, as_cmap=True)
    
    heatmap = sns.heatmap(jaccard_matrix, 
                         annot=True, 
                         fmt='.3f',
                         cmap=cmap,
                         square=True,
                         linewidths=0.5,
                         cbar_kws={"shrink": .8, "label": "Jaccard Similarity"},
                         xticklabels=[name[:8] for name in db_names],
                         yticklabels=[name[:8] for name in db_names],
                         ax=ax1,
                         mask=mask)
    
    ax1.set_title('Database Similarity (Jaccard Index)', fontweight='bold', fontsize=14)
    ax1.set_xlabel('Databases', fontweight='bold')
    ax1.set_ylabel('Databases', fontweight='bold')
    
    # Rotate labels for better readability
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
    ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0)
    
    # Plot 2: Database sizes
    db_sizes = [len(gene_sets[db]) for db in db_names]
    bars = ax2.bar(range(len(db_names)), db_sizes, 
                   color=sns.color_palette("viridis", len(db_names)), alpha=0.8)
    ax2.set_title('Database Sizes', fontweight='bold', fontsize=14)
    ax2.set_xlabel('Databases', fontweight='bold')
    ax2.set_ylabel('Number of Genes', fontweight='bold')
    ax2.set_xticks(range(len(db_names)))
    ax2.set_xticklabels([name[:8] for name in db_names], rotation=45, ha='right')
    ax2.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + max(db_sizes)*0.01,
                f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    # Plot 3: Top pairwise overlaps
    overlap_data = []
    for i, db1 in enumerate(db_names):
        for j, db2 in enumerate(db_names):
            if i < j:  # Only upper triangle
                intersection = len(gene_sets[db1] & gene_sets[db2])
                overlap_data.append({
                    'Pair': f'{db1[:6]}−{db2[:6]}',
                    'Overlap': intersection,
                    'Jaccard': jaccard_matrix[i, j]
                })
    
    overlap_df = pd.DataFrame(overlap_data)
    overlap_df = overlap_df.sort_values('Overlap', ascending=True).tail(10)  # Top 10
    
    if len(overlap_df) > 0:
        bars3 = ax3.barh(range(len(overlap_df)), overlap_df['Overlap'],
                         color=sns.color_palette("plasma", len(overlap_df)))
        ax3.set_title('Top 10 Pairwise Overlaps', fontweight='bold', fontsize=14)
        ax3.set_yticks(range(len(overlap_df)))
        ax3.set_yticklabels(overlap_df['Pair'])
        ax3.set_xlabel('Number of Shared Genes', fontweight='bold')
        ax3.grid(True, alpha=0.3)
        
        # Add value labels
        for i, bar in enumerate(bars3):
            width = bar.get_width()
            ax3.text(width + max(overlap_df['Overlap'])*0.01, bar.get_y() + bar.get_height()/2.,
                    f'{int(width)}', ha='left', va='center', fontweight='bold')
    
    # Plot 4: Clustering dendrogram
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform
    
    # Convert similarity to distance for clustering
    distance_matrix = 1 - jaccard_matrix
    np.fill_diagonal(distance_matrix, 0)
    
    # Convert to condensed form for linkage
    condensed_distances = squareform(distance_matrix)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_distances, method='average')
    
    # Create dendrogram
    dendro = dendrogram(linkage_matrix, 
                       labels=[name[:8] for name in db_names],
                       ax=ax4,
                       orientation='top',
                       leaf_rotation=45)
    
    ax4.set_title('Database Clustering (Based on Gene Overlap)', fontweight='bold', fontsize=14)
    ax4.set_ylabel('Distance (1 - Jaccard)', fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    # Find genes present in multiple databases
    all_genes = set()
    for gene_set in gene_sets.values():
        all_genes.update(gene_set)
    
    gene_presence = {}
    for gene in all_genes:
        count = sum(1 for gene_set in gene_sets.values() if gene in gene_set)
        gene_presence[gene] = count
    
    plt.tight_layout()
    
    # Save the main plot
    output_path = Path('GenePlots/plots')
    output_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path / f'{phenotype}_Analysis3_DatabaseOverlap.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(output_path / f'{phenotype}_Analysis3_DatabaseOverlap.pdf', 
                bbox_inches='tight')
    plt.close()
    
    # Create a separate publication-quality heatmap
    create_similarity_heatmap(jaccard_matrix, db_names, phenotype)
    
    # Generate additional Venn diagrams for top pairs
    create_venn_diagrams(gene_sets, overlap_df, phenotype)
    
    # Create comprehensive Venn diagrams for all databases
    create_comprehensive_venn_analysis(gene_sets, phenotype)
    
    # Save Jaccard matrix to CSV
    report_path = Path('GenePlots/reports')
    report_path.mkdir(parents=True, exist_ok=True)
    
    # Create a DataFrame for the Jaccard matrix
    jaccard_df = pd.DataFrame(jaccard_matrix, 
                             index=db_names, 
                             columns=db_names)
    jaccard_df.to_csv(report_path / f'{phenotype}_Analysis3_JaccardMatrix.csv')
    
    # Save overlap analysis results
    overlap_df.to_csv(report_path / f'{phenotype}_Analysis3_PairwiseOverlaps.csv', index=False)
    
    # Create summary statistics
    summary_stats = {
        'Total_Databases': len(gene_sets),
        'Total_Unique_Genes': len(all_genes),
        'Average_Database_Size': np.mean([len(s) for s in gene_sets.values()]),
        'Average_Jaccard_Similarity': np.mean(jaccard_matrix[np.triu_indices_from(jaccard_matrix, k=1)]),
        'Max_Jaccard_Similarity': np.max(jaccard_matrix[np.triu_indices_from(jaccard_matrix, k=1)]),
        'Min_Jaccard_Similarity': np.min(jaccard_matrix[np.triu_indices_from(jaccard_matrix, k=1)]),
        'Most_Similar_Pair': f"{db_names[np.unravel_index(np.argmax(jaccard_matrix + np.eye(len(db_names))*-1), jaccard_matrix.shape)[0]]}-{db_names[np.unravel_index(np.argmax(jaccard_matrix + np.eye(len(db_names))*-1), jaccard_matrix.shape)[1]]}",
        'Core_Genes_Count': sum(1 for count in gene_presence.values() if count >= len(gene_sets) * 0.5)
    }
    
    summary_df = pd.DataFrame([summary_stats])
    summary_df.to_csv(report_path / f'{phenotype}_Analysis3_OverlapSummary.csv', index=False)
    
    print(f"✅ Database overlap analysis completed")
    print(f"📊 Plot saved: GenePlots/plots/{phenotype}_Analysis3_DatabaseOverlap.png")
    print(f"🔥 Heatmap saved: GenePlots/plots/{phenotype}_Analysis3_SimilarityHeatmap.png")
    print(f"📄 Reports saved: GenePlots/reports/{phenotype}_Analysis3_*.csv")
    print(f"🎯 Key findings:")
    print(f"   • Average Jaccard similarity: {summary_stats['Average_Jaccard_Similarity']:.3f}")
    print(f"   • Most similar pair: {summary_stats['Most_Similar_Pair']}")
    print(f"   • Core genes (in ≥50% databases): {summary_stats['Core_Genes_Count']}")

def create_similarity_heatmap(jaccard_matrix, db_names, phenotype):
    """Create a standalone publication-quality similarity heatmap"""
    # Create a large, detailed heatmap
    plt.figure(figsize=(12, 10))
    
    # Create a custom colormap for better visualization
    cmap = sns.diverging_palette(250, 10, as_cmap=True)
    
    # Create the heatmap with enhanced styling
    mask = np.triu(np.ones_like(jaccard_matrix), k=1)
    
    ax = sns.heatmap(jaccard_matrix, 
                     annot=True, 
                     fmt='.3f',
                     cmap=cmap,
                     square=True,
                     linewidths=0.8,
                     linecolor='white',
                     cbar_kws={
                         "shrink": 0.8, 
                         "label": "Jaccard Similarity Coefficient",
                         "orientation": "vertical"
                     },
                     xticklabels=[name[:10] for name in db_names],
                     yticklabels=[name[:10] for name in db_names],
                     mask=mask,
                     annot_kws={'size': 10, 'weight': 'bold'})
    
    # Enhance the plot styling for publication quality
    plt.title(f'Database Similarity Matrix\n{phenotype.title()} Gene Overlap (Jaccard Index)', 
              fontsize=16, fontweight='bold', pad=20)
    
    plt.xlabel('Databases', fontsize=14, fontweight='bold')
    plt.ylabel('Databases', fontsize=14, fontweight='bold')
    
    # Rotate labels for better readability
    plt.xticks(rotation=45, ha='right', fontsize=11)
    plt.yticks(rotation=0, fontsize=11)
    
    # Add grid for better readability
    ax.grid(False)
    
    # Adjust layout for publication quality
    plt.tight_layout()
    
    # Save with high quality
    output_path = Path('GenePlots/plots')
    output_path.mkdir(parents=True, exist_ok=True)
    
    plt.savefig(output_path / f'{phenotype}_Analysis3_SimilarityHeatmap.png', 
                dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(output_path / f'{phenotype}_Analysis3_SimilarityHeatmap.pdf', 
                bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def create_venn_diagrams(gene_sets, overlap_df, phenotype):
    """Create publication-quality Venn diagrams for top database pairs"""
    db_names = list(gene_sets.keys())
    
    # Get top 6 pairs for 2x3 subplot layout
    top_pairs = overlap_df.nlargest(6, 'Overlap') if len(overlap_df) >= 6 else overlap_df
    
    if len(top_pairs) == 0:
        return
    
    # Create Venn diagrams
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'Top Database Pair Overlaps - {phenotype.title()}', 
                fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    for idx, (_, row) in enumerate(top_pairs.iterrows()):
        if idx >= 6:
            break
            
        ax = axes[idx]
        pair_name = row['Pair']
        db1_name, db2_name = pair_name.split('−')
        
        # Find full database names
        db1_full = next((name for name in db_names if name.startswith(db1_name)), db1_name)
        db2_full = next((name for name in db_names if name.startswith(db2_name)), db2_name)
        
        if db1_full in gene_sets and db2_full in gene_sets:
            set1 = gene_sets[db1_full]
            set2 = gene_sets[db2_full]
            
            # Create Venn diagram
            venn = venn2([set1, set2], 
                        (db1_full[:10], db2_full[:10]), 
                        ax=ax)
            
            if venn:
                # Style the Venn diagram
                venn2_circles([set1, set2], alpha=0.6, ax=ax)
                
                # Add styling
                if venn.get_patch_by_id('10'):
                    venn.get_patch_by_id('10').set_color('lightblue')
                if venn.get_patch_by_id('01'):
                    venn.get_patch_by_id('01').set_color('lightcoral')
                if venn.get_patch_by_id('11'):
                    venn.get_patch_by_id('11').set_color('lightgreen')
            
            # Add title with statistics
            overlap_count = len(set1 & set2)
            jaccard_sim = row['Jaccard']
            ax.set_title(f'{db1_full[:8]} ∩ {db2_full[:8]}\n'
                        f'Overlap: {overlap_count} genes\n'
                        f'Jaccard: {jaccard_sim:.3f}', 
                        fontweight='bold', fontsize=12)
    
    # Hide empty subplots
    for idx in range(len(top_pairs), 6):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    
    # Save the Venn diagrams
    output_path = Path('GenePlots/plots')
    plt.savefig(output_path / f'{phenotype}_Analysis3_VennDiagrams.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(output_path / f'{phenotype}_Analysis3_VennDiagrams.pdf', 
                bbox_inches='tight')
    plt.close()

def create_comprehensive_venn_analysis(gene_sets, phenotype):
    """Create comprehensive Venn diagram analysis for all databases"""
    db_names = list(gene_sets.keys())
    
    # Create multiple Venn diagram layouts based on number of databases
    if len(db_names) < 2:
        print("Need at least 2 databases for Venn diagrams")
        return
    
    # Strategy: Create multiple plots showing different combinations
    # 1. Top databases (largest gene sets)
    # 2. All pairwise combinations
    # 3. Three-way combinations for top databases
    
    # Sort databases by size for better visualization
    db_sizes = [(name, len(gene_sets[name])) for name in db_names]
    db_sizes.sort(key=lambda x: x[1], reverse=True)
    
    # Create comprehensive Venn analysis
    create_top_database_venns(gene_sets, db_sizes, phenotype)
    create_all_pairwise_venns(gene_sets, db_names, phenotype)
    if len(db_names) >= 3:
        create_threeway_venns(gene_sets, db_sizes, phenotype)

def create_top_database_venns(gene_sets, db_sizes, phenotype):
    """Create Venn diagrams for top databases"""
    # Take top 6 databases by size
    top_dbs = [name for name, _ in db_sizes[:6]]
    
    if len(top_dbs) < 2:
        return
    
    # Create 2-way Venn diagrams for top databases
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    fig.suptitle(f'Top Databases Venn Analysis - {phenotype.title()}', 
                fontsize=18, fontweight='bold', y=0.98)
    
    axes = axes.flatten()
    
    # Show top database vs each other
    main_db = top_dbs[0]
    comparisons = 0
    
    for i, other_db in enumerate(top_dbs[1:6]):
        if comparisons >= 6:
            break
            
        ax = axes[comparisons]
        
        set1 = gene_sets[main_db]
        set2 = gene_sets[other_db]
        
        # Create Venn diagram
        venn = venn2([set1, set2], 
                    (main_db[:12], other_db[:12]), 
                    ax=ax)
        
        if venn:
            # Style the Venn diagram
            venn2_circles([set1, set2], alpha=0.7, ax=ax)
            
            # Custom colors
            colors = ['lightblue', 'lightcoral', 'lightgreen']
            if venn.get_patch_by_id('10'):
                venn.get_patch_by_id('10').set_color(colors[0])
                venn.get_patch_by_id('10').set_alpha(0.7)
            if venn.get_patch_by_id('01'):
                venn.get_patch_by_id('01').set_color(colors[1])
                venn.get_patch_by_id('01').set_alpha(0.7)
            if venn.get_patch_by_id('11'):
                venn.get_patch_by_id('11').set_color(colors[2])
                venn.get_patch_by_id('11').set_alpha(0.8)
        
        # Calculate statistics
        overlap_count = len(set1 & set2)
        jaccard_sim = overlap_count / len(set1 | set2) if len(set1 | set2) > 0 else 0
        
        # Add detailed title with statistics
        ax.set_title(f'{main_db[:10]} ∩ {other_db[:10]}\n'
                    f'Overlap: {overlap_count} genes\n'
                    f'Jaccard: {jaccard_sim:.3f}\n'
                    f'Coverage: {overlap_count/len(set1)*100:.1f}% | {overlap_count/len(set2)*100:.1f}%', 
                    fontweight='bold', fontsize=11, pad=10)
        
        comparisons += 1
    
    # Hide unused subplots
    for i in range(comparisons, 6):
        axes[i].set_visible(False)
    
    # Add legend
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor='lightblue', alpha=0.7, label=f'{main_db} only'),
        plt.Rectangle((0,0),1,1, facecolor='lightcoral', alpha=0.7, label='Other DB only'),
        plt.Rectangle((0,0),1,1, facecolor='lightgreen', alpha=0.8, label='Shared genes')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=12, 
              bbox_to_anchor=(0.5, 0.02))
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1)
    
    # Save the plot
    output_path = Path('GenePlots/plots')
    plt.savefig(output_path / f'{phenotype}_Analysis3_TopDatabaseVenns.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(output_path / f'{phenotype}_Analysis3_TopDatabaseVenns.pdf', 
                bbox_inches='tight')
    plt.close()

def create_all_pairwise_venns(gene_sets, db_names, phenotype):
    """Create a comprehensive view of all pairwise overlaps"""
    from itertools import combinations
    
    # Get all pairwise combinations
    pairs = list(combinations(db_names, 2))
    
    if len(pairs) == 0:
        return
    
    # Determine grid size
    n_pairs = len(pairs)
    if n_pairs <= 6:
        rows, cols = 2, 3
    elif n_pairs <= 12:
        rows, cols = 3, 4
    elif n_pairs <= 20:
        rows, cols = 4, 5
    else:
        rows, cols = 5, 6
        pairs = pairs[:30]  # Limit to 30 for readability
    
    fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 3*rows))
    fig.suptitle(f'All Database Pairwise Overlaps - {phenotype.title()}', 
                fontsize=16, fontweight='bold')
    
    if rows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()
    
    # Color palette for databases
    colors = plt.cm.Set3(np.linspace(0, 1, len(db_names)))
    db_colors = dict(zip(db_names, colors))
    
    for i, (db1, db2) in enumerate(pairs):
        if i >= len(axes):
            break
            
        ax = axes[i]
        
        set1 = gene_sets[db1]
        set2 = gene_sets[db2]
        
        # Create Venn diagram
        venn = venn2([set1, set2], 
                    (db1[:8], db2[:8]), 
                    ax=ax)
        
        if venn:
            # Use database-specific colors
            venn2_circles([set1, set2], alpha=0.6, ax=ax)
            
            if venn.get_patch_by_id('10'):
                venn.get_patch_by_id('10').set_color(db_colors[db1])
                venn.get_patch_by_id('10').set_alpha(0.6)
            if venn.get_patch_by_id('01'):
                venn.get_patch_by_id('01').set_color(db_colors[db2])
                venn.get_patch_by_id('01').set_alpha(0.6)
            if venn.get_patch_by_id('11'):
                venn.get_patch_by_id('11').set_color('gold')
                venn.get_patch_by_id('11').set_alpha(0.8)
        
        # Calculate overlap
        overlap_count = len(set1 & set2)
        jaccard_sim = overlap_count / len(set1 | set2) if len(set1 | set2) > 0 else 0
        
        ax.set_title(f'{db1[:6]} ∩ {db2[:6]}\n{overlap_count} genes\nJ={jaccard_sim:.2f}', 
                    fontweight='bold', fontsize=9)
    
    # Hide unused subplots
    for i in range(len(pairs), len(axes)):
        axes[i].set_visible(False)
    
    # Create comprehensive legend
    legend_elements = []
    for db_name in db_names[:10]:  # Show up to 10 databases in legend
        legend_elements.append(
            plt.Rectangle((0,0),1,1, facecolor=db_colors[db_name], alpha=0.6, 
                         label=db_name[:10])
        )
    legend_elements.append(
        plt.Rectangle((0,0),1,1, facecolor='gold', alpha=0.8, label='Overlap')
    )
    
    # Place legend outside the plot
    fig.legend(handles=legend_elements, loc='center right', 
              bbox_to_anchor=(1.02, 0.5), fontsize=10)
    
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
    
    # Save the plot
    output_path = Path('GenePlots/plots')
    plt.savefig(output_path / f'{phenotype}_Analysis3_AllPairwiseVenns.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(output_path / f'{phenotype}_Analysis3_AllPairwiseVenns.pdf', 
                bbox_inches='tight')
    plt.close()

def create_threeway_venns(gene_sets, db_sizes, phenotype):
    """Create three-way Venn diagrams for top database combinations"""
    # Take top databases for 3-way analysis
    top_dbs = [name for name, _ in db_sizes[:min(6, len(db_sizes))]]
    
    if len(top_dbs) < 3:
        return
    
    from itertools import combinations
    
    # Get all 3-way combinations from top databases
    triplets = list(combinations(top_dbs, 3))
    
    # Limit to 6 triplets for readability
    triplets = triplets[:6]
    
    if len(triplets) == 0:
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'Three-Way Database Overlaps - {phenotype.title()}', 
                fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    # Color palette
    colors = ['lightblue', 'lightcoral', 'lightgreen', 'yellow', 'orange', 'purple', 'pink']
    
    for i, (db1, db2, db3) in enumerate(triplets):
        if i >= 6:
            break
            
        ax = axes[i]
        
        set1 = gene_sets[db1]
        set2 = gene_sets[db2]
        set3 = gene_sets[db3]
        
        # Create 3-way Venn diagram
        try:
            venn = venn3([set1, set2, set3], 
                        (db1[:8], db2[:8], db3[:8]), 
                        ax=ax)
            
            if venn:
                venn3_circles([set1, set2, set3], alpha=0.6, ax=ax)
                
                # Apply colors to different regions
                region_colors = {
                    '100': colors[0], '010': colors[1], '001': colors[2],
                    '110': colors[3], '101': colors[4], '011': colors[5],
                    '111': colors[6]
                }
                
                for region_id, color in region_colors.items():
                    if venn.get_patch_by_id(region_id):
                        venn.get_patch_by_id(region_id).set_color(color)
                        venn.get_patch_by_id(region_id).set_alpha(0.7)
            
        except Exception as e:
            # Fallback if 3-way Venn fails
            ax.text(0.5, 0.5, f'3-way Venn\n{db1[:6]} ∩ {db2[:6]} ∩ {db3[:6]}\n'
                             f'Common: {len(set1 & set2 & set3)} genes', 
                   ha='center', va='center', transform=ax.transAxes,
                   fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
        
        # Calculate three-way overlap
        triple_overlap = len(set1 & set2 & set3)
        
        ax.set_title(f'{db1[:6]} ∩ {db2[:6]} ∩ {db3[:6]}\n'
                    f'Triple overlap: {triple_overlap} genes', 
                    fontweight='bold', fontsize=10)
    
    # Hide unused subplots
    for i in range(len(triplets), 6):
        axes[i].set_visible(False)
    
    # Add legend for 3-way Venn
    legend_labels = [
        'DB1 only', 'DB2 only', 'DB3 only',
        'DB1∩DB2', 'DB1∩DB3', 'DB2∩DB3', 
        'All three'
    ]
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor=colors[i], alpha=0.7, label=legend_labels[i])
        for i in range(7)
    ]
    
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=10,
              bbox_to_anchor=(0.5, 0.02))
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)
    
    # Save the plot
    output_path = Path('GenePlots/plots')
    plt.savefig(output_path / f'{phenotype}_Analysis3_ThreeWayVenns.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(output_path / f'{phenotype}_Analysis3_ThreeWayVenns.pdf', 
                bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    print("This module should be run through main_analysis.py")
