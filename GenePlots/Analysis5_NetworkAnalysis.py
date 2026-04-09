#!/usr/bin/env python3
"""
Analysis 5: Network Analysis
============================
Creates protein-protein interaction networks and analyzes network properties
of the genes associated with the phenotype.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter
import networkx as nx

def run_analysis(gene_data, phenotype):
    """
    Run network analysis
    
    Args:
        gene_data (dict): Dictionary of DataFrames with gene data
        phenotype (str): Phenotype name
    """
    print("🕸️ Network Analysis")
    
    # Get all unique genes
    all_genes = set()
    gene_databases = defaultdict(list)  # Track which databases each gene comes from
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
            all_genes.update(genes)
            
            # Track database sources for each gene
            for gene in genes:
                gene_databases[gene].append(db_name)
    
    if len(all_genes) < 2:
        print("❌ Need at least 2 genes for network analysis")
        return
    
    gene_list = list(all_genes)
    print(f"🧬 Analyzing network for {len(gene_list)} genes")
    
    # Create mock protein interaction network (in real scenario, use STRING or other PPI databases)
    G = create_mock_network(gene_list, gene_databases)
    
    # Create comprehensive network visualization
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'Network Analysis for {phenotype.title()}', fontsize=16, fontweight='bold')
    
    # Plot 1: Main network visualization
    ax1 = axes[0, 0]
    pos = nx.spring_layout(G, k=1, iterations=50, seed=42)
    
    # Color nodes by database frequency
    node_colors = []
    for node in G.nodes():
        db_count = len(gene_databases.get(node, []))
        node_colors.append(db_count)
    
    nodes = nx.draw_networkx_nodes(G, pos, node_color=node_colors, 
                                  cmap='viridis', node_size=50, ax=ax1)
    nx.draw_networkx_edges(G, pos, alpha=0.3, width=0.5, ax=ax1)
    
    # Add labels for high-degree nodes
    high_degree_nodes = [node for node, degree in G.degree() if degree >= 5]
    labels = {node: node for node in high_degree_nodes}
    nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax1)
    
    ax1.set_title('Protein Interaction Network', fontweight='bold')
    ax1.axis('off')
    
    # Add colorbar
    plt.colorbar(nodes, ax=ax1, label='Database Count')
    
    # Plot 2: Degree distribution
    ax2 = axes[0, 1]
    degrees = [G.degree(node) for node in G.nodes()]
    
    ax2.hist(degrees, bins=min(20, max(degrees)), alpha=0.7, 
             color='skyblue', edgecolor='black')
    ax2.set_title('Degree Distribution', fontweight='bold')
    ax2.set_xlabel('Node Degree')
    ax2.set_ylabel('Frequency')
    ax2.axvline(x=np.mean(degrees), color='red', linestyle='--', 
               label=f'Mean: {np.mean(degrees):.1f}')
    ax2.legend()
    
    # Plot 3: Top hub genes
    ax3 = axes[0, 2]
    degree_dict = dict(G.degree())
    top_hubs = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)[:15]
    
    if top_hubs:
        hub_genes = [gene for gene, degree in top_hubs]
        hub_degrees = [degree for gene, degree in top_hubs]
        
        bars = ax3.barh(range(len(hub_genes)), hub_degrees,
                       color=sns.color_palette("plasma", len(hub_genes)))
        ax3.set_title('Top Hub Genes', fontweight='bold')
        ax3.set_xlabel('Number of Connections')
        ax3.set_ylabel('Genes')
        ax3.set_yticks(range(len(hub_genes)))
        ax3.set_yticklabels(hub_genes)
        ax3.invert_yaxis()
    
    # Plot 4: Network communities
    ax4 = axes[1, 0]
    
    try:
        # Find communities using Louvain algorithm
        communities = nx.algorithms.community.louvain_communities(G, seed=42)
        
        # Color nodes by community
        community_colors = plt.cm.Set3(np.linspace(0, 1, len(communities)))
        node_colors_comm = {}
        
        for i, community in enumerate(communities):
            for node in community:
                node_colors_comm[node] = community_colors[i]
        
        node_color_list = [node_colors_comm.get(node, 'gray') for node in G.nodes()]
        
        nx.draw(G, pos, node_color=node_color_list, node_size=30, 
               with_labels=False, edge_color='gray', alpha=0.7, ax=ax4)
        ax4.set_title(f'Network Communities ({len(communities)} found)', fontweight='bold')
        
        # Create community size chart
        community_sizes = [len(comm) for comm in communities]
        ax4.text(0.02, 0.98, f'Sizes: {community_sizes}', 
                transform=ax4.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
    except:
        ax4.text(0.5, 0.5, 'Community detection\nnot available', 
                ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title('Communities (NetworkX Required)', fontweight='bold')
    
    # Plot 5: Network metrics over database representation
    ax5 = axes[1, 1]
    
    # Calculate metrics for genes by database count
    db_counts = [len(gene_databases.get(node, [])) for node in G.nodes()]
    degrees = [G.degree(node) for node in G.nodes()]
    
    # Create scatter plot
    scatter = ax5.scatter(db_counts, degrees, alpha=0.6, s=30)
    ax5.set_title('Database Count vs Network Degree', fontweight='bold')
    ax5.set_xlabel('Number of Databases')
    ax5.set_ylabel('Network Degree')
    
    # Add correlation coefficient
    if len(db_counts) > 1 and len(set(db_counts)) > 1:
        correlation = np.corrcoef(db_counts, degrees)[0, 1]
        ax5.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                transform=ax5.transAxes, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Plot 6: Centrality measures
    ax6 = axes[1, 2]
    
    # Calculate different centrality measures for top nodes
    try:
        betweenness = nx.betweenness_centrality(G)
        closeness = nx.closeness_centrality(G)
        eigenvector = nx.eigenvector_centrality(G, max_iter=1000)
        
        # Get top 10 nodes by betweenness centrality
        top_betweenness = sorted(betweenness.items(), key=lambda x: x[1], reverse=True)[:10]
        
        if top_betweenness:
            genes = [gene for gene, _ in top_betweenness]
            betw_values = [betweenness[gene] for gene in genes]
            close_values = [closeness[gene] for gene in genes]
            eigen_values = [eigenvector[gene] for gene in genes]
            
            x = np.arange(len(genes))
            width = 0.25
            
            ax6.bar(x - width, betw_values, width, label='Betweenness', alpha=0.8)
            ax6.bar(x, close_values, width, label='Closeness', alpha=0.8)
            ax6.bar(x + width, eigen_values, width, label='Eigenvector', alpha=0.8)
            
            ax6.set_title('Centrality Measures (Top 10)', fontweight='bold')
            ax6.set_xlabel('Genes')
            ax6.set_ylabel('Centrality Score')
            ax6.set_xticks(x)
            ax6.set_xticklabels(genes, rotation=45, ha='right')
            ax6.legend()
        
    except:
        ax6.text(0.5, 0.5, 'Centrality measures\nnot available', 
                ha='center', va='center', transform=ax6.transAxes)
        ax6.set_title('Centrality Measures', fontweight='bold')
    
    plt.tight_layout()
    
    # Save the plot
    output_path = Path('GenePlots/plots')
    output_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path / f'{phenotype}_Analysis5_NetworkAnalysis.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(output_path / f'{phenotype}_Analysis5_NetworkAnalysis.pdf', 
                bbox_inches='tight')
    plt.close()
    
    # Save network data and statistics
    report_path = Path('GenePlots/reports')
    report_path.mkdir(parents=True, exist_ok=True)
    
    # Network basic statistics
    network_stats = {
        'Number_of_Nodes': G.number_of_nodes(),
        'Number_of_Edges': G.number_of_edges(),
        'Density': nx.density(G),
        'Average_Degree': np.mean([G.degree(node) for node in G.nodes()]),
        'Average_Clustering': nx.average_clustering(G),
        'Number_of_Connected_Components': nx.number_connected_components(G),
        'Largest_Component_Size': len(max(nx.connected_components(G), key=len)) if G.number_of_nodes() > 0 else 0
    }
    
    # Try to calculate additional metrics
    try:
        if nx.is_connected(G):
            network_stats['Average_Shortest_Path_Length'] = nx.average_shortest_path_length(G)
            network_stats['Diameter'] = nx.diameter(G)
        else:
            largest_cc = max(nx.connected_components(G), key=len)
            largest_subgraph = G.subgraph(largest_cc)
            network_stats['Average_Shortest_Path_Length'] = nx.average_shortest_path_length(largest_subgraph)
            network_stats['Diameter'] = nx.diameter(largest_subgraph)
    except:
        network_stats['Average_Shortest_Path_Length'] = 'N/A'
        network_stats['Diameter'] = 'N/A'
    
    stats_df = pd.DataFrame([network_stats])
    stats_df.to_csv(report_path / f'{phenotype}_Analysis5_NetworkStatistics.csv', index=False)
    
    # Node-level statistics
    node_stats = []
    try:
        betweenness = nx.betweenness_centrality(G)
        closeness = nx.closeness_centrality(G)
        eigenvector = nx.eigenvector_centrality(G, max_iter=1000)
    except:
        betweenness = {node: 0 for node in G.nodes()}
        closeness = {node: 0 for node in G.nodes()}
        eigenvector = {node: 0 for node in G.nodes()}
    
    for node in G.nodes():
        node_stats.append({
            'Gene': node,
            'Degree': G.degree(node),
            'Betweenness_Centrality': betweenness.get(node, 0),
            'Closeness_Centrality': closeness.get(node, 0),
            'Eigenvector_Centrality': eigenvector.get(node, 0),
            'Clustering_Coefficient': nx.clustering(G, node),
            'Database_Count': len(gene_databases.get(node, [])),
            'Databases': ';'.join(gene_databases.get(node, []))
        })
    
    node_stats_df = pd.DataFrame(node_stats)
    node_stats_df = node_stats_df.sort_values('Degree', ascending=False)
    node_stats_df.to_csv(report_path / f'{phenotype}_Analysis5_NodeStatistics.csv', index=False)
    
    # Save edge list
    edge_list = []
    for edge in G.edges(data=True):
        edge_list.append({
            'Gene1': edge[0],
            'Gene2': edge[1],
            'Weight': edge[2].get('weight', 1.0)
        })
    
    edge_df = pd.DataFrame(edge_list)
    edge_df.to_csv(report_path / f'{phenotype}_Analysis5_EdgeList.csv', index=False)
    
    print(f"✅ Network analysis completed")
    print(f"📊 Plot saved: GenePlots/plots/{phenotype}_Analysis5_NetworkAnalysis.png")
    print(f"📄 Reports saved: GenePlots/reports/{phenotype}_Analysis5_*.csv")
    print(f"🕸️ Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"📈 Average degree: {np.mean([G.degree(node) for node in G.nodes()]):.2f}")

def create_mock_network(gene_list, gene_databases):
    """
    Create a mock protein interaction network
    In a real implementation, this would query actual PPI databases
    """
    G = nx.Graph()
    
    # Add all genes as nodes
    for gene in gene_list:
        G.add_node(gene)
    
    # Create edges based on various criteria
    np.random.seed(42)  # For reproducible results
    
    # Strategy 1: Genes from same database are more likely to interact
    for gene1 in gene_list:
        for gene2 in gene_list:
            if gene1 < gene2:  # Avoid duplicates
                # Higher probability if genes share databases
                shared_dbs = set(gene_databases[gene1]) & set(gene_databases[gene2])
                
                if shared_dbs:
                    # Base probability increases with shared databases
                    prob = 0.1 + 0.05 * len(shared_dbs)
                else:
                    prob = 0.02  # Low baseline probability
                
                # Add some randomness based on gene names (for consistency)
                hash_val = hash(gene1 + gene2) % 1000 / 1000.0
                prob += hash_val * 0.05
                
                if np.random.random() < prob:
                    # Calculate edge weight based on confidence
                    weight = 0.3 + len(shared_dbs) * 0.2 + np.random.random() * 0.5
                    G.add_edge(gene1, gene2, weight=weight)
    
    # Strategy 2: Create some hub genes (highly connected genes)
    # Select top genes by database frequency as potential hubs
    gene_db_counts = [(gene, len(dbs)) for gene, dbs in gene_databases.items()]
    gene_db_counts.sort(key=lambda x: x[1], reverse=True)
    
    top_genes = [gene for gene, count in gene_db_counts[:min(5, len(gene_list)//10)]]
    
    for hub_gene in top_genes:
        # Connect hub to random other genes
        potential_partners = [g for g in gene_list if g != hub_gene and not G.has_edge(hub_gene, g)]
        n_connections = min(np.random.randint(3, 10), len(potential_partners))
        
        if potential_partners:
            partners = np.random.choice(potential_partners, size=n_connections, replace=False)
            for partner in partners:
                weight = 0.4 + np.random.random() * 0.4
                G.add_edge(hub_gene, partner, weight=weight)
    
    # Remove isolated nodes (genes with no interactions)
    isolated_nodes = [node for node in G.nodes() if G.degree(node) == 0]
    G.remove_nodes_from(isolated_nodes)
    
    print(f"🔗 Created network with {G.number_of_nodes()} connected genes and {G.number_of_edges()} interactions")
    
    return G

if __name__ == "__main__":
    print("This module should be run through main_analysis.py")
