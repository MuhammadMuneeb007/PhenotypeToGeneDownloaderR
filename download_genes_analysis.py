#!/usr/bin/env python3
"""
Gene Analysis Pipeline - Main Coordinator
==========================================
Runs comprehensive gene analysis for a given phenotype using data from AllPackagesGenes/
Saves all outputs to AllAnalysisGene/

Usage: python download_genes_analysis.py <phenotype>
Example: python download_genes_analysis.py migraine
"""

import sys
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from datetime import datetime

# Import analysis modules
sys.path.append('GenePlots')
try:
    import GenePlots.Analysis1_SourceComparison as analysis1
    import GenePlots.Analysis2_GeneFrequency as analysis2
    import GenePlots.Analysis3_DatabaseOverlap as analysis3
    import GenePlots.Analysis4_GeneSetEnrichment as analysis4
    import GenePlots.Analysis6_StatisticalSummary as analysis6

    run_source_comparison = analysis1.run_analysis
    run_gene_frequency    = analysis2.run_analysis
    run_database_overlap  = analysis3.run_analysis
    run_gene_enrichment   = analysis4.run_analysis
    run_statistical_summary = analysis6.run_analysis

except ImportError as e:
    print(f"⚠️  Warning: Could not import analysis modules: {e}")
    print("Make sure GenePlots/ directory contains all analysis modules")
    run_source_comparison   = None
    run_gene_frequency      = None
    run_database_overlap    = None
    run_gene_enrichment     = None
    run_statistical_summary = None


def clean_phenotype_for_filename(phenotype: str) -> str:
    """
    Reproduce the same filename-cleaning logic used by the R scripts:
      gsub("[^a-zA-Z0-9_-]", "_", phenotype)
    so that Python looks for the same filenames R actually created.
    """
    return re.sub(r'[^a-zA-Z0-9_\-]', '_', phenotype)


def load_gene_data(phenotype: str) -> dict:
    """
    Load gene data from AllPackagesGenes directory.
    Matches files with pattern: {clean_phenotype}_*_genes.csv
    Excludes summary files (ALL_SOURCES, SOURCES_SUMMARY).
    Removes duplicate gene entries within each source.

    Returns:
        dict: {database_name: DataFrame}
    """
    print(f"🔍 Loading gene data for phenotype: {phenotype}")

    data_dir = Path("AllPackagesGenes")
    if not data_dir.exists():
        print("❌ Error: AllPackagesGenes directory not found!")
        print(f"Please run the R pipeline first: Rscript download_genes.R {phenotype}")
        return {}

    # Use the same filename cleaning as R
    clean_pheno = clean_phenotype_for_filename(phenotype)
    print(f"   Filename prefix used: {clean_pheno}_")

    all_csv_files = list(data_dir.glob("*.csv"))
    csv_files = []

    for file_path in all_csv_files:
        stem = file_path.stem  # filename without extension, original case

        # Must start with the cleaned phenotype prefix
        if not stem.lower().startswith(f"{clean_pheno.lower()}_"):
            continue

        # Must end with _genes
        if not stem.lower().endswith("_genes"):
            continue

        # Exclude combined summary files
        stem_lower = stem.lower()
        if any(excl in stem_lower for excl in ["all_sources", "sources_summary"]):
            continue

        csv_files.append(file_path)

    if not csv_files:
        print(f"❌ No CSV files matching '{clean_pheno}_*_genes.csv' found.")
        print("Expected files like:")
        print(f"  - {clean_pheno}_pubmed_genes.csv")
        print(f"  - {clean_pheno}_omim_genes.csv")
        print(f"  - {clean_pheno}_opentargets_genes.csv")
        print(f"\nActual files found:")
        for f in sorted(all_csv_files):
            print(f"  - {f.name}")
        return {}

    print(f"📁 Found {len(csv_files)} gene data files:")

    gene_data    = {}
    gene_columns = ['Gene', 'GeneSymbol', 'Gene_Symbol', 'gene_symbol', 'symbol']

    for file_path in sorted(csv_files):
        try:
            df = pd.read_csv(file_path)

            # Extract database name: strip leading "{clean_pheno}_" and trailing "_genes"
            stem = file_path.stem
            # Remove phenotype prefix (case-insensitive)
            prefix = f"{clean_pheno}_"
            if stem.lower().startswith(prefix.lower()):
                remainder = stem[len(prefix):]
            else:
                remainder = stem

            # Remove trailing _genes
            if remainder.lower().endswith("_genes"):
                database_name = remainder[:-len("_genes")]
            else:
                database_name = remainder

            print(f"   📄 Processing {database_name}: {len(df)} initial records")

            if df.empty:
                print(f"   ⚠️  {database_name}: Empty file — skipping")
                continue

            # Find the gene column
            gene_col = None
            for col in gene_columns:
                if col in df.columns:
                    gene_col = col
                    break

            if gene_col is None:
                print(f"   ⚠️  {database_name}: No recognised gene column found")
                print(f"      Available columns: {list(df.columns)}")
                continue

            # Remove empty/null gene values
            df = df[df[gene_col].notna() & (df[gene_col] != '') & (df[gene_col] != 'NA')]
            original_count = len(df)

            # Remove duplicate genes
            df = df.drop_duplicates(subset=[gene_col], keep='first')
            deduplicated_count = len(df)

            if deduplicated_count < original_count:
                print(f"      🧹 Removed {original_count - deduplicated_count} duplicate genes")

            print(f"      ✅ Final: {deduplicated_count} records, {df[gene_col].nunique()} unique genes")
            gene_data[database_name] = df

        except Exception as e:
            print(f"   ❌ Error loading {file_path.name}: {e}")

    if not gene_data:
        print(f"❌ No valid gene data loaded for phenotype: {phenotype}")
        print("Make sure your files:")
        print(f"  1. Follow pattern: {clean_pheno}_<database>_genes.csv")
        print("  2. Contain a gene column (Gene, GeneSymbol, Gene_Symbol, etc.)")
        return {}

    total_records = sum(len(df) for df in gene_data.values())
    print(f"✅ Loaded data from {len(gene_data)} databases")
    print(f"📊 Total records after deduplication: {total_records}")
    print(f"🏷️  Databases extracted: {sorted(gene_data.keys())}")
    return gene_data


def setup_output_directories() -> Path:
    """Create output directory structure for AllAnalysisGene."""
    base_dir   = Path("AllAnalysisGene")
    for subdir in [base_dir, base_dir / "plots", base_dir / "reports", base_dir / "data"]:
        subdir.mkdir(parents=True, exist_ok=True)
    return base_dir


def run_comprehensive_analysis(phenotype: str):
    """Run the complete gene analysis pipeline."""
    print(f"\n🚀 Starting comprehensive gene analysis for: {phenotype.upper()}")
    print("=" * 60)

    output_dir = setup_output_directories()
    print(f"📁 Output directory: {output_dir.absolute()}")

    gene_data = load_gene_data(phenotype)
    if not gene_data:
        return

    plots_dir   = output_dir / "plots"
    reports_dir = output_dir / "reports"
    data_dir    = output_dir / "data"

    # ── Save combined data ────────────────────────────────────────────────
    print(f"\n💾 Saving combined data...")
    gene_columns = ['Gene', 'GeneSymbol', 'Gene_Symbol', 'gene_symbol', 'symbol']
    all_genes    = set()

    for db_name, df in gene_data.items():
        gene_col = next((c for c in gene_columns if c in df.columns), None)
        if gene_col:
            genes = {g for g in df[gene_col].dropna().unique()
                     if g != '' and pd.notna(g)}
            all_genes.update(genes)

    all_genes_df = pd.DataFrame({'Gene': sorted(all_genes)})
    all_genes_df['Phenotype']     = phenotype
    all_genes_df['Analysis_Date'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    all_genes_df.to_csv(data_dir / f'{phenotype}_all_genes.csv', index=False)

    # Summary statistics per database
    summary_rows = []
    for db_name, df in gene_data.items():
        gene_col = next((c for c in gene_columns if c in df.columns), None)
        unique_count = len(set(df[gene_col].dropna().unique())) if gene_col else 0
        summary_rows.append({
            'Database':     db_name,
            'Records':      len(df),
            'Unique_Genes': unique_count
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(data_dir / f'{phenotype}_summary_stats.csv', index=False)

    print(f"✅ Combined data saved: {len(all_genes)} unique genes across {len(gene_data)} databases")

    # ── Run analyses ──────────────────────────────────────────────────────
    # Ensure GenePlots output directories exist so analysis modules can write there
    original_plots_dir   = Path('GenePlots/plots')
    original_reports_dir = Path('GenePlots/reports')
    original_plots_dir.mkdir(parents=True, exist_ok=True)
    original_reports_dir.mkdir(parents=True, exist_ok=True)

    analyses = [
        ("Source Comparison",   run_source_comparison),
        ("Gene Frequency",      run_gene_frequency),
        ("Database Overlap",    run_database_overlap),
        ("Gene Set Enrichment", run_gene_enrichment),
        ("Statistical Summary", run_statistical_summary)
    ]

    completed_analyses = []

    for analysis_name, analysis_func in analyses:
        if analysis_func is None:
            print(f"⚠️  Skipping {analysis_name} — module not available")
            continue

        try:
            print(f"\n🔬 Running {analysis_name} Analysis...")
            analysis_func(gene_data, phenotype)

            # Move plots
            for ext in ("*.png", "*.pdf"):
                for plot_file in original_plots_dir.glob(f'{phenotype}_*{ext[1:]}'):
                    dest = plots_dir / plot_file.name
                    if plot_file.exists():
                        plot_file.rename(dest)

            # Move reports
            for ext in ("*.csv", "*.txt"):
                for report_file in original_reports_dir.glob(f'{phenotype}_*{ext[1:]}'):
                    dest = reports_dir / report_file.name
                    if report_file.exists():
                        report_file.rename(dest)

            completed_analyses.append(analysis_name)
            print(f"✅ {analysis_name} analysis completed")

        except Exception as e:
            print(f"❌ Error in {analysis_name} analysis: {e}")
            import traceback
            traceback.print_exc()

    # ── Final summary ─────────────────────────────────────────────────────
    print(f"\n🎉 Analysis Complete!")
    print("=" * 60)
    print(f"📊 Phenotype analysed:   {phenotype}")
    print(f"🗃️  Databases processed:  {len(gene_data)}")
    print(f"🧬 Total unique genes:   {len(all_genes)}")
    print(f"✅ Completed analyses:   {len(completed_analyses)}/{len(analyses)}")
    print(f"📁 All outputs saved to: {output_dir.absolute()}")

    if completed_analyses:
        print("\n📋 Completed analyses:")
        for a in completed_analyses:
            print(f"   ✓ {a}")

    print("\n📂 Output structure:")
    print("   📊 Plots:   AllAnalysisGene/plots/")
    print("   📄 Reports: AllAnalysisGene/reports/")
    print("   💾 Data:    AllAnalysisGene/data/")

    exec_summary = reports_dir / f'{phenotype}_ExecutiveSummary.txt'
    if exec_summary.exists():
        print(f"\n📋 Executive Summary: {exec_summary}")
        print("🔍 Review this file for key findings and recommendations")


def main():
    """Handle command line arguments and run analysis."""
    if len(sys.argv) != 2:
        print("Usage: python download_genes_analysis.py <phenotype>")
        print("Example: python download_genes_analysis.py migraine")
        sys.exit(1)

    phenotype = sys.argv[1].strip()

    if not phenotype:
        print("❌ Error: Please provide a valid phenotype name")
        sys.exit(1)

    print(f"🧬 Gene Analysis Pipeline")
    print(f"📅 Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"🎯 Phenotype: {phenotype}")

    if not Path("AllPackagesGenes").exists():
        print("\n❌ Error: AllPackagesGenes directory not found!")
        print("Please run the R pipeline first:")
        print(f"   Rscript download_genes.R {phenotype}")
        sys.exit(1)

    run_comprehensive_analysis(phenotype)


if __name__ == "__main__":
    main()