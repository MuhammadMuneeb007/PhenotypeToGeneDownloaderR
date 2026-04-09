#!/usr/bin/env python3
"""
Analysis 5+6 - Biological Enrichment and Robustness
=====================================================
Produces supplementary tables for Sections 5 and 6.

Section 5 - Biological Enrichment:
  Table S18: Pathway enrichment of final validated gene sets
  Table S19: Top enriched terms per phenotype
  Table S20: Cross-phenotype enrichment summary

Section 6 - Robustness:
  Table S21: Cross-phenotype comparison summary
  Table S22: Source contribution breakdown by evidence category
  Table S23: Pipeline reproducibility metrics

Usage:
    pip install gprofiler-official
    python Analysis_5_6_Enrichment_Robustness.py
"""

import re
import sys
import warnings
from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np

warnings.filterwarnings("ignore")

try:
    from gprofiler import GProfiler
    GPROFILER_AVAILABLE = True
except ImportError:
    GPROFILER_AVAILABLE = False
    print("⚠️  gprofiler-official not installed.")
    print("   Install with: pip install gprofiler-official")
    print("   Enrichment tables will be skipped.\n")

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

PHENOTYPES = {
    "asthma":                          "asthma",
    "blood pressure medication":       "blood_pressure_medication",
    "body mass index":                 "body_mass_index",
    "cholesterol lowering medication": "cholesterol_lowering_medication",
    "depression":                      "depression",
    "gastro-oesophageal reflux":       "gastro-oesophageal_reflux",
    "allergic rhinitis":               "allergic_rhinitis",
    "high cholesterol":                "high_cholesterol",
    "hypertension":                    "hypertension",
    "hypothyroidism":                  "hypothyroidism",
    "irritable bowel syndrome":        "irritable_bowel_syndrome",
    "migraine":                        "migraine",
    "osteoarthritis":                  "osteoarthritis",
}

DATABASES = {
    "clinvar":           "ClinVar",
    "gtex":              "GTEx",
    "gwasrapidd":        "GWAS Catalog",
    "hpo":               "HPO",
    "kegg":              "KEGG",
    "omim":              "OMIM",
    "opentargets":       "Open Targets",
    "pubmed":            "PubMed",
    "reactome_pathways": "Reactome",
    "string_db":         "STRING-DB",
    "uniprot":           "UniProt",
    "disgenet":          "DisGeNET",
    "gene_ontology":     "Gene Ontology",
}

# Evidence category groupings for Section 6
EVIDENCE_CATEGORIES = {
    "Clinical / Variant":   ["clinvar", "omim", "hpo"],
    "GWAS / Association":   ["gwasrapidd"],
    "Expression":           ["gtex"],
    "Pathway / Function":   ["kegg", "reactome_pathways",
                             "gene_ontology"],
    "Literature":           ["pubmed"],
    "Target / Drug":        ["opentargets"],
    "Protein / Network":    ["uniprot", "string_db", "disgenet"],
}

# g:Profiler sources to query
ENRICHMENT_SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "HP"]

# FDR threshold
FDR_THRESHOLD = 0.05

# Top N enriched terms to report per phenotype
TOP_N_TERMS = 10

DATA_DIR = Path("AllPackagesGenes")
OUT_DIR  = Path("AllAnalysisGene/Analysis5_6")


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def print_table(title: str, df: pd.DataFrame):
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)
    print(df.to_string(index=True))
    print("=" * 80)


def load_gene_set(path: Path) -> set:
    if not path.exists():
        return set()
    try:
        df = pd.read_csv(path)
        for col in ["Gene", "Official_Symbol",
                    "Gene_Symbol", "GeneSymbol"]:
            if col in df.columns:
                return set(
                    df[col].dropna().astype(str)
                    .str.strip().str.upper()
                )
        return set()
    except Exception:
        return set()


def load_valid_set(pheno_clean: str) -> set:
    valid_path = DATA_DIR / \
        f"{pheno_clean}_ALL_SOURCES_GENES_valid_genes.csv"
    if valid_path.exists():
        s = load_gene_set(valid_path)
        if s:
            return s
    combined = DATA_DIR / f"{pheno_clean}_ALL_SOURCES_GENES.csv"
    return load_gene_set(combined)


def load_ranked_genes(pheno_clean: str,
                      valid_set: set,
                      top_n: int = None) -> list:
    """
    Return genes ranked by number of supporting databases.
    Optionally limit to top_n.
    """
    gene_counts: dict[str, int] = {}
    for db_key in DATABASES:
        path  = DATA_DIR / f"{pheno_clean}_{db_key}_genes.csv"
        genes = load_gene_set(path)
        for g in genes & valid_set:
            gene_counts[g] = gene_counts.get(g, 0) + 1

    ranked = sorted(gene_counts.items(),
                    key=lambda x: (-x[1], x[0]))
    if top_n:
        ranked = ranked[:top_n]
    return [g for g, _ in ranked]


# ─────────────────────────────────────────────────────────────────────────────
# Section 5 — Biological Enrichment
# ─────────────────────────────────────────────────────────────────────────────

def run_enrichment(gene_list: list,
                   organism: str = "hsapiens") -> pd.DataFrame:
    """
    Run g:Profiler enrichment on a gene list.
    Returns a DataFrame of significant results.
    """
    if not GPROFILER_AVAILABLE or not gene_list:
        return pd.DataFrame()

    try:
        gp     = GProfiler(return_dataframe=True)
        result = gp.profile(
            organism   = organism,
            query      = gene_list,
            sources    = ENRICHMENT_SOURCES,
            user_threshold = FDR_THRESHOLD,
            significance_threshold_method = "fdr",
            no_iea     = True,
        )
        return result
    except Exception as e:
        print(f"      ⚠️  g:Profiler error: {e}")
        return pd.DataFrame()


def table_s18_enrichment_summary(enrich_results: dict) -> pd.DataFrame:
    """
    Per phenotype: number of significant terms per source category.
    """
    source_labels = {
        "GO:BP": "GO Biol. Process",
        "GO:MF": "GO Mol. Function",
        "GO:CC": "GO Cell. Component",
        "KEGG":  "KEGG",
        "REAC":  "Reactome",
        "HP":    "HPO",
    }

    rows = []
    for pheno_label in PHENOTYPES:
        df  = enrich_results.get(pheno_label, pd.DataFrame())
        row = {"Phenotype": pheno_label}

        if df.empty:
            for sl in source_labels.values():
                row[sl] = 0
            row["Total_Significant"] = 0
        else:
            total = 0
            for src_key, src_label in source_labels.items():
                n = len(df[df["source"] == src_key]) \
                    if "source" in df.columns else 0
                row[src_label] = n
                total         += n
            row["Total_Significant"] = total

        rows.append(row)

    return pd.DataFrame(rows).set_index("Phenotype")


def table_s19_top_terms(enrich_results: dict) -> pd.DataFrame:
    """
    Top N enriched terms per phenotype ranked by p-value.
    """
    rows = []
    for pheno_label in PHENOTYPES:
        df = enrich_results.get(pheno_label, pd.DataFrame())

        if df.empty:
            rows.append({
                "Phenotype":   pheno_label,
                "Rank":        1,
                "Source":      "—",
                "Term_ID":     "—",
                "Term_Name":   "No enrichment results",
                "p_value":     "—",
                "FDR":         "—",
                "Gene_Count":  "—",
            })
            continue

        # Sort by p-value
        df_sorted = df.sort_values("p_value").head(TOP_N_TERMS)

        for rank, (_, r) in enumerate(df_sorted.iterrows(), 1):
            rows.append({
                "Phenotype":   pheno_label,
                "Rank":        rank,
                "Source":      r.get("source",      "—"),
                "Term_ID":     r.get("native",       "—"),
                "Term_Name":   r.get("name",         "—"),
                "p_value":     f"{r.get('p_value', 0):.2e}",
                "FDR":         f"{r.get('p_value', 0):.2e}",
                "Gene_Count":  r.get("intersection_size", "—"),
            })

    return pd.DataFrame(rows)


def table_s20_cross_phenotype_enrichment(
        enrich_results: dict) -> pd.DataFrame:
    """
    Terms enriched in 3+ phenotypes simultaneously.
    """
    term_phenos: dict[str, list] = defaultdict(list)

    for pheno_label in PHENOTYPES:
        df = enrich_results.get(pheno_label, pd.DataFrame())
        if df.empty or "name" not in df.columns:
            continue
        for _, r in df.iterrows():
            term_name = r.get("name", "")
            if term_name:
                term_phenos[term_name].append(pheno_label)

    rows = []
    for term, phenos in term_phenos.items():
        if len(phenos) >= 3:
            rows.append({
                "Term":              term,
                "N_Phenotypes":      len(phenos),
                "Phenotypes":        "; ".join(sorted(phenos)),
            })

    if not rows:
        return pd.DataFrame(columns=["Term", "N_Phenotypes",
                                     "Phenotypes"])

    df = pd.DataFrame(rows).sort_values(
        "N_Phenotypes", ascending=False).reset_index(drop=True)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Section 6 — Robustness
# ─────────────────────────────────────────────────────────────────────────────

def table_s21_cross_phenotype_summary() -> pd.DataFrame:
    """
    Cross-phenotype comparison of key pipeline metrics.
    """
    # Load gold standard info from Analysis 4 if available
    gold_path = Path("AllAnalysisGene/Analysis4/TableS15_Recall_Precision.csv")
    gold_df   = pd.DataFrame()
    if gold_path.exists():
        try:
            gold_df = pd.read_csv(gold_path)
            if "Phenotype" in gold_df.columns:
                gold_df = gold_df.set_index("Phenotype")
        except Exception:
            pass

    # Load validation info from Analysis 1
    val_path = Path("AllAnalysisGene/Analysis1/Table3_Validation_Summary.csv")
    val_df   = pd.DataFrame()
    if val_path.exists():
        try:
            val_df = pd.read_csv(val_path)
            if "Phenotype" in val_df.columns:
                val_df = val_df.set_index("Phenotype")
        except Exception:
            pass

    rows = []
    for pheno_label, pheno_clean in PHENOTYPES.items():
        valid_set = load_valid_set(pheno_clean)
        ranked    = load_ranked_genes(pheno_clean, valid_set)

        # Count sources per phenotype
        n_sources = sum(
            1 for db_key in DATABASES
            if len(load_gene_set(
                DATA_DIR / f"{pheno_clean}_{db_key}_genes.csv")) > 0
        )

        # Genes supported by 3+ sources
        gene_counts: dict[str, int] = {}
        for db_key in DATABASES:
            path  = DATA_DIR / f"{pheno_clean}_{db_key}_genes.csv"
            genes = load_gene_set(path)
            for g in genes & valid_set:
                gene_counts[g] = gene_counts.get(g, 0) + 1

        n3plus = sum(1 for c in gene_counts.values() if c >= 3)
        total  = len(gene_counts)
        mean_s = round(
            sum(gene_counts.values()) / total, 2) if total else 0.0

        # Validation rate
        val_rate = "—"
        if not val_df.empty and pheno_label in val_df.index:
            val_rate = val_df.loc[pheno_label].get(
                "Validation_Rate", "—")

        # Recall
        recall = "—"
        if not gold_df.empty and pheno_label in gold_df.index:
            recall = gold_df.loc[pheno_label].get("Recall_%", "—")

        rows.append({
            "Phenotype":              pheno_label,
            "Sources_Success":        n_sources,
            "Total_Validated_Genes":  len(valid_set),
            "Validation_Rate":        val_rate,
            "Genes_3plus_Sources":    n3plus,
            "Genes_3plus_%":          round(100*n3plus/total, 1)
                                      if total else 0.0,
            "Mean_Sources_per_Gene":  mean_s,
            "Gold_Recall_%":          recall,
        })

    return pd.DataFrame(rows).set_index("Phenotype")


def table_s22_evidence_category_breakdown() -> pd.DataFrame:
    """
    For each phenotype: how many validated genes come from each
    evidence category (Clinical, GWAS, Expression, etc.)
    """
    rows = []
    for pheno_label, pheno_clean in PHENOTYPES.items():
        valid_set = load_valid_set(pheno_clean)
        row       = {"Phenotype": pheno_label}

        for cat_label, db_keys in EVIDENCE_CATEGORIES.items():
            cat_genes = set()
            for db_key in db_keys:
                path  = DATA_DIR / \
                    f"{pheno_clean}_{db_key}_genes.csv"
                cat_genes |= load_gene_set(path)
            if valid_set:
                cat_genes &= valid_set
            row[cat_label] = len(cat_genes)

        # Unique to single category
        all_cat_genes: dict[str, set] = {}
        for cat_label, db_keys in EVIDENCE_CATEGORIES.items():
            gs = set()
            for db_key in db_keys:
                path = DATA_DIR / \
                    f"{pheno_clean}_{db_key}_genes.csv"
                gs |= load_gene_set(path)
            if valid_set:
                gs &= valid_set
            all_cat_genes[cat_label] = gs

        exclusive = 0
        for cat_label, gs in all_cat_genes.items():
            others = set()
            for other_cat, other_gs in all_cat_genes.items():
                if other_cat != cat_label:
                    others |= other_gs
            exclusive += len(gs - others)

        row["Exclusive_to_1_category"] = exclusive
        rows.append(row)

    return pd.DataFrame(rows).set_index("Phenotype")


def table_s23_reproducibility() -> pd.DataFrame:
    """
    Report pipeline consistency metrics:
    - File presence check per phenotype/database
    - Consistency of gene counts with expected patterns
    """
    rows = []
    for pheno_label, pheno_clean in PHENOTYPES.items():
        expected_files = 0
        present_files  = 0
        empty_files    = 0

        for db_key in DATABASES:
            genes_path = DATA_DIR / \
                f"{pheno_clean}_{db_key}_genes.csv"
            expected_files += 1
            if genes_path.exists():
                present_files += 1
                gs = load_gene_set(genes_path)
                if len(gs) == 0:
                    empty_files += 1

        combined = DATA_DIR / \
            f"{pheno_clean}_ALL_SOURCES_GENES.csv"
        valid    = DATA_DIR / \
            f"{pheno_clean}_ALL_SOURCES_GENES_valid_genes.csv"

        rows.append({
            "Phenotype":            pheno_label,
            "DB_Files_Expected":    expected_files,
            "DB_Files_Present":     present_files,
            "DB_Files_Empty":       empty_files,
            "DB_Files_Missing":     expected_files - present_files,
            "Combined_File":        "Yes" if combined.exists() else "No",
            "Valid_File":           "Yes" if valid.exists()    else "No",
            "File_Completeness_%":  round(
                100 * present_files / expected_files, 1)
                if expected_files else 0.0,
        })

    return pd.DataFrame(rows).set_index("Phenotype")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n🚀 Analysis 5+6 — Enrichment and Robustness")
    print(f"   Data dir : {DATA_DIR.absolute()}")
    print(f"   Output   : {OUT_DIR.absolute()}")

    # ══════════════════════════════════════════════════════════════════════
    # SECTION 5 — ENRICHMENT
    # ══════════════════════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("  SECTION 5 — BIOLOGICAL ENRICHMENT")
    print("─" * 60)

    enrich_results = {}

    if GPROFILER_AVAILABLE:
        for pheno_label, pheno_clean in PHENOTYPES.items():
            print(f"\n⏳ Running enrichment for: {pheno_label}")
            valid_set = load_valid_set(pheno_clean)

            if not valid_set:
                print(f"   ⚠️  No validated genes — skipping")
                enrich_results[pheno_label] = pd.DataFrame()
                continue

            # Use top 500 genes ranked by source frequency
            gene_list = load_ranked_genes(
                pheno_clean, valid_set, top_n=500)

            print(f"   Querying g:Profiler with {len(gene_list)} genes...")
            result = run_enrichment(gene_list)

            if not result.empty:
                sig = result[result["p_value"] <= FDR_THRESHOLD] \
                    if "p_value" in result.columns else result
                print(f"   ✅ {len(sig)} significant terms (FDR≤{FDR_THRESHOLD})")
                enrich_results[pheno_label] = sig
            else:
                print(f"   ⚠️  No significant terms")
                enrich_results[pheno_label] = pd.DataFrame()

            # Save per-phenotype enrichment
            if not result.empty:
                result.to_csv(
                    OUT_DIR / f"{pheno_clean}_enrichment.csv",
                    index=False)
    else:
        print("\n⚠️  Skipping enrichment — gprofiler-official not available")
        for pheno_label in PHENOTYPES:
            enrich_results[pheno_label] = pd.DataFrame()

    # Table S18
    print("\n⏳ Building Table S18 (enrichment summary)...")
    t18 = table_s18_enrichment_summary(enrich_results)
    print_table(
        "TABLE S18 — Significant Enrichment Terms per Phenotype "
        f"(FDR ≤ {FDR_THRESHOLD})",
        t18
    )
    t18.to_csv(OUT_DIR / "TableS18_Enrichment_Summary.csv")
    print("   💾 Saved TableS18_Enrichment_Summary.csv")

    # Table S19
    print("\n⏳ Building Table S19 (top enriched terms per phenotype)...")
    t19 = table_s19_top_terms(enrich_results)
    print_table(
        f"TABLE S19 — Top {TOP_N_TERMS} Enriched Terms per Phenotype",
        t19
    )
    t19.to_csv(OUT_DIR / "TableS19_Top_Terms.csv", index=False)
    print("   💾 Saved TableS19_Top_Terms.csv")

    # Table S20
    print("\n⏳ Building Table S20 (cross-phenotype enrichment)...")
    t20 = table_s20_cross_phenotype_enrichment(enrich_results)
    print_table(
        "TABLE S20 — Terms Enriched in 3+ Phenotypes Simultaneously",
        t20
    )
    t20.to_csv(OUT_DIR / "TableS20_CrossPhenotype_Terms.csv",
               index=False)
    print("   💾 Saved TableS20_CrossPhenotype_Terms.csv")

    # ══════════════════════════════════════════════════════════════════════
    # SECTION 6 — ROBUSTNESS
    # ══════════════════════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("  SECTION 6 — ROBUSTNESS AND REPRODUCIBILITY")
    print("─" * 60)

    # Table S21
    print("\n⏳ Building Table S21 (cross-phenotype summary)...")
    t21 = table_s21_cross_phenotype_summary()
    print_table(
        "TABLE S21 — Cross-Phenotype Pipeline Performance Summary",
        t21
    )
    t21.to_csv(OUT_DIR / "TableS21_CrossPhenotype_Summary.csv")
    print("   💾 Saved TableS21_CrossPhenotype_Summary.csv")

    # Table S22
    print("\n⏳ Building Table S22 (evidence category breakdown)...")
    t22 = table_s22_evidence_category_breakdown()
    print_table(
        "TABLE S22 — Validated Gene Counts by Evidence Category "
        "per Phenotype",
        t22
    )
    t22.to_csv(OUT_DIR / "TableS22_Evidence_Category_Breakdown.csv")
    print("   💾 Saved TableS22_Evidence_Category_Breakdown.csv")

    # Table S23
    print("\n⏳ Building Table S23 (reproducibility / file completeness)...")
    t23 = table_s23_reproducibility()
    print_table(
        "TABLE S23 — Pipeline File Completeness and Reproducibility "
        "Metrics",
        t23
    )
    t23.to_csv(OUT_DIR / "TableS23_Reproducibility.csv")
    print("   💾 Saved TableS23_Reproducibility.csv")

    # ── Overall summary ───────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  FINAL OVERALL SUMMARY")
    print("=" * 80)

    total_sig = sum(
        len(df) for df in enrich_results.values()
        if not df.empty
    )
    phenos_enriched = sum(
        1 for df in enrich_results.values()
        if not df.empty
    )

    print(f"  Phenotypes with significant enrichment : "
          f"{phenos_enriched}/{len(PHENOTYPES)}")
    print(f"  Total significant enrichment terms     : {total_sig:,}")

    # Cross-phenotype stats from S21
    try:
        mean_sources = t21["Sources_Success"].mean()
        mean_valid   = t21["Total_Validated_Genes"].mean()
        mean_3plus   = t21["Genes_3plus_%"].mean()
        print(f"  Mean databases succeeding per phenotype: "
              f"{round(mean_sources,1)}")
        print(f"  Mean validated genes per phenotype     : "
              f"{int(mean_valid):,}")
        print(f"  Mean genes with 3+ source support      : "
              f"{round(mean_3plus,1)}%")
    except Exception:
        pass

    print("=" * 80)
    print(f"\n✅ Analysis 5+6 complete.")
    print(f"📁 All tables saved to: {OUT_DIR.absolute()}\n")


if __name__ == "__main__":
    main()