#!/usr/bin/env python3
"""
Analysis 3 - Cross-Source Overlap and Complementarity
======================================================
Produces supplementary tables for Section 3.

Tables produced:
  Table S10: Pairwise Jaccard similarity matrix (averaged across phenotypes)
  Table S11: Gene frequency distribution (how many sources support each gene)
  Table S12: Top recurrent genes per phenotype (supported by most sources)
  Table S13: Source complementarity — unique contribution per database

Usage:
    python Analysis_3_Overlap.py
"""

import re
import sys
from pathlib import Path
from itertools import combinations
import pandas as pd
import numpy as np

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

DATA_DIR = Path("AllPackagesGenes")
OUT_DIR  = Path("AllAnalysisGene/Analysis3")

TOP_N_GENES = 20   # top recurrent genes to report per phenotype


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
    """Load a set of uppercase gene symbols from a genes-only CSV."""
    if not path.exists():
        return set()
    try:
        df = pd.read_csv(path)
        for col in ["Gene", "Gene_Symbol", "GeneSymbol", "gene_symbol"]:
            if col in df.columns:
                return set(
                    df[col].dropna().astype(str).str.strip().str.upper()
                )
        return set()
    except Exception:
        return set()


def load_valid_gene_set(pheno_clean: str) -> set:
    """
    Load validated gene set for a phenotype.
    Prefers the valid_genes file; falls back to ALL_SOURCES_GENES.
    """
    valid_path = DATA_DIR / \
        f"{pheno_clean}_ALL_SOURCES_GENES_valid_genes.csv"
    if valid_path.exists():
        try:
            df = pd_read(valid_path)
            for col in ["Official_Symbol", "Gene_Symbol", "Gene"]:
                if col in df.columns:
                    return set(
                        df[col].dropna().astype(str).str.strip().str.upper()
                    )
        except Exception:
            pass

    combined_path = DATA_DIR / f"{pheno_clean}_ALL_SOURCES_GENES.csv"
    return load_gene_set(combined_path)


def pd_read(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame()


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return round(inter / union, 4) if union > 0 else 0.0


# ─────────────────────────────────────────────────────────────────────────────
# Load all per-database gene sets for all phenotypes
# ─────────────────────────────────────────────────────────────────────────────

def load_all_db_sets() -> dict:
    """
    Returns:
        {pheno_label: {db_label: set_of_genes}}
    Only includes validated genes where available.
    """
    # Build validated sets per phenotype
    valid_sets = {}
    for pheno_label, pheno_clean in PHENOTYPES.items():
        valid_sets[pheno_label] = load_valid_gene_set(pheno_clean)

    data = {}
    for pheno_label, pheno_clean in PHENOTYPES.items():
        data[pheno_label] = {}
        vset = valid_sets[pheno_label]

        for db_key, db_label in DATABASES.items():
            genes_path = DATA_DIR / f"{pheno_clean}_{db_key}_genes.csv"
            raw        = load_gene_set(genes_path)

            if len(vset) > 0 and len(raw) > 0:
                # Intersect with validated set
                data[pheno_label][db_label] = raw & vset
            else:
                data[pheno_label][db_label] = raw

    return data


# ─────────────────────────────────────────────────────────────────────────────
# Table S10 — pairwise Jaccard similarity matrix
# ─────────────────────────────────────────────────────────────────────────────

def table_s10_jaccard(all_db_sets: dict) -> pd.DataFrame:
    """
    Average Jaccard similarity between every pair of databases
    across all phenotypes where both databases returned genes.
    """
    db_labels = list(DATABASES.values())
    matrix    = pd.DataFrame(0.0, index=db_labels, columns=db_labels)
    counts    = pd.DataFrame(0,   index=db_labels, columns=db_labels)

    for pheno_label, db_sets in all_db_sets.items():
        for db_a, db_b in combinations(db_labels, 2):
            a = db_sets.get(db_a, set())
            b = db_sets.get(db_b, set())
            if len(a) > 0 and len(b) > 0:
                j = jaccard(a, b)
                matrix.loc[db_a, db_b] += j
                matrix.loc[db_b, db_a] += j
                counts.loc[db_a, db_b] += 1
                counts.loc[db_b, db_a] += 1

    # Diagonal = 1
    for db in db_labels:
        matrix.loc[db, db] = 1.0
        counts.loc[db, db] = 1

    # Average
    avg = matrix.copy()
    for db_a in db_labels:
        for db_b in db_labels:
            if db_a != db_b and counts.loc[db_a, db_b] > 0:
                avg.loc[db_a, db_b] = round(
                    matrix.loc[db_a, db_b] / counts.loc[db_a, db_b], 3)
            elif db_a != db_b:
                avg.loc[db_a, db_b] = 0.0

    return avg


# ─────────────────────────────────────────────────────────────────────────────
# Table S11 — gene frequency distribution
# ─────────────────────────────────────────────────────────────────────────────

def table_s11_frequency_distribution(all_db_sets: dict) -> pd.DataFrame:
    """
    For each phenotype: how many genes are supported by exactly
    1, 2, 3, 4, 5+ databases.
    """
    rows = []
    for pheno_label, db_sets in all_db_sets.items():
        # Count how many databases support each gene
        gene_counts: dict[str, int] = {}
        for db_label, genes in db_sets.items():
            for g in genes:
                gene_counts[g] = gene_counts.get(g, 0) + 1

        if not gene_counts:
            rows.append({
                "Phenotype": pheno_label,
                "1 source":  0,
                "2 sources": 0,
                "3 sources": 0,
                "4 sources": 0,
                "5+ sources":0,
                "Total_Unique": 0,
                "Mean_Sources": 0.0,
            })
            continue

        counts = list(gene_counts.values())
        rows.append({
            "Phenotype":   pheno_label,
            "1 source":    sum(1 for c in counts if c == 1),
            "2 sources":   sum(1 for c in counts if c == 2),
            "3 sources":   sum(1 for c in counts if c == 3),
            "4 sources":   sum(1 for c in counts if c == 4),
            "5+ sources":  sum(1 for c in counts if c >= 5),
            "Total_Unique":len(gene_counts),
            "Mean_Sources":round(sum(counts) / len(counts), 2),
        })

    df = pd.DataFrame(rows).set_index("Phenotype")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Table S12 — top recurrent genes per phenotype
# ─────────────────────────────────────────────────────────────────────────────

def table_s12_top_recurrent(all_db_sets: dict) -> pd.DataFrame:
    """
    Top N genes by number of supporting databases, per phenotype.
    Columns: Phenotype, Gene, N_Sources, Supporting_Databases
    """
    rows = []
    for pheno_label, db_sets in all_db_sets.items():
        gene_support: dict[str, list] = {}
        for db_label, genes in db_sets.items():
            for g in genes:
                if g not in gene_support:
                    gene_support[g] = []
                gene_support[g].append(db_label)

        if not gene_support:
            continue

        sorted_genes = sorted(
            gene_support.items(),
            key=lambda x: (-len(x[1]), x[0])
        )

        for gene, sources in sorted_genes[:TOP_N_GENES]:
            rows.append({
                "Phenotype":           pheno_label,
                "Gene":                gene,
                "N_Sources":           len(sources),
                "Supporting_Databases":"; ".join(sorted(sources)),
            })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Table S13 — unique contribution per database
# ─────────────────────────────────────────────────────────────────────────────

def table_s13_unique_contribution(all_db_sets: dict) -> pd.DataFrame:
    """
    For each database: how many genes does it contribute that are NOT
    found in any other database, summed across all phenotypes.
    Also reports total genes and fraction that are unique to that source.
    """
    db_labels = list(DATABASES.values())
    rows = []

    for db_label in db_labels:
        total_genes  = 0
        unique_genes = 0

        for pheno_label, db_sets in all_db_sets.items():
            focal = db_sets.get(db_label, set())
            if not focal:
                continue

            others = set()
            for other_db, other_genes in db_sets.items():
                if other_db != db_label:
                    others |= other_genes

            total_genes  += len(focal)
            unique_genes += len(focal - others)

        unique_rate = round(100 * unique_genes / total_genes, 1) \
            if total_genes > 0 else 0.0

        rows.append({
            "Database":           db_label,
            "Total_Genes_All":    total_genes,
            "Unique_Only_Genes":  unique_genes,
            "Unique_Rate_%":      unique_rate,
        })

    df = pd.DataFrame(rows).sort_values(
        "Unique_Only_Genes", ascending=False
    ).reset_index(drop=True)

    return df


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n🚀 Analysis 3 — Cross-Source Overlap and Complementarity")
    print(f"   Data dir : {DATA_DIR.absolute()}")
    print(f"   Output   : {OUT_DIR.absolute()}")

    print("\n⏳ Loading all per-database gene sets...")
    all_db_sets = load_all_db_sets()
    print("   ✅ Loaded")

    # ── Table S10 ─────────────────────────────────────────────────────────
    print("\n⏳ Building Table S10 (Jaccard similarity matrix)...")
    t10 = table_s10_jaccard(all_db_sets)
    print_table(
        "TABLE S10 — Average Pairwise Jaccard Similarity between Databases\n"
        "  (averaged across phenotypes where both databases returned genes)",
        t10
    )
    t10.to_csv(OUT_DIR / "TableS10_Jaccard_Matrix.csv")
    print("   💾 Saved TableS10_Jaccard_Matrix.csv")

    # ── Table S11 ─────────────────────────────────────────────────────────
    print("\n⏳ Building Table S11 (gene frequency distribution)...")
    t11 = table_s11_frequency_distribution(all_db_sets)
    print_table(
        "TABLE S11 — Gene Frequency Distribution\n"
        "  (number of genes supported by exactly N databases per phenotype)",
        t11
    )
    t11.to_csv(OUT_DIR / "TableS11_Frequency_Distribution.csv")
    print("   💾 Saved TableS11_Frequency_Distribution.csv")

    # ── Table S12 ─────────────────────────────────────────────────────────
    print("\n⏳ Building Table S12 (top recurrent genes per phenotype)...")
    t12 = table_s12_top_recurrent(all_db_sets)
    print_table(
        f"TABLE S12 — Top {TOP_N_GENES} Recurrent Genes per Phenotype\n"
        "  (ranked by number of supporting databases)",
        t12
    )
    t12.to_csv(OUT_DIR / "TableS12_Top_Recurrent_Genes.csv", index=False)
    print("   💾 Saved TableS12_Top_Recurrent_Genes.csv")

    # ── Table S13 ─────────────────────────────────────────────────────────
    print("\n⏳ Building Table S13 (unique contribution per database)...")
    t13 = table_s13_unique_contribution(all_db_sets)
    print_table(
        "TABLE S13 — Unique Gene Contribution per Database\n"
        "  (genes found exclusively in that database, not in any other)",
        t13
    )
    t13.to_csv(OUT_DIR / "TableS13_Unique_Contribution.csv", index=False)
    print("   💾 Saved TableS13_Unique_Contribution.csv")

    # ── Summary stats ──────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  CROSS-SOURCE OVERLAP SUMMARY")
    print("=" * 80)

    # Average Jaccard excluding diagonal and zeros
    jac_vals = []
    db_labels = list(DATABASES.values())
    for i, db_a in enumerate(db_labels):
        for j, db_b in enumerate(db_labels):
            if i < j:
                val = t10.loc[db_a, db_b]
                if val > 0:
                    jac_vals.append(val)

    if jac_vals:
        print(f"  Mean pairwise Jaccard (non-zero pairs) : "
              f"{round(np.mean(jac_vals), 3)}")
        print(f"  Max pairwise Jaccard                   : "
              f"{round(max(jac_vals), 3)}")
        print(f"  Min pairwise Jaccard (non-zero)        : "
              f"{round(min(jac_vals), 3)}")

    # Genes supported by 3+ sources
    for pheno_label, db_sets in all_db_sets.items():
        gene_counts: dict[str, int] = {}
        for db_label, genes in db_sets.items():
            for g in genes:
                gene_counts[g] = gene_counts.get(g, 0) + 1
        n3plus = sum(1 for c in gene_counts.values() if c >= 3)
        total  = len(gene_counts)
        if total > 0:
            print(f"  {pheno_label:35s}: "
                  f"{n3plus}/{total} genes supported by 3+ sources "
                  f"({round(100*n3plus/total,1)}%)")

    print("=" * 80)
    print(f"\n✅ Analysis 3 complete.")
    print(f"📁 All tables saved to: {OUT_DIR.absolute()}\n")


if __name__ == "__main__":
    main()