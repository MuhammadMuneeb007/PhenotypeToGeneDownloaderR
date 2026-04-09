#!/usr/bin/env python3
"""
Analysis 4 - Recovery of Known Phenotype-Associated Genes
==========================================================
Validates biological relevance by measuring how well the pipeline
recovers genes already known to be associated with each phenotype.

Gold standard approach:
  - Build a curated reference set per phenotype from three
    independent curated sources: HPO + ClinVar + OMIM
  - Treat genes supported by at least 2 of these 3 curated sources
    as the gold standard positive set
  - Evaluate recall, precision@k, and leave-one-source-out recovery

Tables produced:
  Table S14: Gold standard set sizes per phenotype
  Table S15: Recall and Precision@k for the full pipeline output
  Table S16: Leave-one-source-out recovery analysis
  Table S17: Top-ranked genes and their gold standard status

Usage:
    python Analysis_4_KnownGenes.py
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

# Curated sources used to build the gold standard
CURATED_SOURCES = {
    "hpo":     "HPO",
    "clinvar": "ClinVar",
    "omim":    "OMIM",
}

# Minimum number of curated sources a gene must appear in
# to be included in the gold standard
GOLD_MIN_SOURCES = 2

# All databases for leave-one-out
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
OUT_DIR  = Path("AllAnalysisGene/Analysis4")

K_VALUES = [10, 20, 50, 100]


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
        for col in ["Gene", "Gene_Symbol", "Official_Symbol",
                    "GeneSymbol", "gene_symbol"]:
            if col in df.columns:
                return set(
                    df[col].dropna().astype(str)
                    .str.strip().str.upper()
                )
        return set()
    except Exception:
        return set()


def load_valid_set(pheno_clean: str) -> set:
    """Load validated combined gene set for a phenotype."""
    valid_path = DATA_DIR / \
        f"{pheno_clean}_ALL_SOURCES_GENES_valid_genes.csv"
    if valid_path.exists():
        s = load_gene_set(valid_path)
        if s:
            return s
    combined = DATA_DIR / f"{pheno_clean}_ALL_SOURCES_GENES.csv"
    return load_gene_set(combined)


# ─────────────────────────────────────────────────────────────────────────────
# Build gold standard per phenotype
# ─────────────────────────────────────────────────────────────────────────────

def build_gold_standards() -> dict:
    """
    For each phenotype, build a gold standard from genes appearing
    in at least GOLD_MIN_SOURCES of the three curated databases
    (HPO, ClinVar, OMIM).

    Returns:
        {pheno_label: {
            'gold':     set of gold standard genes,
            'hpo':      set,
            'clinvar':  set,
            'omim':     set,
        }}
    """
    gold = {}
    for pheno_label, pheno_clean in PHENOTYPES.items():
        sets = {}
        for db_key, db_label in CURATED_SOURCES.items():
            path = DATA_DIR / f"{pheno_clean}_{db_key}_genes.csv"
            sets[db_key] = load_gene_set(path)

        # Count source support per gene
        gene_support: dict[str, int] = {}
        for db_key, gene_set in sets.items():
            for g in gene_set:
                gene_support[g] = gene_support.get(g, 0) + 1

        gold_set = {
            g for g, n in gene_support.items()
            if n >= GOLD_MIN_SOURCES
        }

        gold[pheno_label] = {
            "gold":    gold_set,
            "hpo":     sets.get("hpo",     set()),
            "clinvar": sets.get("clinvar", set()),
            "omim":    sets.get("omim",    set()),
        }

    return gold


# ─────────────────────────────────────────────────────────────────────────────
# Rank genes by number of supporting databases
# ─────────────────────────────────────────────────────────────────────────────

def rank_genes(pheno_clean: str, valid_set: set) -> list:
    """
    Rank genes in the combined valid set by how many databases
    retrieved them. Returns list of (gene, n_sources) tuples,
    highest first.
    """
    gene_counts: dict[str, int] = {}
    for db_key in DATABASES:
        path  = DATA_DIR / f"{pheno_clean}_{db_key}_genes.csv"
        genes = load_gene_set(path)
        # Only count genes that passed validation
        for g in genes & valid_set:
            gene_counts[g] = gene_counts.get(g, 0) + 1

    ranked = sorted(gene_counts.items(),
                    key=lambda x: (-x[1], x[0]))
    return ranked


# ─────────────────────────────────────────────────────────────────────────────
# Table S14 — gold standard set sizes
# ─────────────────────────────────────────────────────────────────────────────

def table_s14_gold_standard(gold: dict) -> pd.DataFrame:
    rows = []
    for pheno_label in PHENOTYPES:
        g     = gold[pheno_label]
        hpo   = g["hpo"]
        clv   = g["clinvar"]
        omim  = g["omim"]
        gs    = g["gold"]

        rows.append({
            "Phenotype":          pheno_label,
            "HPO_Genes":          len(hpo),
            "ClinVar_Genes":      len(clv),
            "OMIM_Genes":         len(omim),
            "Union_All_3":        len(hpo | clv | omim),
            f"Gold_Standard (≥{GOLD_MIN_SOURCES} sources)": len(gs),
            "HPO∩ClinVar":        len(hpo & clv),
            "HPO∩OMIM":           len(hpo & omim),
            "ClinVar∩OMIM":       len(clv & omim),
            "All_3_Agree":        len(hpo & clv & omim),
        })

    df = pd.DataFrame(rows).set_index("Phenotype")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Table S15 — recall and precision@k
# ─────────────────────────────────────────────────────────────────────────────

def table_s15_recall_precision(gold: dict) -> pd.DataFrame:
    rows = []
    for pheno_label, pheno_clean in PHENOTYPES.items():
        gs        = gold[pheno_label]["gold"]
        valid_set = load_valid_set(pheno_clean)
        ranked    = rank_genes(pheno_clean, valid_set)

        if not gs:
            row = {"Phenotype": pheno_label,
                   "Gold_Standard_Size": 0,
                   "Pipeline_Output_Size": len(valid_set),
                   "Recall_%": "—"}
            for k in K_VALUES:
                row[f"Recall@{k}_%"]    = "—"
                row[f"Precision@{k}_%"] = "—"
            rows.append(row)
            continue

        all_genes       = {g for g, _ in ranked}
        recovered       = gs & all_genes
        recall          = round(100 * len(recovered) / len(gs), 1) \
                          if gs else 0.0

        row = {
            "Phenotype":           pheno_label,
            "Gold_Standard_Size":  len(gs),
            "Pipeline_Output_Size":len(valid_set),
            "Recall_%":            recall,
        }

        ranked_genes = [g for g, _ in ranked]
        for k in K_VALUES:
            topk     = set(ranked_genes[:k])
            rec_k    = len(gs & topk)
            row[f"Recall@{k}_%"]    = round(
                100 * rec_k / len(gs), 1) if gs else 0.0
            row[f"Precision@{k}_%"] = round(
                100 * rec_k / k, 1) if k > 0 else 0.0

        rows.append(row)

    df = pd.DataFrame(rows).set_index("Phenotype")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Table S16 — leave-one-source-out recovery
# ─────────────────────────────────────────────────────────────────────────────

def table_s16_leave_one_out(gold: dict) -> pd.DataFrame:
    """
    For each database: remove it from the combined set and measure
    how many gold standard genes are still recovered by the remaining
    sources. Averaged across phenotypes where gold standard exists.
    """
    db_labels = list(DATABASES.values())
    rows      = []

    for db_key, db_label in DATABASES.items():
        total_gs       = 0
        total_recovered= 0
        n_phenotypes   = 0

        for pheno_label, pheno_clean in PHENOTYPES.items():
            gs = gold[pheno_label]["gold"]
            if not gs:
                continue

            # Build combined set excluding this database
            combined = set()
            for other_key, other_label in DATABASES.items():
                if other_key == db_key:
                    continue
                path  = DATA_DIR / f"{pheno_clean}_{other_key}_genes.csv"
                combined |= load_gene_set(path)

            recovered      = gs & combined
            total_gs       += len(gs)
            total_recovered+= len(recovered)
            n_phenotypes   += 1

        recall = round(100 * total_recovered / total_gs, 1) \
                 if total_gs > 0 else 0.0

        # How many gold standard genes does this DB uniquely contribute
        unique_gs = 0
        for pheno_label, pheno_clean in PHENOTYPES.items():
            gs = gold[pheno_label]["gold"]
            if not gs:
                continue
            focal  = load_gene_set(
                DATA_DIR / f"{pheno_clean}_{db_key}_genes.csv")
            others = set()
            for ok in DATABASES:
                if ok != db_key:
                    others |= load_gene_set(
                        DATA_DIR / f"{pheno_clean}_{ok}_genes.csv")
            unique_gs += len((gs & focal) - others)

        rows.append({
            "Database_Removed":          db_label,
            "Recall_Without_DB_%":       recall,
            "Unique_GoldStandard_Genes": unique_gs,
        })

    df = pd.DataFrame(rows).sort_values(
        "Recall_Without_DB_%").reset_index(drop=True)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Table S17 — top ranked genes with gold standard annotation
# ─────────────────────────────────────────────────────────────────────────────

def table_s17_top_genes_annotated(gold: dict,
                                   top_n: int = 15) -> pd.DataFrame:
    rows = []
    for pheno_label, pheno_clean in PHENOTYPES.items():
        gs        = gold[pheno_label]["gold"]
        valid_set = load_valid_set(pheno_clean)
        ranked    = rank_genes(pheno_clean, valid_set)

        for rank, (gene, n_sources) in enumerate(ranked[:top_n], 1):
            in_gold = gene in gs if gs else None
            rows.append({
                "Phenotype":        pheno_label,
                "Rank":             rank,
                "Gene":             gene,
                "N_Sources":        n_sources,
                "In_GoldStandard":  "Yes" if in_gold else (
                                    "No"  if in_gold is False else "—"),
            })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n🚀 Analysis 4 — Recovery of Known Phenotype-Associated Genes")
    print(f"   Data dir        : {DATA_DIR.absolute()}")
    print(f"   Output dir      : {OUT_DIR.absolute()}")
    print(f"   Gold standard   : genes in ≥{GOLD_MIN_SOURCES} of "
          f"{{HPO, ClinVar, OMIM}}")

    print("\n⏳ Building gold standard sets...")
    gold = build_gold_standards()
    print("   ✅ Done")

    # ── Table S14 ─────────────────────────────────────────────────────────
    print("\n⏳ Building Table S14 (gold standard sizes)...")
    t14 = table_s14_gold_standard(gold)
    print_table(
        f"TABLE S14 — Gold Standard Gene Set Sizes per Phenotype\n"
        f"  (genes in ≥{GOLD_MIN_SOURCES} of HPO, ClinVar, OMIM)",
        t14
    )
    t14.to_csv(OUT_DIR / "TableS14_Gold_Standard_Sizes.csv")
    print("   💾 Saved TableS14_Gold_Standard_Sizes.csv")

    # ── Table S15 ─────────────────────────────────────────────────────────
    print("\n⏳ Building Table S15 (recall and precision@k)...")
    t15 = table_s15_recall_precision(gold)
    print_table(
        "TABLE S15 — Recall and Precision@k for Pipeline Output\n"
        "  Genes ranked by number of supporting databases",
        t15
    )
    t15.to_csv(OUT_DIR / "TableS15_Recall_Precision.csv")
    print("   💾 Saved TableS15_Recall_Precision.csv")

    # ── Table S16 ─────────────────────────────────────────────────────────
    print("\n⏳ Building Table S16 (leave-one-source-out)...")
    t16 = table_s16_leave_one_out(gold)
    print_table(
        "TABLE S16 — Leave-One-Source-Out Gold Standard Recovery\n"
        "  Recall of gold standard genes when each database is removed",
        t16
    )
    t16.to_csv(OUT_DIR / "TableS16_LeaveOneOut.csv", index=False)
    print("   💾 Saved TableS16_LeaveOneOut.csv")

    # ── Table S17 ─────────────────────────────────────────────────────────
    print("\n⏳ Building Table S17 (top ranked genes annotated)...")
    t17 = table_s17_top_genes_annotated(gold, top_n=15)
    print_table(
        "TABLE S17 — Top 15 Ranked Genes per Phenotype with "
        "Gold Standard Annotation",
        t17
    )
    t17.to_csv(OUT_DIR / "TableS17_TopGenes_Annotated.csv", index=False)
    print("   💾 Saved TableS17_TopGenes_Annotated.csv")

    # ── Overall summary ───────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  OVERALL GOLD STANDARD RECOVERY SUMMARY")
    print("=" * 80)

    total_gs   = sum(len(v["gold"]) for v in gold.values())
    total_rec  = 0
    precisions = []

    for pheno_label, pheno_clean in PHENOTYPES.items():
        gs        = gold[pheno_label]["gold"]
        valid_set = load_valid_set(pheno_clean)
        ranked    = rank_genes(pheno_clean, valid_set)
        all_genes = {g for g, _ in ranked}

        if gs:
            total_rec  += len(gs & all_genes)
            top20       = {g for g, _ in ranked[:20]}
            p20         = round(100 * len(gs & top20) / 20, 1)
            precisions.append(p20)

    overall_recall = round(100 * total_rec / total_gs, 1) \
                     if total_gs > 0 else 0.0
    mean_p20       = round(np.mean(precisions), 1) \
                     if precisions else 0.0

    print(f"  Total gold standard genes (all phenotypes) : {total_gs:,}")
    print(f"  Total recovered by pipeline                : {total_rec:,}")
    print(f"  Overall recall                             : {overall_recall}%")
    print(f"  Mean Precision@20 across phenotypes        : {mean_p20}%")
    print("=" * 80)

    print(f"\n✅ Analysis 4 complete.")
    print(f"📁 All tables saved to: {OUT_DIR.absolute()}\n")


if __name__ == "__main__":
    main()