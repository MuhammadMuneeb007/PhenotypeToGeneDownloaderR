#!/usr/bin/env python3
"""
Analysis 1 - Pipeline Execution and Source Coverage
=====================================================
Prints tables to screen and saves as CSV.

Usage:
    python Analysis_1_Coverage.py
"""

import re
import sys
from pathlib import Path
import pandas as pd

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
OUT_DIR  = Path("AllAnalysisGene/Analysis1")
RUNTIME  = Path("runtime_summary.tsv")


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def count_genes(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        df = pd.read_csv(path)
        for col in ['Gene', 'GeneSymbol', 'Gene_Symbol', 'gene_symbol']:
            if col in df.columns:
                return int(df[col].dropna().nunique())
        return 0
    except Exception:
        return 0


def print_table(title: str, df: pd.DataFrame):
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)
    print(df.to_string(index=True))
    print("=" * 80)


# ─────────────────────────────────────────────────────────────────────────────
# Table 1 — per-phenotype × per-source gene counts
# ─────────────────────────────────────────────────────────────────────────────

def table1_count_matrix() -> pd.DataFrame:
    rows = []
    for pheno_label, pheno_clean in PHENOTYPES.items():
        row = {"Phenotype": pheno_label}
        for db_key, db_label in DATABASES.items():
            f = DATA_DIR / f"{pheno_clean}_{db_key}_genes.csv"
            row[db_label] = count_genes(f)
        rows.append(row)
    df = pd.DataFrame(rows).set_index("Phenotype")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Table 2 — per-database summary across all phenotypes
# ─────────────────────────────────────────────────────────────────────────────

def table2_source_summary(count_matrix: pd.DataFrame) -> pd.DataFrame:
    n = len(count_matrix)
    rows = []
    for db_label in count_matrix.columns:
        col     = count_matrix[db_label]
        success = col[col > 0]
        rows.append({
            "Database":            db_label,
            "Phenotypes_Success":  f"{len(success)}/{n}",
            "Success_Rate_%":      round(100 * len(success) / n, 1),
            "Mean_Genes":          round(success.mean(), 0) if len(success) else 0,
            "Median_Genes":        round(success.median(), 0) if len(success) else 0,
            "Max_Genes":           int(success.max()) if len(success) else 0,
            "Total_Genes_All":     int(col.sum()),
        })
    df = pd.DataFrame(rows).sort_values("Total_Genes_All", ascending=False)
    df = df.reset_index(drop=True)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Table 3 — combined + validated gene counts per phenotype
# ─────────────────────────────────────────────────────────────────────────────

def table3_validation_summary() -> pd.DataFrame:
    rows = []
    for pheno_label, pheno_clean in PHENOTYPES.items():
        combined  = DATA_DIR / f"{pheno_clean}_ALL_SOURCES_GENES.csv"
        valid_f   = DATA_DIR / f"{pheno_clean}_ALL_SOURCES_GENES_valid_genes.csv"
        invalid_f = DATA_DIR / f"{pheno_clean}_ALL_SOURCES_GENES_invalid_genes.csv"

        raw     = count_genes(combined)
        valid   = count_genes(valid_f)
        invalid = count_genes(invalid_f)

        if valid == 0 and invalid == 0:
            valid_str   = "not run"
            invalid_str = "not run"
            rate_str    = "not run"
        else:
            valid_str   = str(valid)
            invalid_str = str(invalid)
            rate_str    = f"{round(100 * valid / raw, 1)}%" if raw > 0 else "0%"

        rows.append({
            "Phenotype":          pheno_label,
            "Raw_Combined_Genes": raw,
            "Valid_Genes":        valid_str,
            "Invalid_Genes":      invalid_str,
            "Validation_Rate":    rate_str,
        })
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Table 4 — runtime per phenotype
# ─────────────────────────────────────────────────────────────────────────────

def table4_runtime() -> pd.DataFrame:
    if not RUNTIME.exists():
        print(f"   ⚠️  Runtime file not found: {RUNTIME}")
        return pd.DataFrame()

    try:
        df = pd.read_csv(RUNTIME, sep='\t')
        df.columns = df.columns.str.strip()

        def to_minutes(t):
            try:
                parts = str(t).strip().split(':')
                if len(parts) == 2:
                    return round(int(parts[0]) + int(parts[1]) / 60, 2)
                elif len(parts) == 3:
                    return round(int(parts[0]) * 60 + int(parts[1]) + int(parts[2]) / 60, 2)
            except Exception:
                return None

        if 'Time' in df.columns:
            df['Runtime_Minutes'] = df['Time'].apply(to_minutes)
        if 'MaxRSS_KB' in df.columns:
            df['MaxRSS_GB'] = df['MaxRSS_KB'].apply(
                lambda x: round(int(x) / 1e6, 3) if str(x).strip().isdigit() else None)

        return df
    except Exception as e:
        print(f"   ❌ Could not parse runtime file: {e}")
        return pd.DataFrame()


# ─────────────────────────────────────────────────────────────────────────────
# Table 5 — overall summary statistics
# ─────────────────────────────────────────────────────────────────────────────

def table5_overall(count_matrix: pd.DataFrame,
                   val_summary: pd.DataFrame) -> pd.DataFrame:
    total_cells   = count_matrix.size
    success_cells = (count_matrix > 0).sum().sum()
    sources_any   = int((count_matrix > 0).any(axis=0).sum())
    phenos_any    = int((count_matrix > 0).any(axis=1).sum())

    # Validated stats if available
    try:
        val_done = val_summary[val_summary['Valid_Genes'] != 'not run']
        total_valid = val_done['Valid_Genes'].astype(int).sum()
        total_raw   = val_done['Raw_Combined_Genes'].sum()
        val_rate    = f"{round(100 * total_valid / total_raw, 1)}%"
        val_genes   = str(int(total_valid))
    except Exception:
        val_rate  = "not run"
        val_genes = "not run"

    rows = [
        ("Phenotypes attempted",              len(PHENOTYPES)),
        ("Databases attempted",               len(DATABASES)),
        ("Phenotypes with ≥1 gene any source", f"{phenos_any}/{len(PHENOTYPES)}"),
        ("Databases with ≥1 gene any phenotype", f"{sources_any}/{len(DATABASES)}"),
        ("Overall source×phenotype success rate",
            f"{round(100 * success_cells / total_cells, 1)}%"),
        ("Total raw gene retrievals (sum)",
            f"{int(count_matrix.values.sum()):,}"),
        ("Overall symbol validation rate",    val_rate),
        ("Total validated unique genes",      val_genes),
    ]
    return pd.DataFrame(rows, columns=["Metric", "Value"])


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n🚀 Analysis 1 — Pipeline Execution and Source Coverage")
    print(f"   Data dir : {DATA_DIR.absolute()}")
    print(f"   Output   : {OUT_DIR.absolute()}")

    # Table 1
    print("\n⏳ Building Table 1 (count matrix)...")
    t1 = table1_count_matrix()
    print_table("TABLE 1 — Gene Counts per Phenotype × Database", t1)
    t1.to_csv(OUT_DIR / "Table1_Count_Matrix.csv")
    print(f"   💾 Saved Table1_Count_Matrix.csv")

    # Table 2
    print("\n⏳ Building Table 2 (source summary)...")
    t2 = table2_source_summary(t1)
    print_table("TABLE 2 — Per-Database Retrieval Summary (all phenotypes)", t2)
    t2.to_csv(OUT_DIR / "Table2_Source_Summary.csv", index=False)
    print(f"   💾 Saved Table2_Source_Summary.csv")

    # Table 3
    print("\n⏳ Building Table 3 (validation summary)...")
    t3 = table3_validation_summary()
    print_table("TABLE 3 — Combined Gene Counts and Validation per Phenotype", t3)
    t3.to_csv(OUT_DIR / "Table3_Validation_Summary.csv", index=False)
    print(f"   💾 Saved Table3_Validation_Summary.csv")

    # Table 4
    print("\n⏳ Building Table 4 (runtime)...")
    t4 = table4_runtime()
    if not t4.empty:
        print_table("TABLE 4 — Pipeline Runtime per Phenotype", t4)
        t4.to_csv(OUT_DIR / "Table4_Runtime.csv", index=False)
        print(f"   💾 Saved Table4_Runtime.csv")

    # Table 5
    print("\n⏳ Building Table 5 (overall summary)...")
    t5 = table5_overall(t1, t3)
    print_table("TABLE 5 — Overall Pipeline Statistics", t5)
    t5.to_csv(OUT_DIR / "Table5_Overall_Summary.csv", index=False)
    print(f"   💾 Saved Table5_Overall_Summary.csv")

    print("\n✅ Analysis 1 complete.")
    print(f"📁 All tables saved to: {OUT_DIR.absolute()}\n")


if __name__ == "__main__":
    main()