#!/usr/bin/env python3
"""
Analysis 2 - Gene Symbol Validation Quality
============================================
Produces supplementary tables for Section 2.

Tables produced:
  Table S6:  Per-database validation rate across all phenotypes
  Table S7:  Per-phenotype per-database valid/invalid breakdown
  Table S8:  Invalid symbol categories (artefact classification)
  Table S9:  Synonym rescue summary per phenotype

Usage:
    python Analysis_2_Validation.py
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
OUT_DIR  = Path("AllAnalysisGene/Analysis2")


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def print_table(title: str, df: pd.DataFrame):
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)
    print(df.to_string(index=True))
    print("=" * 80)


def load_csv_genes(path: Path, col: str = "Gene") -> list:
    """Load gene list from a CSV. Returns empty list if file missing."""
    if not path.exists():
        return []
    try:
        df = pd.read_csv(path)
        for c in [col, "Gene_Symbol", "GeneSymbol", "gene_symbol"]:
            if c in df.columns:
                return df[c].dropna().astype(str).str.strip().tolist()
        return []
    except Exception:
        return []


def load_validated(path: Path) -> pd.DataFrame:
    """Load a valid/invalid CSV produced by FilterFinalGeneList.py."""
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame()


def classify_invalid(symbol: str) -> str:
    """
    Classify an invalid gene symbol into a category:
      - OMIM_code      : short numeric or disorder codes e.g. A81, 100050
      - lncRNA         : Ensembl lncRNA IDs e.g. AC005154.6, LINC01234
      - Compound       : compound strings with semicolons or brackets
      - Numeric        : purely numeric
      - Short          : 1-2 characters
      - Other_artefact : everything else
    """
    s = str(symbol).strip()
    if re.match(r'^\d+$', s):
        return "Numeric"
    if len(s) <= 2:
        return "Short"
    if re.match(r'^AC\d+\.\d+$', s) or re.match(r'^AL\d+\.\d+$', s) or \
       re.match(r'^LINC\d+', s) or re.match(r'^MIR\d+', s):
        return "lncRNA / non-coding"
    if re.match(r'^[A-Z]{1,3}\d+$', s) and len(s) <= 5:
        return "OMIM disorder code"
    if ';' in s or '[' in s or '(' in s:
        return "Compound string"
    if re.match(r'^[A-Z][A-Z0-9\-]{0,3}$', s):
        return "OMIM disorder code"
    return "Other artefact"


# ─────────────────────────────────────────────────────────────────────────────
# Table S6 — per-database validation rate across all phenotypes
# ─────────────────────────────────────────────────────────────────────────────

def table_s6_per_database_validation() -> pd.DataFrame:
    """
    For each database and phenotype, load the per-database genes file
    and the per-database valid/invalid files if they exist.
    Falls back to estimating from the ALL_SOURCES valid set.
    """
    # Build a set of validated official symbols per phenotype
    # from the ALL_SOURCES valid files
    valid_sets = {}
    for pheno_label, pheno_clean in PHENOTYPES.items():
        valid_path = DATA_DIR / f"{pheno_clean}_ALL_SOURCES_GENES_valid_genes.csv"
        valid_df   = load_validated(valid_path)
        if not valid_df.empty:
            for col in ["Official_Symbol", "Gene_Symbol", "Gene"]:
                if col in valid_df.columns:
                    valid_sets[pheno_label] = set(
                        valid_df[col].dropna().astype(str).str.strip().str.upper()
                    )
                    break
        else:
            valid_sets[pheno_label] = None

    rows = []
    for db_key, db_label in DATABASES.items():
        db_row = {"Database": db_label}
        totals_raw   = 0
        totals_valid = 0
        n_phenotypes = 0

        for pheno_label, pheno_clean in PHENOTYPES.items():
            genes_path = DATA_DIR / f"{pheno_clean}_{db_key}_genes.csv"
            genes      = load_csv_genes(genes_path)

            if len(genes) == 0:
                db_row[pheno_label] = "—"
                continue

            n_phenotypes += 1
            totals_raw   += len(genes)

            vset = valid_sets.get(pheno_label)
            if vset is not None:
                n_valid = sum(1 for g in genes
                              if str(g).strip().upper() in vset)
                rate    = round(100 * n_valid / len(genes), 1)
                db_row[pheno_label] = f"{n_valid}/{len(genes)} ({rate}%)"
                totals_valid += n_valid
            else:
                db_row[pheno_label] = f"{len(genes)} (val. pending)"

        if totals_raw > 0 and totals_valid > 0:
            overall = round(100 * totals_valid / totals_raw, 1)
            db_row["Overall_Rate_%"] = overall
        else:
            db_row["Overall_Rate_%"] = "—"

        rows.append(db_row)

    df = pd.DataFrame(rows).set_index("Database")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Table S7 — compact per-phenotype per-database valid gene counts
# ─────────────────────────────────────────────────────────────────────────────

def table_s7_compact_valid_matrix() -> pd.DataFrame:
    """
    Matrix of validated gene counts per phenotype × database.
    Uses the ALL_SOURCES valid set to estimate per-database validated counts.
    """
    valid_sets = {}
    for pheno_label, pheno_clean in PHENOTYPES.items():
        valid_path = DATA_DIR / f"{pheno_clean}_ALL_SOURCES_GENES_valid_genes.csv"
        valid_df   = load_validated(valid_path)
        if not valid_df.empty:
            for col in ["Official_Symbol", "Gene_Symbol", "Gene"]:
                if col in valid_df.columns:
                    valid_sets[pheno_label] = set(
                        valid_df[col].dropna().astype(str).str.strip().str.upper()
                    )
                    break

    rows = []
    for pheno_label, pheno_clean in PHENOTYPES.items():
        row  = {"Phenotype": pheno_label}
        vset = valid_sets.get(pheno_label)

        for db_key, db_label in DATABASES.items():
            genes_path = DATA_DIR / f"{pheno_clean}_{db_key}_genes.csv"
            genes      = load_csv_genes(genes_path)

            if len(genes) == 0:
                row[db_label] = 0
            elif vset is not None:
                row[db_label] = sum(
                    1 for g in genes
                    if str(g).strip().upper() in vset
                )
            else:
                row[db_label] = len(genes)  # unvalidated — show raw

        rows.append(row)

    df = pd.DataFrame(rows).set_index("Phenotype")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Table S8 — invalid symbol categories
# ─────────────────────────────────────────────────────────────────────────────

def table_s8_invalid_categories() -> pd.DataFrame:
    """
    Classify all invalid symbols per phenotype and summarise by category.
    """
    category_counts = {pheno: {} for pheno in PHENOTYPES}

    for pheno_label, pheno_clean in PHENOTYPES.items():
        invalid_path = DATA_DIR / \
            f"{pheno_clean}_ALL_SOURCES_GENES_invalid_genes.csv"
        invalid_df = load_validated(invalid_path)

        if invalid_df.empty:
            continue

        sym_col = None
        for col in ["Gene_Symbol", "Gene", "GeneSymbol"]:
            if col in invalid_df.columns:
                sym_col = col
                break

        if sym_col is None:
            continue

        for sym in invalid_df[sym_col].dropna().astype(str):
            cat = classify_invalid(sym)
            category_counts[pheno_label][cat] = \
                category_counts[pheno_label].get(cat, 0) + 1

    all_cats = [
        "OMIM disorder code",
        "lncRNA / non-coding",
        "Compound string",
        "Numeric",
        "Short",
        "Other artefact"
    ]

    rows = []
    for pheno_label in PHENOTYPES:
        row = {"Phenotype": pheno_label}
        total = 0
        for cat in all_cats:
            n = category_counts[pheno_label].get(cat, 0)
            row[cat] = n
            total   += n
        row["Total_Invalid"] = total
        rows.append(row)

    df = pd.DataFrame(rows).set_index("Phenotype")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Table S9 — synonym rescue summary
# ─────────────────────────────────────────────────────────────────────────────

def table_s9_synonym_rescue() -> pd.DataFrame:
    """
    Count how many symbols were rescued via synonym mapping.
    A rescued symbol has Gene_Symbol != Official_Symbol in the valid file.
    """
    rows = []
    for pheno_label, pheno_clean in PHENOTYPES.items():
        valid_path = DATA_DIR / \
            f"{pheno_clean}_ALL_SOURCES_GENES_valid_genes.csv"
        valid_df = load_validated(valid_path)

        if valid_df.empty:
            rows.append({
                "Phenotype":       pheno_label,
                "Total_Valid":     "—",
                "Direct_Match":    "—",
                "Synonym_Rescued": "—",
                "Rescue_Rate_%":   "—",
            })
            continue

        # Check if both columns present
        has_input    = "Gene_Symbol"     in valid_df.columns
        has_official = "Official_Symbol" in valid_df.columns

        total = len(valid_df)

        if has_input and has_official:
            rescued = (
                valid_df["Gene_Symbol"].astype(str).str.strip().str.upper() !=
                valid_df["Official_Symbol"].astype(str).str.strip().str.upper()
            ).sum()
            direct  = total - rescued
            rate    = round(100 * rescued / total, 1) if total > 0 else 0.0
        else:
            rescued = "—"
            direct  = "—"
            rate    = "—"

        rows.append({
            "Phenotype":       pheno_label,
            "Total_Valid":     total,
            "Direct_Match":    direct,
            "Synonym_Rescued": rescued,
            "Rescue_Rate_%":   rate,
        })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n🚀 Analysis 2 — Gene Symbol Validation Quality")
    print(f"   Data dir : {DATA_DIR.absolute()}")
    print(f"   Output   : {OUT_DIR.absolute()}")

    # ── Table S6 ──────────────────────────────────────────────────────────
    print("\n⏳ Building Table S6 (per-database validation rates)...")
    t6 = table_s6_per_database_validation()
    print_table(
        "TABLE S6 — Per-Database Validation Rate per Phenotype\n"
        "  Format: valid/raw (rate%) or — if database returned no genes",
        t6
    )
    t6.to_csv(OUT_DIR / "TableS6_PerDatabase_ValidationRate.csv")
    print("   💾 Saved TableS6_PerDatabase_ValidationRate.csv")

    # ── Table S7 ──────────────────────────────────────────────────────────
    print("\n⏳ Building Table S7 (validated gene count matrix)...")
    t7 = table_s7_compact_valid_matrix()
    print_table(
        "TABLE S7 — Validated Gene Counts per Phenotype × Database\n"
        "  (symbols confirmed as official NCBI human gene entries)",
        t7
    )
    t7.to_csv(OUT_DIR / "TableS7_Validated_Count_Matrix.csv")
    print("   💾 Saved TableS7_Validated_Count_Matrix.csv")

    # ── Table S8 ──────────────────────────────────────────────────────────
    print("\n⏳ Building Table S8 (invalid symbol categories)...")
    t8 = table_s8_invalid_categories()
    print_table(
        "TABLE S8 — Classification of Invalid Gene Symbols per Phenotype",
        t8
    )
    t8.to_csv(OUT_DIR / "TableS8_Invalid_Categories.csv")
    print("   💾 Saved TableS8_Invalid_Categories.csv")

    # ── Table S9 ──────────────────────────────────────────────────────────
    print("\n⏳ Building Table S9 (synonym rescue summary)...")
    t9 = table_s9_synonym_rescue()
    print_table(
        "TABLE S9 — Synonym Rescue Summary per Phenotype\n"
        "  Rescued = input symbol differed from official symbol but mapped via synonym",
        t9
    )
    t9.to_csv(OUT_DIR / "TableS9_Synonym_Rescue.csv", index=False)
    print("   💾 Saved TableS9_Synonym_Rescue.csv")

    # ── Overall validation stats ──────────────────────────────────────────
    print("\n" + "=" * 80)
    print("  OVERALL VALIDATION SUMMARY")
    print("=" * 80)

    try:
        total_raw   = 0
        total_valid = 0
        total_inv   = 0

        for pheno_label, pheno_clean in PHENOTYPES.items():
            valid_path   = DATA_DIR / \
                f"{pheno_clean}_ALL_SOURCES_GENES_valid_genes.csv"
            invalid_path = DATA_DIR / \
                f"{pheno_clean}_ALL_SOURCES_GENES_invalid_genes.csv"

            nv = len(load_validated(valid_path))
            ni = len(load_validated(invalid_path))

            if nv > 0 or ni > 0:
                total_valid += nv
                total_inv   += ni
                total_raw   += nv + ni

        if total_raw > 0:
            print(f"  Total raw symbols (validated phenotypes) : {total_raw:,}")
            print(f"  Total valid                              : {total_valid:,} "
                  f"({round(100*total_valid/total_raw,1)}%)")
            print(f"  Total invalid                            : {total_inv:,} "
                  f"({round(100*total_inv/total_raw,1)}%)")
    except Exception as e:
        print(f"  Could not compute overall stats: {e}")

    print("=" * 80)
    print(f"\n✅ Analysis 2 complete.")
    print(f"📁 All tables saved to: {OUT_DIR.absolute()}\n")


if __name__ == "__main__":
    main()