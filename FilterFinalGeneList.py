#!/usr/bin/env python3
"""
NCBI Gene Validator
====================
Validates gene symbols against the official NCBI human gene database.
Resolves synonyms to current official symbols and adds Ensembl/Entrez IDs.

Usage:
    python FilterFinalGeneList.py <phenotype_or_csv> [output_valid] [output_invalid]

Examples:
    python FilterFinalGeneList.py migraine
    python FilterFinalGeneList.py migraine_all_genes.csv
    python FilterFinalGeneList.py migraine_all_genes.csv valid.csv invalid.csv
"""

import pandas as pd
import gzip
import requests
import sys
import os
import re
from pathlib import Path
from tqdm import tqdm

# ─────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────
GENE_INFO_URL  = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
GENE_INFO_FILE = "Homo_sapiens.gene_info.gz"
CACHE_FILE     = "human_genes_cache.pkl"
# Separate cache for the lookup dicts — avoids the attrs deepcopy crash
LOOKUP_CACHE   = "human_genes_lookups.pkl"


# ─────────────────────────────────────────────
# Resolve input: phenotype name or CSV path
# ─────────────────────────────────────────────
def resolve_input_file(arg: str) -> str:
    """
    Accept either a direct CSV path or a phenotype name.
    If the argument is not an existing file, search AllAnalysisGene/data/
    and AllPackagesGenes/ for a matching genes CSV.
    """
    # Direct file path
    if os.path.isfile(arg):
        return arg

    # Treat as phenotype — clean the same way R does
    clean = re.sub(r'[^a-zA-Z0-9_\-]', '_', arg)

    search_dirs = [
        Path("AllAnalysisGene") / "data",
        Path("AllPackagesGenes"),
        Path("AllAnalysisGene") / "reports",
        Path(".")
    ]

    # Priority order: all_genes > ALL_SOURCES > any _genes file
    candidates = []
    for d in search_dirs:
        if not d.exists():
            continue
        for f in d.glob(f"{clean}*genes*.csv"):
            candidates.append(f)
        for f in d.glob(f"{clean}*ALL_SOURCES*.csv"):
            candidates.append(f)

    if candidates:
        # Prefer the combined all_genes file
        for c in candidates:
            if "all_genes" in c.name.lower():
                print(f"   Auto-detected input file: {c}")
                return str(c)
        for c in candidates:
            if "all_sources" in c.name.lower():
                print(f"   Auto-detected input file: {c}")
                return str(c)
        print(f"   Auto-detected input file: {candidates[0]}")
        return str(candidates[0])

    print(f"❌ Could not find a CSV file for phenotype '{arg}'.")
    print(f"   Searched for '{clean}*genes*.csv' in:")
    for d in search_dirs:
        print(f"     {d}/")
    print("   Pass the full CSV path directly instead.")
    sys.exit(1)


# ─────────────────────────────────────────────
# Download and parse NCBI gene database
# ─────────────────────────────────────────────
def download_gene_database():
    """
    Download and parse NCBI human gene info.
    Returns (gene_db DataFrame, symbol_lookup dict, synonym_index dict).
    Lookup dicts are cached separately to avoid the attrs deepcopy crash.
    """
    # Try loading everything from cache
    if os.path.exists(CACHE_FILE) and os.path.exists(LOOKUP_CACHE):
        print(f"📦 Loading cached gene database from {CACHE_FILE}...")
        gene_db = pd.read_pickle(CACHE_FILE)
        print(f"✅ Loaded {len(gene_db):,} human genes from cache")

        print(f"📦 Loading cached lookup dicts from {LOOKUP_CACHE}...")
        lookups = pd.read_pickle(LOOKUP_CACHE)
        symbol_lookup  = lookups['symbol_lookup']
        synonym_index  = lookups['synonym_index']
        print(f"✅ Loaded lookup dicts "
              f"({len(symbol_lookup):,} symbols, {len(synonym_index):,} synonyms)\n")
        return gene_db, symbol_lookup, synonym_index

    # ── Download raw file if needed ───────────────────────────────────────
    if not os.path.exists(GENE_INFO_FILE):
        print("⬇️  Downloading NCBI gene database (~50 MB)...")
        response   = requests.get(GENE_INFO_URL, stream=True, timeout=120)
        total_size = int(response.headers.get('content-length', 0))

        with open(GENE_INFO_FILE, 'wb') as fh:
            with tqdm(total=total_size, unit='B', unit_scale=True,
                      desc="Downloading") as pbar:
                for chunk in response.iter_content(chunk_size=8192):
                    fh.write(chunk)
                    pbar.update(len(chunk))
        print("✅ Download complete\n")

    # ── Parse ─────────────────────────────────────────────────────────────
    print("📖 Parsing gene database...")
    gene_data = []

    with gzip.open(GENE_INFO_FILE, 'rt') as fh:
        fh.readline()  # skip header
        for line in tqdm(fh, desc="Reading genes", unit=" genes"):
            fields = line.strip().split('\t')
            if len(fields) < 15:
                continue

            gene_id  = fields[1]
            symbol   = fields[2]
            synonyms = fields[4].split('|') if fields[4] != '-' else []
            desc     = fields[8]

            ensembl_id = ''
            for ref in (fields[5].split('|') if fields[5] != '-' else []):
                if ref.startswith('Ensembl:'):
                    ensembl_id = ref.replace('Ensembl:', '')
                    break

            gene_data.append({
                'Symbol':          symbol.upper(),
                'Official_Symbol': symbol,
                'Gene_ID':         gene_id,
                'Description':     desc,
                'Ensembl_ID':      ensembl_id,
                'Synonyms':        '|'.join(synonyms).upper()
            })

    gene_db = pd.DataFrame(gene_data)

    # ── Build lookup dicts using vectorised pandas — no iterrows() ────────
    print("🔧 Building symbol lookup (vectorised)...")
    # symbol_lookup: uppercase symbol -> positional index
    symbol_lookup = dict(zip(gene_db['Symbol'], range(len(gene_db))))

    print("🔧 Building synonym index (vectorised)...")
    synonym_index: dict[str, int] = {}
    for idx, syn_str in enumerate(tqdm(gene_db['Synonyms'],
                                       desc="Indexing synonyms",
                                       unit=" genes")):
        if not syn_str:
            continue
        for syn in syn_str.split('|'):
            syn = syn.strip()
            if syn and syn not in synonym_index:
                synonym_index[syn] = idx

    # ── Save caches ───────────────────────────────────────────────────────
    # Save DataFrame WITHOUT attrs (keeps pickle small and fast)
    gene_db_clean = gene_db.copy()
    gene_db_clean.attrs = {}
    gene_db_clean.to_pickle(CACHE_FILE)

    pd.to_pickle({'symbol_lookup': symbol_lookup,
                  'synonym_index': synonym_index}, LOOKUP_CACHE)

    print(f"✅ Parsed {len(gene_db):,} human genes")
    print(f"💾 DataFrame cached:    {CACHE_FILE}")
    print(f"💾 Lookup dicts cached: {LOOKUP_CACHE}")
    print("   (Delete both files to re-download latest genes)\n")

    return gene_db, symbol_lookup, synonym_index


# ─────────────────────────────────────────────
# Validate a single gene symbol
# ─────────────────────────────────────────────
def validate_gene(symbol: str, gene_db: pd.DataFrame,
                  symbol_lookup: dict, synonym_index: dict):
    """
    Returns: (is_valid, official_symbol, gene_id, ensembl_id, description)
    """
    if not symbol or not str(symbol).strip():
        return False, None, None, None, None

    symbol_upper = str(symbol).strip().upper()

    idx = symbol_lookup.get(symbol_upper)
    if idx is not None:
        row = gene_db.iloc[idx]
        return (True, row['Official_Symbol'], row['Gene_ID'],
                row['Ensembl_ID'], row['Description'])

    idx = synonym_index.get(symbol_upper)
    if idx is not None:
        row = gene_db.iloc[idx]
        return (True, row['Official_Symbol'], row['Gene_ID'],
                row['Ensembl_ID'], row['Description'])

    return False, None, None, None, None


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main():
    print("=" * 80)
    print("NCBI Gene Validator - Local Database Mode")
    print("=" * 80 + "\n")

    if len(sys.argv) < 2:
        print("Usage: python FilterFinalGeneList.py <phenotype_or_csv> "
              "[output_valid] [output_invalid]")
        print()
        print("Examples:")
        print("  python FilterFinalGeneList.py migraine")
        print("  python FilterFinalGeneList.py migraine_all_genes.csv")
        print("  python FilterFinalGeneList.py migraine_all_genes.csv "
              "valid.csv invalid.csv")
        sys.exit(1)

    # ── Resolve input file ────────────────────────────────────────────────
    csv_file = resolve_input_file(sys.argv[1])

    stem           = Path(csv_file).stem
    output_valid   = sys.argv[2] if len(sys.argv) > 2 else f"{stem}_valid_genes.csv"
    output_invalid = sys.argv[3] if len(sys.argv) > 3 else f"{stem}_invalid_genes.csv"

    # ── Load gene database ────────────────────────────────────────────────
    gene_db, symbol_lookup, synonym_index = download_gene_database()

    # ── Load input CSV ────────────────────────────────────────────────────
    print(f"📄 Loading gene list from: {csv_file}")
    try:
        df = pd.read_csv(csv_file)
        print(f"✅ Loaded {len(df):,} rows\n")
    except FileNotFoundError:
        print(f"❌ Error: {csv_file} not found!")
        sys.exit(1)

    # Find gene column
    gene_columns = ['Gene', 'GeneSymbol', 'Gene_Symbol', 'gene_symbol', 'symbol']
    gene_col = next((c for c in gene_columns if c in df.columns), None)

    if gene_col is None:
        print("❌ No recognised gene column found.")
        print(f"   Available columns: {list(df.columns)}")
        sys.exit(1)

    print(f"✅ Using gene column: '{gene_col}'\n")

    # ── Sanity check ──────────────────────────────────────────────────────
    print("🧪 Testing with known gene symbols...")
    for gene in ['BRCA1', 'TP53', 'GLUT1', 'SLC2A1', 'NDUFA4L2', 'NOTCH3']:
        is_valid, official, gid, ensembl, desc = validate_gene(
            gene, gene_db, symbol_lookup, synonym_index)
        status       = "✓" if is_valid else "✗"
        official_str = official if official else 'NOT FOUND'
        desc_str     = desc[:40] if desc else 'N/A'
        print(f"  {status} {gene:15} → {official_str:15} | {desc_str}")

    print("\n" + "=" * 80)
    print(f"🚀 Validating {len(df):,} genes...\n")

    # ── Validate ──────────────────────────────────────────────────────────
    valid_rows   = []
    invalid_rows = []

    for _, row in tqdm(df.iterrows(), total=len(df),
                       desc="Validating", unit="gene"):
        gene = str(row[gene_col]).strip()

        (is_valid, official_symbol,
         gene_id, ensembl_id, description) = validate_gene(
            gene, gene_db, symbol_lookup, synonym_index)

        result = {
            'Gene_Symbol':     gene,
            'Is_Valid':        is_valid,
            'Official_Symbol': official_symbol if official_symbol else 'N/A',
            'Gene_Name':       description     if description     else 'N/A',
            'NCBI_Gene_ID':    gene_id         if gene_id         else 'N/A',
            'Ensembl_ID':      ensembl_id      if ensembl_id      else 'N/A'
        }

        for col in df.columns:
            if col not in result:
                result[col] = row[col]

        if is_valid:
            valid_rows.append(result)
        else:
            invalid_rows.append(result)

    print("\n" + "=" * 80)

    # ── Save ──────────────────────────────────────────────────────────────
    total = len(df)

    if valid_rows:
        valid_df = pd.DataFrame(valid_rows)
        valid_df.to_csv(output_valid, index=False)
        pct = 100 * len(valid_df) / total
        print(f"✅ Valid genes:   {len(valid_df):,} ({pct:.1f}%) → {output_valid}")
    else:
        print("❌ No valid genes found!")

    if invalid_rows:
        invalid_df = pd.DataFrame(invalid_rows)
        invalid_df.to_csv(output_invalid, index=False)
        pct = 100 * len(invalid_df) / total
        print(f"❌ Invalid genes: {len(invalid_df):,} ({pct:.1f}%) → {output_invalid}")

    if valid_rows:
        print(f"\n📋 Sample VALID genes:")
        for r in valid_rows[:5]:
            print(f"  ✓ {r['Gene_Symbol']:15} → "
                  f"{r['Official_Symbol']:15} | {r['Gene_Name'][:50]}")

    if invalid_rows:
        print(f"\n📋 Sample INVALID genes:")
        for r in invalid_rows[:10]:
            print(f"  ✗ {r['Gene_Symbol']}")

    print("\n" + "=" * 80)
    print("✨ DONE!")
    print(f"\n💡 Caches: {CACHE_FILE}  |  {LOOKUP_CACHE}")
    print("   (Delete both to re-download latest genes)")


if __name__ == "__main__":
    main()