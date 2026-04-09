#!/bin/bash
#SBATCH --job-name=Phenotype2Gene
#SBATCH --nodes=1
#SBATCH --partition=ascher
#SBATCH --time=24:00:00
#SBATCH --output=download_genes.%A_%a.out
#SBATCH --error=download_genes.%A_%a.err
#SBATCH --array=1-13
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1

set -euo pipefail

PHENO_FILE="PhenotypeNames.txt"

if [[ ! -f "$PHENO_FILE" ]]; then
    echo "Error: $PHENO_FILE not found"
    exit 1
fi

phenotype=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PHENO_FILE")

if [[ -z "${phenotype:-}" ]]; then
    echo "Error: No phenotype found for task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

safe_name=$(echo "$phenotype" | tr ' /()' '_' | tr -s '_' | sed 's/^_//; s/_$//')
outfile="runtime_${safe_name}.tsv"

echo -e "TaskID\tPhenotype\tElapsed\tMaxRSS_KB\tExitStatus" > "$outfile"

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Phenotype: $phenotype"
echo "Start: $(date)"

/usr/bin/time -f "${SLURM_ARRAY_TASK_ID}\t${phenotype}\t%E\t%M\t%x" \
    -a -o "$outfile" \
    Rscript download_genes.R "$phenotype"

echo "End: $(date)"