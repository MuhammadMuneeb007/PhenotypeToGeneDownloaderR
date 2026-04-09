#!/usr/bin/env Rscript
# KEGG phenotype -> matched pathway titles -> genes
# No hardcoded disease expansions
# No separate batched KEGG gene symbol conversion
#
# Usage:
#   Rscript kegg.R migraine
#   Rscript kegg.R cancer
#   Rscript kegg.R "breast cancer"

if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  stop("Package 'KEGGREST' is not installed. Install it first with: BiocManager::install('KEGGREST')")
}
suppressPackageStartupMessages(library(KEGGREST))

empty_result <- function() {
  data.frame(
    Phenotype = character(),
    Pathway_ID = character(),
    Pathway_Name = character(),
    KEGG_Gene_ID = character(),
    Gene = character(),
    Source = character(),
    stringsAsFactors = FALSE
  )
}

safe_trim <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  x <- as.character(x[1])
  if (is.na(x)) return(NA_character_)
  trimws(x)
}

is_blank <- function(x) {
  if (is.null(x) || length(x) == 0) return(TRUE)
  x <- as.character(x[1])
  if (is.na(x)) return(TRUE)
  !nzchar(trimws(x))
}

clean_filename <- function(x) {
  x <- gsub("[^[:alnum:]_ -]", "", x)
  x <- gsub("\\s+", "_", trimws(x))
  if (!nzchar(x)) x <- "phenotype"
  x
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

split_words <- function(x) {
  x <- tolower(trimws(x))
  words <- unlist(strsplit(x, "\\s+"))
  words <- words[!is.na(words) & nzchar(words)]
  unique(words)
}

score_pathway_match <- function(phenotype, pathway_name) {
  phenotype <- tolower(trimws(phenotype))
  pathway_name <- tolower(trimws(pathway_name))

  if (!nzchar(phenotype) || !nzchar(pathway_name)) return(0)

  score <- 0

  # full phrase match gets highest weight
  if (grepl(phenotype, pathway_name, fixed = TRUE)) {
    score <- score + 100
  }

  # word overlap
  pwords <- split_words(phenotype)
  if (length(pwords) > 0) {
    for (w in pwords) {
      pattern <- paste0("\\b", escape_regex(w), "\\b")
      if (grepl(pattern, pathway_name, perl = TRUE)) {
        score <- score + 1
      }
    }
  }

  score
}

find_matching_pathways <- function(phenotype, min_score = 1) {
  cat("🔍 Searching KEGG human pathway titles for:", phenotype, "\n")
  cat("   Fetching all human KEGG pathways...\n")

  all_pathways <- tryCatch(
    KEGGREST::keggList("pathway", "hsa"),
    error = function(e) {
      cat("❌ Failed to fetch KEGG pathway list:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(all_pathways) || length(all_pathways) == 0) {
    return(data.frame(
      Pathway_ID = character(),
      Pathway_Name = character(),
      Score = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  pathway_ids <- names(all_pathways)
  pathway_names <- as.character(unname(all_pathways))
  pathway_names[is.na(pathway_names)] <- ""

  scores <- vapply(
    pathway_names,
    function(x) score_pathway_match(phenotype, x),
    numeric(1)
  )

  matched <- data.frame(
    Pathway_ID = pathway_ids,
    Pathway_Name = pathway_names,
    Score = scores,
    stringsAsFactors = FALSE
  )

  matched <- matched[matched$Score >= min_score, , drop = FALSE]

  if (nrow(matched) == 0) {
    cat("   No pathway titles matched the phenotype text.\n")
    return(matched)
  }

  matched <- matched[order(-matched$Score, matched$Pathway_Name), , drop = FALSE]

  cat("   Matched", nrow(matched), "pathways\n")
  for (i in seq_len(nrow(matched))) {
    cat("   ", matched$Pathway_ID[i], " | score=", matched$Score[i], " | ", matched$Pathway_Name[i], "\n", sep = "")
  }

  matched
}

extract_gene_symbol_from_desc <- function(desc) {
  desc <- safe_trim(desc)
  if (is_blank(desc)) return(NA_character_)

  # Example:
  # "AKT1, AKT serine/threonine kinase 1"
  # "TP53, tumor protein p53"
  parts <- strsplit(desc, ",")[[1]]
  if (length(parts) == 0) return(NA_character_)

  symbol <- safe_trim(parts[1])
  if (is_blank(symbol)) return(NA_character_)

  symbol
}

parse_pathway_genes <- function(pathway_entry, phenotype, pathway_id, pathway_name) {
  gene_field <- pathway_entry$GENE

  if (is.null(gene_field) || length(gene_field) == 0) {
    return(empty_result())
  }

  gene_field <- as.character(gene_field)

  # KEGG pathway GENE field is alternating:
  # [1] "10458" [2] "BAIAP2, BAR/IMD domain containing adaptor protein 2"
  # [3] "10743" [4] "AKT1, AKT serine/threonine kinase 1"
  #
  # So we parse pairs: odd=index gene id, even=index description
  n <- length(gene_field)
  if (n < 2) {
    return(empty_result())
  }

  out <- list()
  idx <- seq(1, n - 1, by = 2)

  for (i in idx) {
    gene_id_raw <- safe_trim(gene_field[i])
    gene_desc   <- safe_trim(gene_field[i + 1])

    if (is_blank(gene_id_raw) || is_blank(gene_desc)) next

    gene_symbol <- extract_gene_symbol_from_desc(gene_desc)
    if (is.na(gene_symbol) || !nzchar(gene_symbol)) next

    kegg_gene_id <- paste0("hsa:", gene_id_raw)

    out[[length(out) + 1]] <- data.frame(
      Phenotype = phenotype,
      Pathway_ID = pathway_id,
      Pathway_Name = pathway_name,
      KEGG_Gene_ID = kegg_gene_id,
      Gene = gene_symbol,
      Source = "KEGG",
      stringsAsFactors = FALSE
    )
  }

  if (length(out) == 0) {
    return(empty_result())
  }

  do.call(rbind, out)
}

download_kegg_genes_from_pathways <- function(phenotype) {
  matched_pathways <- find_matching_pathways(phenotype, min_score = 1)

  if (nrow(matched_pathways) == 0) {
    return(empty_result())
  }

  all_rows <- list()

  cat("🧬 Retrieving genes directly from matched KEGG pathway entries...\n")

  for (i in seq_len(nrow(matched_pathways))) {
    pid <- matched_pathways$Pathway_ID[i]
    pname <- matched_pathways$Pathway_Name[i]

    pathway_info <- tryCatch(
      KEGGREST::keggGet(pid),
      error = function(e) {
        cat("   ⚠️ Failed to fetch pathway", pid, ":", conditionMessage(e), "\n")
        NULL
      }
    )

    if (is.null(pathway_info) || length(pathway_info) == 0 || is.null(pathway_info[[1]])) {
      next
    }

    one_df <- parse_pathway_genes(
      pathway_entry = pathway_info[[1]],
      phenotype = phenotype,
      pathway_id = pid,
      pathway_name = pname
    )

    if (nrow(one_df) > 0) {
      cat("   ", pid, " — ", nrow(one_df), " genes — ", pname, "\n", sep = "")
      all_rows[[length(all_rows) + 1]] <- one_df
    } else {
      cat("   ", pid, " — 0 genes parsed — ", pname, "\n", sep = "")
    }

    Sys.sleep(0.15)
  }

  if (length(all_rows) == 0) {
    return(empty_result())
  }

  result <- do.call(rbind, all_rows)

  if (is.null(result) || !is.data.frame(result) || nrow(result) == 0) {
    return(empty_result())
  }

  result <- result[!duplicated(result[, c("Pathway_ID", "KEGG_Gene_ID")]), , drop = FALSE]
  result <- result[order(result$Pathway_ID, result$Gene), , drop = FALSE]

  cat("   ✅ Total pathway-gene rows:", nrow(result), "\n")
  cat("   ✅ Total unique genes:", length(unique(result$Gene)), "\n")

  result
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("Usage: Rscript kegg.R <phenotype>\n")
    cat("Examples:\n")
    cat("  Rscript kegg.R migraine\n")
    cat("  Rscript kegg.R cancer\n")
    cat("  Rscript kegg.R \"breast cancer\"\n")
    quit(status = 1)
  }

  phenotype <- paste(args, collapse = " ")

  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  clean_phenotype <- clean_filename(phenotype)
  output_file <- file.path(output_dir, paste0(clean_phenotype, "_kegg.csv"))
  genes_only_file <- file.path(output_dir, paste0(clean_phenotype, "_kegg_genes.csv"))
  pathways_only_file <- file.path(output_dir, paste0(clean_phenotype, "_kegg_pathways.csv"))

  cat("🎯 Phenotype:           ", phenotype, "\n")
  cat("📁 Output directory:    ", output_dir, "\n")
  cat("⏰ Start time:          ", format(Sys.time()), "\n\n")

  results <- tryCatch(
    download_kegg_genes_from_pathways(phenotype),
    error = function(e) {
      cat("❌ Unexpected error:", conditionMessage(e), "\n")
      empty_result()
    }
  )

  if (is.null(results) || !is.data.frame(results)) {
    results <- empty_result()
  }

  if (nrow(results) > 0) {
    write.csv(results, output_file, row.names = FALSE)

    genes_only <- data.frame(
      Gene = sort(unique(results$Gene)),
      stringsAsFactors = FALSE
    )
    write.csv(genes_only, genes_only_file, row.names = FALSE)

    pathways_only <- unique(results[, c("Phenotype", "Pathway_ID", "Pathway_Name"), drop = FALSE])
    pathways_only <- pathways_only[order(pathways_only$Pathway_ID), , drop = FALSE]
    write.csv(pathways_only, pathways_only_file, row.names = FALSE)

    cat("\n✅ SUCCESS!\n")
    cat("💾 Full results:        ", output_file, "\n")
    cat("🧬 Genes-only file:     ", genes_only_file, "\n")
    cat("🛤️  Pathways-only file: ", pathways_only_file, "\n")
    cat("📊 Pathway-gene rows:   ", nrow(results), "\n")
    cat("🧬 Unique genes:        ", length(unique(results$Gene)), "\n")
    cat("🛤️  Matched pathways:   ", nrow(pathways_only), "\n\n")

    cat("🧬 Example genes:\n")
    ug <- sort(unique(results$Gene))
    for (i in seq(1, min(length(ug), 40), by = 8)) {
      j <- min(i + 7, length(ug), 40)
      cat("   ", paste(ug[i:j], collapse = ", "), "\n")
    }
  } else {
    cat("❌ No genes found for:", phenotype, "\n")
    cat("This means either:\n")
    cat("   1. no KEGG human pathway title matched the phenotype text, or\n")
    cat("   2. matched pathways did not return a parseable GENE field.\n")
  }

  cat("\n⏰ End time: ", format(Sys.time()), "\n")
}

if (!interactive()) {
  main()
}