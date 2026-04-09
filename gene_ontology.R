#!/usr/bin/env Rscript

# ============================================================
# Gene Ontology phenotype-to-gene downloader
# ------------------------------------------------------------
# Usage:
#   Rscript gene_ontology.R migraine
#   Rscript gene_ontology.R "breast cancer"
#
# Output:
#   AllPackagesGenes/<phenotype>_gene_ontology_full.csv
#   AllPackagesGenes/<phenotype>_gene_ontology_genes.csv
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

install_and_load <- function(pkg, bioc = FALSE) {
  if (!suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE))) {
    cat("Installing", pkg, "...\n")
    if (bioc) {
      if (!suppressWarnings(require("BiocManager", quietly = TRUE))) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
}

# Required packages
install_and_load("AnnotationDbi", bioc = TRUE)
install_and_load("GO.db", bioc = TRUE)
install_and_load("org.Hs.eg.db", bioc = TRUE)

safe_cat <- function(...) {
  cat(..., "\n", sep = "")
}

clean_filename <- function(x) {
  x <- gsub("[^[:alnum:]_ -]", "", x)   # hyphen at end
  x <- gsub("\\s+", "_", trimws(x))
  if (nchar(x) == 0) x <- "phenotype"
  x
}

build_search_patterns <- function(phenotype) {
  phenotype <- tolower(trimws(phenotype))
  phenotype_escaped <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", phenotype)

  words <- unlist(strsplit(phenotype, "[[:space:]/_-]+"))
  words <- words[nchar(words) > 1]

  patterns <- c(
    paste0("\\b", phenotype_escaped, "\\b"),
    paste0("\\b", gsub("[[:space:]]+", ".*", phenotype_escaped), "\\b")
  )

  if (length(words) > 1) {
    patterns <- c(
      patterns,
      paste(words, collapse = ".*"),
      paste(rev(words), collapse = ".*")
    )
  }

  patterns <- c(
    patterns,
    paste0("\\b", phenotype_escaped, " disease\\b"),
    paste0("\\b", phenotype_escaped, " disorder\\b"),
    paste0("\\b", phenotype_escaped, " syndrome\\b"),
    paste0("\\b", phenotype_escaped, " process\\b"),
    paste0("\\b", phenotype_escaped, " pathway\\b"),
    paste0("\\b", phenotype_escaped, " signaling\\b"),
    paste0("\\b", phenotype_escaped, " response\\b")
  )

  unique(patterns)
}

get_go_metadata <- function() {
  go_ids <- AnnotationDbi::keys(GO.db, keytype = "GOID")

  term_vec <- AnnotationDbi::Term(GOTERM[go_ids])
  ont_vec  <- AnnotationDbi::Ontology(GOTERM[go_ids])

  def_vec <- vapply(
    go_ids,
    function(id) {
      obj <- GOTERM[[id]]
      if (is.null(obj)) {
        return(NA_character_)
      }
      d <- tryCatch(Definition(obj), error = function(e) NA_character_)
      if (length(d) == 0) NA_character_ else d
    },
    FUN.VALUE = character(1)
  )

  data.frame(
    GOID = go_ids,
    TERM = unname(term_vec),
    DEFINITION = unname(def_vec),
    ONTOLOGY = unname(ont_vec),
    stringsAsFactors = FALSE
  )
}

find_relevant_go_terms <- function(phenotype, go_meta, max_terms = 100) {
  patterns <- build_search_patterns(phenotype)

  term_text <- tolower(ifelse(is.na(go_meta$TERM), "", go_meta$TERM))
  def_text  <- tolower(ifelse(is.na(go_meta$DEFINITION), "", go_meta$DEFINITION))

  matched_rows <- integer(0)
  match_reason <- character(0)

  for (pat in patterns) {
    idx_term <- which(grepl(pat, term_text, perl = TRUE))
    idx_def  <- which(grepl(pat, def_text, perl = TRUE))

    if (length(idx_term) > 0) {
      matched_rows <- c(matched_rows, idx_term)
      match_reason <- c(match_reason, rep(paste0("TERM:", pat), length(idx_term)))
    }

    if (length(idx_def) > 0) {
      matched_rows <- c(matched_rows, idx_def)
      match_reason <- c(match_reason, rep(paste0("DEFINITION:", pat), length(idx_def)))
    }
  }

  if (length(matched_rows) == 0) {
    return(data.frame())
  }

  matched_df <- data.frame(
    row_id = matched_rows,
    Match_Reason = match_reason,
    stringsAsFactors = FALSE
  )

  matched_df <- merge(
    matched_df,
    cbind(row_id = seq_len(nrow(go_meta)), go_meta),
    by = "row_id",
    all.x = TRUE
  )

  matched_df <- matched_df[!duplicated(matched_df$GOID), ]
  matched_df <- matched_df[order(matched_df$ONTOLOGY, matched_df$TERM), ]

  if (nrow(matched_df) > max_terms) {
    matched_df <- matched_df[seq_len(max_terms), ]
  }

  matched_df
}

get_genes_for_go_terms <- function(go_ids, phenotype, go_meta_subset) {
  if (length(go_ids) == 0) {
    return(data.frame())
  }

  safe_cat("   Mapping GO terms to human genes using org.Hs.eg.db ...")

  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(go_ids),
    keytype = "GOALL",
    columns = c("SYMBOL", "GENENAME", "ENTREZID", "GOALL", "ONTOLOGYALL")
  )

  if (is.null(gene_map) || nrow(gene_map) == 0) {
    return(data.frame())
  }

  gene_map <- gene_map[!is.na(gene_map$SYMBOL) & gene_map$SYMBOL != "", , drop = FALSE]
  gene_map <- gene_map[!is.na(gene_map$GOALL) & gene_map$GOALL != "", , drop = FALSE]

  if (nrow(gene_map) == 0) {
    return(data.frame())
  }

  colnames(go_meta_subset)[colnames(go_meta_subset) == "GOID"] <- "GOALL"

  merged <- merge(
    gene_map,
    go_meta_subset[, c("GOALL", "TERM", "DEFINITION", "ONTOLOGY", "Match_Reason")],
    by = "GOALL",
    all.x = TRUE
  )

  if (nrow(merged) == 0) {
    return(data.frame())
  }

  out <- data.frame(
    Gene = merged$SYMBOL,
    Gene_Name = merged$GENENAME,
    Entrez_ID = merged$ENTREZID,
    GO_ID = merged$GOALL,
    GO_Term = merged$TERM,
    GO_Definition = merged$DEFINITION,
    GO_Ontology = merged$ONTOLOGY,
    Match_Reason = merged$Match_Reason,
    Source = "Gene_Ontology",
    Phenotype = phenotype,
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$Gene) & out$Gene != "", , drop = FALSE]
  out <- out[order(out$GO_Ontology, out$GO_Term, out$Gene), , drop = FALSE]

  out
}

download_go_genes <- function(phenotype, max_terms = 100) {
  safe_cat("🔍 Searching Gene Ontology for phenotype: ", phenotype)

  go_meta <- get_go_metadata()
  safe_cat("   Total GO terms available: ", nrow(go_meta))

  relevant_go <- find_relevant_go_terms(phenotype, go_meta, max_terms = max_terms)

  if (nrow(relevant_go) == 0) {
    safe_cat("   No relevant GO terms found for: ", phenotype)
    return(list(full = data.frame(), genes_only = data.frame(), go_terms = data.frame()))
  }

  safe_cat("   Relevant GO terms found: ", nrow(relevant_go))

  full_results <- get_genes_for_go_terms(
    go_ids = relevant_go$GOID,
    phenotype = phenotype,
    go_meta_subset = relevant_go
  )

  if (nrow(full_results) == 0) {
    safe_cat("   No genes found for matched GO terms")
    return(list(full = data.frame(), genes_only = data.frame(), go_terms = relevant_go))
  }

  genes_only <- data.frame(Gene = sort(unique(full_results$Gene)), stringsAsFactors = FALSE)

  safe_cat("   Gene-term records found: ", nrow(full_results))
  safe_cat("   Unique genes found: ", nrow(genes_only))

  list(
    full = full_results,
    genes_only = genes_only,
    go_terms = relevant_go
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    safe_cat("Gene Ontology phenotype-to-gene downloader")
    safe_cat("Usage:   Rscript gene_ontology.R <phenotype>")
    safe_cat("Example: Rscript gene_ontology.R migraine")
    quit(status = 1)
  }

  phenotype <- args[1]

  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  clean_phenotype <- clean_filename(phenotype)

  full_file      <- file.path(output_dir, paste0(clean_phenotype, "_gene_ontology_full.csv"))
  genes_only_file <- file.path(output_dir, paste0(clean_phenotype, "_gene_ontology_genes.csv"))
  terms_file     <- file.path(output_dir, paste0(clean_phenotype, "_gene_ontology_terms.csv"))

  safe_cat("🎯 Phenotype: ", phenotype)
  safe_cat("📁 Output directory: ", output_dir)
  safe_cat("⏰ Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  safe_cat("")

  result <- tryCatch(
    download_go_genes(phenotype, max_terms = 100),
    error = function(e) {
      safe_cat("❌ Error: ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(result)) {
    quit(status = 1)
  }

  if (nrow(result$go_terms) > 0) {
    write.csv(result$go_terms, terms_file, row.names = FALSE)
  }

  if (nrow(result$full) > 0) {
    write.csv(result$full, full_file, row.names = FALSE)
    write.csv(result$genes_only, genes_only_file, row.names = FALSE)

    safe_cat("")
    safe_cat("✅ SUCCESS")
    safe_cat("💾 Full results: ", full_file)
    safe_cat("🧬 Genes-only file: ", genes_only_file)
    safe_cat("📘 Matched GO terms file: ", terms_file)
    safe_cat("📊 Total gene-term records: ", nrow(result$full))
    safe_cat("🧬 Total unique genes: ", nrow(result$genes_only))

    if ("GO_Ontology" %in% colnames(result$full)) {
      ont_summary <- table(result$full$GO_Ontology)
      safe_cat("")
      safe_cat("Genes by GO ontology:")
      for (ont in names(ont_summary)) {
        ont_name <- switch(
          ont,
          BP = "Biological Process",
          MF = "Molecular Function",
          CC = "Cellular Component",
          ont
        )
        safe_cat("   ", ont_name, " (", ont, "): ", as.integer(ont_summary[[ont]]))
      }
    }

    safe_cat("")
    safe_cat("Top genes:")
    preview_genes <- head(result$genes_only$Gene, 30)
    for (i in seq(1, length(preview_genes), by = 10)) {
      j <- min(i + 9, length(preview_genes))
      safe_cat("   ", paste(preview_genes[i:j], collapse = ", "))
    }

  } else {
    safe_cat("❌ No genes found for phenotype: ", phenotype)
    safe_cat("   This usually means GO term text did not match the phenotype well.")
    safe_cat("   Try a related keyword, for example:")
    safe_cat("   - migraine -> neuroinflammatory, pain, vascular")
    safe_cat("   - diabetes -> insulin, glucose, metabolic")
    safe_cat("   - asthma -> immune, inflammatory, airway")
  }

  safe_cat("")
  safe_cat("⏰ End time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

if (!interactive()) {
  main()
}