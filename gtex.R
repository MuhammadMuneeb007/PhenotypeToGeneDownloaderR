#!/usr/bin/env Rscript
# gtex.R
# Dynamic phenotype -> tissue -> gene GTEx script
# Uses PubMed counts to rank tissues dynamically, then retrieves GTEx eGenes.

if (!interactive() && is.null(getOption("repos"))) {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

options(gtexr.itemsPerPage = 100000)
options(gtexr.verbose = FALSE)

required_packages <- c("gtexr", "dplyr", "httr", "jsonlite")

install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

for (pkg in required_packages) {
  install_and_load(pkg)
}

msg <- function(...) cat(..., "\n", sep = "")

clean_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

fetch_available_tissues <- function() {
  msg("📥 Fetching GTEx tissues...")

  out <- tryCatch({
    x <- gtexr::get_tissue_site_detail()

    id_col <- if ("tissueSiteDetailId" %in% names(x)) "tissueSiteDetailId" else names(x)[1]
    label_col <- if ("tissueSiteDetail" %in% names(x)) "tissueSiteDetail" else id_col

    df <- data.frame(
      Tissue_ID = as.character(x[[id_col]]),
      Tissue_Label = as.character(x[[label_col]]),
      stringsAsFactors = FALSE
    )

    df <- df[!is.na(df$Tissue_ID) & nzchar(df$Tissue_ID), , drop = FALSE]
    df$Tissue_Label[is.na(df$Tissue_Label) | !nzchar(df$Tissue_Label)] <- gsub("_", " ", df$Tissue_ID[is.na(df$Tissue_Label) | !nzchar(df$Tissue_Label)])
    df <- df[!duplicated(df$Tissue_ID), , drop = FALSE]
    rownames(df) <- NULL
    df
  }, error = function(e) {
    msg("❌ Could not fetch GTEx tissues: ", conditionMessage(e))
    data.frame()
  })

  out
}

normalize_tissue_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\(.*?\\)", "", x)
  x <- gsub("-", " ", x)
  x <- gsub("\\bBA[0-9]+\\b", "", x, ignore.case = TRUE)
  x <- gsub("\\bbasal ganglia\\b", "", x, ignore.case = TRUE)
  x <- gsub("\\bcervical c[- ]?1\\b", "spinal cord", x, ignore.case = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

build_pubmed_query <- function(phenotype, tissue_label) {
  ph <- gsub('"', "", phenotype)
  ti <- gsub('"', "", tissue_label)

  paste0('("', ph, '"[Title/Abstract]) AND ("', ti, '"[Title/Abstract])')
}

get_pubmed_count <- function(query) {
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

  resp <- tryCatch(
    httr::GET(
      base_url,
      query = list(
        db = "pubmed",
        term = query,
        retmode = "json"
      ),
      httr::timeout(45)
    ),
    error = function(e) NULL
  )

  if (is.null(resp)) return(NA_real_)
  if (httr::status_code(resp) != 200) return(NA_real_)

  txt <- tryCatch(httr::content(resp, "text", encoding = "UTF-8"), error = function(e) "")
  if (!nzchar(txt)) return(NA_real_)

  obj <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(obj)) return(NA_real_)

  count <- obj$esearchresult$count
  safe_num(count)
}

rank_tissues_by_pubmed <- function(phenotype, tissues_df, max_tissues = 8) {
  if (nrow(tissues_df) == 0) return(data.frame())

  msg("📚 Ranking tissues dynamically with PubMed...")

  res <- vector("list", nrow(tissues_df))

  for (i in seq_len(nrow(tissues_df))) {
    tissue_id <- tissues_df$Tissue_ID[i]
    tissue_label <- normalize_tissue_label(tissues_df$Tissue_Label[i])
    query <- build_pubmed_query(phenotype, tissue_label)
    count <- get_pubmed_count(query)

    res[[i]] <- data.frame(
      Tissue_ID = tissue_id,
      Tissue_Label = tissues_df$Tissue_Label[i],
      Tissue_Query = tissue_label,
      Literature_Count = count,
      Query = query,
      stringsAsFactors = FALSE
    )

    msg(sprintf("   [%d/%d] %-40s -> %s",
                i, nrow(tissues_df), tissue_label,
                ifelse(is.na(count), "NA", as.character(count))))
    Sys.sleep(0.2)
  }

  ranked <- do.call(rbind, res)
  ranked <- ranked[!is.na(ranked$Literature_Count), , drop = FALSE]

  if (nrow(ranked) == 0) {
    msg("⚠️ No literature counts were retrieved.")
    return(data.frame())
  }

  ranked <- ranked[order(-ranked$Literature_Count, ranked$Tissue_Label), , drop = FALSE]

  # Prefer positive counts
  pos <- ranked[ranked$Literature_Count > 0, , drop = FALSE]
  if (nrow(pos) == 0) {
    msg("⚠️ No positive PubMed matches found; using top tissues with numeric counts.")
    pos <- ranked
  }

  pos <- head(pos, max_tissues)
  rownames(pos) <- NULL
  pos
}

fetch_tissue_egenes <- function(tissue, q_threshold = 0.05, max_genes_per_tissue = 200) {
  msg("🔬 Processing GTEx tissue: ", tissue)

  out <- tryCatch({
    x <- gtexr::get_eqtl_genes(
      tissueSiteDetailIds = tissue,
      datasetId = "gtex_v8",
      itemsPerPage = 100000,
      .verbose = FALSE
    )

    if (is.null(x) || nrow(x) == 0) {
      msg("   No eGene data")
      return(data.frame())
    }

    needed <- c(
      "tissueSiteDetailId",
      "geneSymbol",
      "gencodeId",
      "empiricalPValue",
      "pValue",
      "pValueThreshold",
      "qValue",
      "log2AllelicFoldChange"
    )

    for (nm in needed) {
      if (!nm %in% names(x)) x[[nm]] <- NA
    }

    x <- x[!is.na(x$geneSymbol) & nzchar(x$geneSymbol), , drop = FALSE]

    if (any(!is.na(x$qValue))) {
      x <- x[!is.na(x$qValue) & x$qValue <= q_threshold, , drop = FALSE]
    } else {
      x <- x[!is.na(x$pValue) & !is.na(x$pValueThreshold) & x$pValue <= x$pValueThreshold, , drop = FALSE]
    }

    if (nrow(x) == 0) {
      msg("   No genes passed threshold")
      return(data.frame())
    }

    ord_q <- ifelse(is.na(x$qValue), Inf, x$qValue)
    ord_ep <- ifelse(is.na(x$empiricalPValue), Inf, x$empiricalPValue)
    ord_p <- ifelse(is.na(x$pValue), Inf, x$pValue)

    x <- x[order(ord_q, ord_ep, ord_p), , drop = FALSE]
    x <- head(x, max_genes_per_tissue)

    data.frame(
      Tissue = x$tissueSiteDetailId,
      Gene = x$geneSymbol,
      Gencode_ID = x$gencodeId,
      QValue = x$qValue,
      Empirical_PValue = x$empiricalPValue,
      PValue = x$pValue,
      PValue_Threshold = x$pValueThreshold,
      Log2_AFC = x$log2AllelicFoldChange,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    msg("   ❌ Error: ", conditionMessage(e))
    data.frame()
  })

  Sys.sleep(0.25)
  out
}

summarize_genes <- function(full_df, phenotype, selected_tissues_df) {
  if (nrow(full_df) == 0) return(data.frame())

  tissue_best <- full_df |>
    dplyr::group_by(Gene, Tissue) |>
    dplyr::arrange(QValue, Empirical_PValue, PValue, .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  summary_df <- tissue_best |>
    dplyr::group_by(Gene) |>
    dplyr::summarise(
      Tissue_Count = dplyr::n_distinct(Tissue),
      Best_QValue = if (all(is.na(QValue))) NA_real_ else min(QValue, na.rm = TRUE),
      Best_Empirical_PValue = if (all(is.na(Empirical_PValue))) NA_real_ else min(Empirical_PValue, na.rm = TRUE),
      Best_PValue = if (all(is.na(PValue))) NA_real_ else min(PValue, na.rm = TRUE),
      Mean_Abs_Log2_AFC = if (all(is.na(Log2_AFC))) NA_real_ else mean(abs(Log2_AFC), na.rm = TRUE),
      Example_Tissues = paste(head(sort(unique(Tissue)), 10), collapse = ";"),
      Gencode_ID = dplyr::first(Gencode_ID),
      Phenotype = phenotype,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(Tissue_Count), Best_QValue, Best_Empirical_PValue, Best_PValue, Gene)

  if (nrow(selected_tissues_df) > 0) {
    summary_df$Selected_Tissues <- paste(
      paste0(selected_tissues_df$Tissue_ID, " (count=", selected_tissues_df$Literature_Count, ")"),
      collapse = "; "
    )
  } else {
    summary_df$Selected_Tissues <- ""
  }

  as.data.frame(summary_df)
}

save_outputs <- function(full_df, summary_df, tissues_df, phenotype) {
  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  prefix <- clean_filename(phenotype)

  tissues_file <- file.path(output_dir, paste0(prefix, "_gtex_ranked_tissues.csv"))
  full_file <- file.path(output_dir, paste0(prefix, "_gtex_tissue_eqtls.csv"))
  summary_file <- file.path(output_dir, paste0(prefix, "_gtex_prioritized_genes.csv"))
  genes_file <- file.path(output_dir, paste0(prefix, "_gtex_genes.csv"))

  write.csv(tissues_df, tissues_file, row.names = FALSE)
  write.csv(full_df, full_file, row.names = FALSE)
  write.csv(summary_df, summary_file, row.names = FALSE)
  write.csv(data.frame(Gene = summary_df$Gene, stringsAsFactors = FALSE), genes_file, row.names = FALSE)

  list(
    tissues_file = tissues_file,
    full_file = full_file,
    summary_file = summary_file,
    genes_file = genes_file
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    msg("Usage: Rscript gtex.R <phenotype> [q_value_threshold] [max_genes_per_tissue] [max_tissues]")
    msg('Example: Rscript gtex.R "migraine"')
    msg('Example: Rscript gtex.R "type 2 diabetes" 0.05 200 8')
    quit(status = 1)
  }

  phenotype <- args[1]
  q_threshold <- if (length(args) >= 2) as.numeric(args[2]) else 0.05
  max_genes_per_tissue <- if (length(args) >= 3) as.integer(args[3]) else 200
  max_tissues <- if (length(args) >= 4) as.integer(args[4]) else 8

  if (is.na(q_threshold) || q_threshold <= 0 || q_threshold > 1) {
    stop("q_value_threshold must be between 0 and 1")
  }
  if (is.na(max_genes_per_tissue) || max_genes_per_tissue <= 0) {
    stop("max_genes_per_tissue must be > 0")
  }
  if (is.na(max_tissues) || max_tissues <= 0) {
    stop("max_tissues must be > 0")
  }

  msg("🧬 DYNAMIC PHENOTYPE → TISSUE → GENE GTEx SCRIPT")
  msg("══════════════════════════════════════════════════")
  msg("🎯 Phenotype: ", phenotype)
  msg("📉 q-value threshold: ", q_threshold)
  msg("📊 max genes per tissue: ", max_genes_per_tissue)
  msg("🫀 max tissues: ", max_tissues)
  msg("⏰ Start time: ", format(Sys.time()))
  msg("")

  tissues_df <- fetch_available_tissues()
  if (nrow(tissues_df) == 0) {
    stop("No GTEx tissues could be retrieved.")
  }

  ranked_tissues <- rank_tissues_by_pubmed(
    phenotype = phenotype,
    tissues_df = tissues_df,
    max_tissues = max_tissues
  )

  if (nrow(ranked_tissues) == 0) {
    msg("")
    msg("❌ No dynamic tissue matches found for phenotype: ", phenotype)
    msg("This usually means PubMed could not be reached from the HPC or all queries failed.")
    quit(status = 0)
  }

  msg("")
  msg("✅ Selected tissues:")
  for (i in seq_len(nrow(ranked_tissues))) {
    msg(sprintf("  %2d. %-35s count=%s",
                i,
                ranked_tissues$Tissue_Label[i],
                ranked_tissues$Literature_Count[i]))
  }
  msg("")

  all_results <- list()

  for (i in seq_len(nrow(ranked_tissues))) {
    msg("[", i, "/", nrow(ranked_tissues), "]")
    res <- fetch_tissue_egenes(
      tissue = ranked_tissues$Tissue_ID[i],
      q_threshold = q_threshold,
      max_genes_per_tissue = max_genes_per_tissue
    )
    if (nrow(res) > 0) {
      all_results[[length(all_results) + 1]] <- res
    }
  }

  if (length(all_results) == 0) {
    msg("")
    msg("❌ No significant GTEx eGenes found for the dynamically selected tissues.")
    msg("Try a less strict threshold, e.g. 0.1")
    quit(status = 0)
  }

  full_df <- do.call(rbind, all_results)
  rownames(full_df) <- NULL

  summary_df <- summarize_genes(full_df, phenotype, ranked_tissues)
  saved <- save_outputs(full_df, summary_df, ranked_tissues, phenotype)

  msg("")
  msg("✅ SUCCESS!")
  msg("📄 Ranked tissues file: ", saved$tissues_file)
  msg("📊 Tissue-level rows: ", nrow(full_df))
  msg("🧬 Prioritized unique genes: ", nrow(summary_df))
  msg("💾 Full results: ", saved$full_file)
  msg("💾 Summary results: ", saved$summary_file)
  msg("🧬 Genes-only file: ", saved$genes_file)
  msg("")

  top_n <- min(20, nrow(summary_df))
  msg("🏆 Top ", top_n, " prioritized genes for ", phenotype, ":")

  top_df <- summary_df[seq_len(top_n), , drop = FALSE]
  for (i in seq_len(nrow(top_df))) {
    msg(sprintf(
      "  %2d. %-15s tissues=%2d best_q=%s best_p=%s",
      i,
      top_df$Gene[i],
      top_df$Tissue_Count[i],
      format(top_df$Best_QValue[i], scientific = TRUE, digits = 3),
      format(top_df$Best_PValue[i], scientific = TRUE, digits = 3)
    ))
  }

  msg("")
  msg("⏰ End time: ", format(Sys.time()))
}

if (!interactive()) {
  main()
}