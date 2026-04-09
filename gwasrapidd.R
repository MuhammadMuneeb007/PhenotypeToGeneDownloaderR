#!/usr/bin/env Rscript
# Smart GWAS Catalog gene downloader using gwasrapidd
# Works better for terms like "migraine" by trying phenotype expansions.

required_cran <- c("gwasrapidd", "dplyr")

for (pkg in required_cran) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

safe_trim <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x)] <- ""
  x
}

clean_filename <- function(x) {
  x <- gsub("[^[:alnum:]_ -]", "", x)   # hyphen at end
  x <- gsub("\\s+", "_", trimws(x))
  if (nchar(x) == 0) x <- "phenotype"
  x
}

empty_gene_df <- function() {
  data.frame(
    Gene = character(0),
    Variant_ID = character(0),
    Chromosome = character(0),
    Position = numeric(0),
    Distance = numeric(0),
    Is_Mapped_Gene = logical(0),
    Is_Closest_Gene = logical(0),
    Is_Intergenic = logical(0),
    Is_Upstream = logical(0),
    Is_Downstream = logical(0),
    Mapping_Source = character(0),
    Mapping_Method = character(0),
    Search_Type = character(0),
    Query_Term = character(0),
    Phenotype = character(0),
    stringsAsFactors = FALSE
  )
}

build_search_terms <- function(phenotype) {
  p <- tolower(safe_trim(phenotype))

  terms <- c(phenotype)

  # Generic expansions
  terms <- c(
    terms,
    paste(phenotype, "disorder"),
    paste(phenotype, "trait"),
    paste(phenotype, "disease")
  )

  # Phenotype-specific expansions
  if (grepl("migraine", p)) {
    terms <- c(
      terms,
      "migraine disorder",
      "migraine with aura",
      "migraine without aura",
      "headache",
      "headache disorder",
      "familial hemiplegic migraine"
    )
  }

  if (grepl("cancer|tumou?r|neoplasm", p)) {
    terms <- c(
      terms,
      "cancer",
      "neoplasm",
      "tumor",
      "tumour",
      "malignancy"
    )
  }

  if (grepl("diabetes", p)) {
    terms <- c(
      terms,
      "diabetes mellitus",
      "type 1 diabetes",
      "type 2 diabetes"
    )
  }

  terms <- unique(safe_trim(terms))
  terms[terms != ""]
}

discover_trait_terms <- function(phenotype) {
  found_terms <- character(0)

  trait_obj <- tryCatch(
    gwasrapidd::get_traits(
      efo_trait = phenotype,
      set_operation = "union",
      interactive = FALSE,
      warnings = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(trait_obj) &&
      methods::is(trait_obj, "traits") &&
      length(trait_obj@traits) > 0 &&
      nrow(trait_obj@traits) > 0) {

    trait_df <- trait_obj@traits

    candidate_cols <- intersect(
      c("efo_trait", "trait", "trait_name", "uri"),
      colnames(trait_df)
    )

    for (cc in candidate_cols) {
      found_terms <- c(found_terms, safe_trim(trait_df[[cc]]))
    }
  }

  found_terms <- unique(found_terms)
  found_terms[found_terms != ""]
}

extract_genes_from_variants <- function(variants_obj, phenotype, search_type, query_term) {
  if (is.null(variants_obj) || !methods::is(variants_obj, "variants")) {
    return(empty_gene_df())
  }

  if (length(variants_obj@genomic_contexts) == 0 || nrow(variants_obj@genomic_contexts) == 0) {
    return(empty_gene_df())
  }

  gc <- variants_obj@genomic_contexts

  if (!all(c("variant_id", "gene_name") %in% colnames(gc))) {
    return(empty_gene_df())
  }

  out <- data.frame(
    Gene = safe_trim(gc$gene_name),
    Variant_ID = safe_trim(gc$variant_id),
    Chromosome = if ("chromosome_name" %in% colnames(gc)) safe_trim(gc$chromosome_name) else "",
    Position = if ("chromosome_position" %in% colnames(gc)) gc$chromosome_position else NA_real_,
    Distance = if ("distance" %in% colnames(gc)) gc$distance else NA_real_,
    Is_Mapped_Gene = if ("is_mapped_gene" %in% colnames(gc)) gc$is_mapped_gene else NA,
    Is_Closest_Gene = if ("is_closest_gene" %in% colnames(gc)) gc$is_closest_gene else NA,
    Is_Intergenic = if ("is_intergenic" %in% colnames(gc)) gc$is_intergenic else NA,
    Is_Upstream = if ("is_upstream" %in% colnames(gc)) gc$is_upstream else NA,
    Is_Downstream = if ("is_downstream" %in% colnames(gc)) gc$is_downstream else NA,
    Mapping_Source = if ("source" %in% colnames(gc)) safe_trim(gc$source) else "",
    Mapping_Method = if ("mapping_method" %in% colnames(gc)) safe_trim(gc$mapping_method) else "",
    Search_Type = search_type,
    Query_Term = query_term,
    Phenotype = phenotype,
    stringsAsFactors = FALSE
  )

  out <- out %>%
    dplyr::filter(!is.na(Gene), Gene != "", !is.na(Variant_ID), Variant_ID != "") %>%
    dplyr::distinct()

  out
}

run_variant_search <- function(query_term, phenotype, search_type) {
  result <- tryCatch(
    {
      if (search_type == "EFO_trait") {
        gwasrapidd::get_variants(
          efo_trait = query_term,
          set_operation = "union",
          interactive = FALSE,
          warnings = FALSE
        )
      } else {
        gwasrapidd::get_variants(
          reported_trait = query_term,
          set_operation = "union",
          interactive = FALSE,
          warnings = FALSE
        )
      }
    },
    error = function(e) {
      cat("   ⚠️", search_type, "search failed for", shQuote(query_term), ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  extract_genes_from_variants(result, phenotype, search_type, query_term)
}

download_gwas_genes <- function(phenotype) {
  cat("🔍 Searching GWAS Catalog for:", phenotype, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  search_terms <- build_search_terms(phenotype)

  # Also try trait discovery from GWAS Catalog itself
  discovered_terms <- discover_trait_terms(phenotype)
  if (length(discovered_terms) > 0) {
    cat("   Trait discovery found", length(discovered_terms), "additional candidate terms\n")
    search_terms <- unique(c(search_terms, discovered_terms))
  }

  cat("   Candidate search terms:\n")
  for (term in search_terms) {
    cat("   -", term, "\n")
  }

  all_rows <- list()
  idx <- 1

  for (term in search_terms) {
    cat("\n   Query term:", term, "\n")

    efo_df <- run_variant_search(term, phenotype, "EFO_trait")
    if (nrow(efo_df) > 0) {
      cat("     EFO trait hits:", nrow(efo_df), "\n")
      all_rows[[idx]] <- efo_df
      idx <- idx + 1
    } else {
      cat("     EFO trait hits: 0\n")
    }

    Sys.sleep(0.3)

    reported_df <- run_variant_search(term, phenotype, "reported_trait")
    if (nrow(reported_df) > 0) {
      cat("     Reported trait hits:", nrow(reported_df), "\n")
      all_rows[[idx]] <- reported_df
      idx <- idx + 1
    } else {
      cat("     Reported trait hits: 0\n")
    }

    Sys.sleep(0.3)
  }

  if (length(all_rows) == 0) {
    cat("\n❌ No GWAS variant gene mappings found for:", phenotype, "\n")
    return(empty_gene_df())
  }

  combined <- dplyr::bind_rows(all_rows)

  final_results <- combined %>%
    dplyr::filter(!is.na(Gene), Gene != "") %>%
    dplyr::arrange(
      dplyr::desc(Is_Mapped_Gene),
      dplyr::desc(Is_Closest_Gene),
      Gene,
      Variant_ID
    ) %>%
    dplyr::distinct(Gene, Variant_ID, Search_Type, Query_Term, .keep_all = TRUE)

  cat("\n🧹 Cleaning results...\n")
  cat("   Raw gene-context rows:", nrow(combined), "\n")
  cat("   Final gene-context rows:", nrow(final_results), "\n")
  cat("   Unique genes:", length(unique(final_results$Gene)), "\n")

  final_results
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("Smart GWAS Catalog Gene Downloader\n")
    cat("Usage:   Rscript gwasrapidd_smart.R <phenotype>\n")
    cat("Example: Rscript gwasrapidd_smart.R migraine\n")
    quit(status = 1)
  }

  phenotype <- args[1]

  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  clean_phenotype <- clean_filename(phenotype)
  output_file <- file.path(output_dir, paste0(clean_phenotype, "_gwasrapidd.csv"))
  genes_only_file <- file.path(output_dir, paste0(clean_phenotype, "_gwasrapidd_genes.csv"))

  cat("🎯 Phenotype:        ", phenotype, "\n")
  cat("📁 Output directory: ", output_dir, "\n")
  cat("📄 Full output:      ", output_file, "\n")
  cat("🧬 Genes-only file:  ", genes_only_file, "\n")
  cat("⏰ Start time:       ", format(Sys.time()), "\n\n")

  results <- download_gwas_genes(phenotype)

  if (nrow(results) > 0) {
    write.csv(results, output_file, row.names = FALSE)

    genes_only <- data.frame(
      Gene = sort(unique(results$Gene)),
      stringsAsFactors = FALSE
    )
    write.csv(genes_only, genes_only_file, row.names = FALSE)

    cat("\n✅ SUCCESS!\n")
    cat("📊 Total records:   ", nrow(results), "\n")
    cat("🧬 Unique genes:    ", nrow(genes_only), "\n")
    cat("💾 Full results:    ", output_file, "\n")
    cat("🧬 Genes-only file: ", genes_only_file, "\n\n")

    cat("🧬 Genes found for", phenotype, ":\n")
    unique_genes <- genes_only$Gene
    for (i in seq(1, length(unique_genes), by = 8)) {
      end_idx <- min(i + 7, length(unique_genes))
      cat("   ", paste(unique_genes[i:end_idx], collapse = ", "), "\n")
    }

    cat("\n📋 Example rows:\n")
    print(utils::head(results, 10))
  } else {
    cat("❌ No genes found for:", phenotype, "\n")
  }

  cat("\n⏰ End time: ", format(Sys.time()), "\n")
}

if (!interactive()) {
  main()
}