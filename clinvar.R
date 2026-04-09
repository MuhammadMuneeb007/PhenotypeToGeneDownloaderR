#!/usr/bin/env Rscript

# ClinVar gene downloader (structured XML version)
# This script searches ClinVar for a phenotype, retrieves matching variation
# records in XML, and extracts gene symbols from structured XML fields
# instead of scraping free text.

required_packages <- c("rentrez", "xml2", "dplyr")

for (pkg in required_packages) {
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
  x <- gsub("[^[:alnum:]_ -]", "", x)
  x <- gsub("\\s+", "_", trimws(x))
  if (nchar(x) == 0) x <- "phenotype"
  x
}

unique_nonempty <- function(x) {
  x <- safe_trim(x)
  x <- unique(x[x != ""])
  x
}

extract_record_id <- function(record_node) {
  # Try several possible places for the variation ID
  id_candidates <- c(
    xml2::xml_attr(record_node, "VariationID"),
    xml2::xml_attr(record_node, "VariationId"),
    xml2::xml_attr(record_node, "AccessionVersion"),
    xml2::xml_text(xml2::xml_find_first(record_node, ".//*[local-name()='VariationID']")),
    xml2::xml_text(xml2::xml_find_first(record_node, ".//*[local-name()='VariationId']"))
  )
  id_candidates <- unique_nonempty(id_candidates)
  if (length(id_candidates) > 0) return(id_candidates[1])
  return(NA_character_)
}

extract_accession <- function(record_node) {
  acc_candidates <- c(
    xml2::xml_attr(record_node, "Accession"),
    xml2::xml_attr(record_node, "VariationName"),
    xml2::xml_text(xml2::xml_find_first(record_node, ".//*[local-name()='Accession']")),
    xml2::xml_text(xml2::xml_find_first(record_node, ".//*[local-name()='VariationName']"))
  )
  acc_candidates <- unique_nonempty(acc_candidates)
  if (length(acc_candidates) > 0) return(acc_candidates[1])
  return(NA_character_)
}

extract_clinical_significance <- function(record_node) {
  sig_candidates <- c(
    xml2::xml_attr(xml2::xml_find_first(record_node, ".//*[local-name()='ClassifiedRecord']/*[local-name()='ClinicalSignificance']"), "Description"),
    xml2::xml_attr(xml2::xml_find_first(record_node, ".//*[local-name()='ClinicalSignificance']"), "Description"),
    xml2::xml_text(xml2::xml_find_first(record_node, ".//*[local-name()='Description']"))
  )
  sig_candidates <- unique_nonempty(sig_candidates)
  if (length(sig_candidates) > 0) return(sig_candidates[1])
  return(NA_character_)
}

extract_conditions <- function(record_node) {
  condition_nodes <- xml2::xml_find_all(
    record_node,
    ".//*[local-name()='Trait']/*[local-name()='Name']/*[local-name()='ElementValue'] |
     .//*[local-name()='ConditionList']//*[local-name()='Name'] |
     .//*[local-name()='TraitSet']//*[local-name()='ElementValue']"
  )
  conditions <- unique_nonempty(xml2::xml_text(condition_nodes))
  if (length(conditions) == 0) return(NA_character_)
  paste(conditions, collapse = " | ")
}

extract_gene_symbols <- function(record_node) {
  # Structured extraction only. No generic uppercase regex scraping.
  gene_nodes <- xml2::xml_find_all(
    record_node,
    ".//*[local-name()='Gene'] |
     .//*[local-name()='GeneSymbol'] |
     .//*[local-name()='MeasureRelationship']/*[local-name()='Symbol']/*[local-name()='ElementValue'] |
     .//*[local-name()='Name']/*[local-name()='ElementValue'][@Type='PreferredSymbol']"
  )

  values <- c(
    xml2::xml_attr(gene_nodes, "Symbol"),
    xml2::xml_attr(gene_nodes, "symbol"),
    xml2::xml_attr(gene_nodes, "GeneSymbol"),
    xml2::xml_text(gene_nodes)
  )

  values <- unique_nonempty(values)

  # Conservative cleanup: keep gene-like symbols only
  values <- values[
    grepl("^[A-Z][A-Z0-9\\-]{1,14}$", values)
  ]

  # Exclude obvious non-gene artifacts
  exclude <- c(
    "DNA", "RNA", "PCR", "SNP", "CNV", "OMIM", "HGVS", "GRCH37", "GRCH38",
    "PATHOGENIC", "BENIGN", "UNCERTAIN", "SIGNIFICANCE", "SOMATIC"
  )
  values <- values[!values %in% exclude]

  unique(values)
}

parse_clinvar_xml_records <- function(xml_text, phenotype, search_term) {
  if (is.null(xml_text) || length(xml_text) == 0 || safe_trim(xml_text) == "") {
    return(data.frame())
  }

  doc <- tryCatch(
    xml2::read_xml(xml_text),
    error = function(e) NULL
  )

  if (is.null(doc)) {
    return(data.frame())
  }

  # Handle several possible record wrappers
  record_nodes <- xml2::xml_find_all(
    doc,
    ".//*[local-name()='VariationArchive'] |
     .//*[local-name()='ClinVarResult-Set']/*[local-name()='VariationArchive'] |
     .//*[local-name()='DocumentSummary']"
  )

  if (length(record_nodes) == 0) {
    return(data.frame())
  }

  out <- list()
  row_idx <- 1

  for (record_node in record_nodes) {
    variation_id <- extract_record_id(record_node)
    accession <- extract_accession(record_node)
    clinical_significance <- extract_clinical_significance(record_node)
    conditions <- extract_conditions(record_node)
    genes <- extract_gene_symbols(record_node)

    if (length(genes) == 0) next

    for (gene in genes) {
      out[[row_idx]] <- data.frame(
        Gene = gene,
        ClinVar_VariationID = ifelse(is.na(variation_id), "", variation_id),
        ClinVar_Accession = ifelse(is.na(accession), "", accession),
        Query_Phenotype = phenotype,
        Matched_Condition = ifelse(is.na(conditions), "", conditions),
        Clinical_Significance = ifelse(is.na(clinical_significance), "", clinical_significance),
        Search_Term = search_term,
        Source = "ClinVar_XML",
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1
    }
  }

  if (length(out) == 0) {
    return(data.frame())
  }

  dplyr::bind_rows(out)
}

download_clinvar_genes <- function(phenotype, max_results = 2000, batch_size = 50) {
  cat("🔍 ClinVar search for:", phenotype, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  all_results <- list()
  result_idx <- 1

  phenotype <- safe_trim(phenotype)
  phenotype_lower <- tolower(phenotype)

  search_terms <- unique(c(
    paste0("\"", phenotype, "\""),
    paste0("\"", phenotype, "\"[disease]"),
    paste0("\"", phenotype, "\"[phenotype]"),
    paste0("\"", phenotype, "\" AND pathogenic"),
    paste0("\"", phenotype, "\" AND \"likely pathogenic\"")
  ))
 

  for (search_num in seq_along(search_terms)) {
    search_term <- search_terms[search_num]
    cat("\n", search_num, ". Searching:", search_term, "\n")

    search_result <- tryCatch(
      rentrez::entrez_search(
        db = "clinvar",
        term = search_term,
        retmax = max_results
      ),
      error = function(e) {
        cat("   ❌ Search error:", conditionMessage(e), "\n")
        return(NULL)
      }
    )

    if (is.null(search_result) || length(search_result$ids) == 0) {
      cat("   No entries found\n")
      next
    }

    ids <- unique(search_result$ids)
    cat("   Found", length(ids), "variation IDs\n")

    n_batches <- ceiling(length(ids) / batch_size)

    for (batch_num in seq_len(n_batches)) {
      start_idx <- (batch_num - 1) * batch_size + 1
      end_idx <- min(batch_num * batch_size, length(ids))
      batch_ids <- ids[start_idx:end_idx]

      cat("     Batch", batch_num, "of", n_batches, "- IDs:", length(batch_ids), "\n")

      xml_text <- tryCatch(
        rentrez::entrez_fetch(
          db = "clinvar",
          id = batch_ids,
          rettype = "vcv",
          retmode = "xml",
          parsed = FALSE,
          config = NULL,
          is_variationid = TRUE
        ),
        error = function(e) {
          cat("       ❌ Fetch error:", conditionMessage(e), "\n")
          return(NULL)
        }
      )

      if (is.null(xml_text)) {
        Sys.sleep(0.34)
        next
      }

      parsed_batch <- tryCatch(
        parse_clinvar_xml_records(xml_text, phenotype = phenotype, search_term = search_term),
        error = function(e) {
          cat("       ❌ Parse error:", conditionMessage(e), "\n")
          return(data.frame())
        }
      )

      if (nrow(parsed_batch) > 0) {
        all_results[[result_idx]] <- parsed_batch
        result_idx <- result_idx + 1
      }

      # Respectful pacing for NCBI
      Sys.sleep(0.34)
    }
  }

  if (length(all_results) == 0) {
    return(data.frame())
  }

  combined <- dplyr::bind_rows(all_results)

  final_results <- combined %>%
    dplyr::filter(!is.na(Gene), Gene != "") %>%
    dplyr::distinct(Gene, ClinVar_VariationID, .keep_all = TRUE) %>%
    dplyr::arrange(Gene, ClinVar_VariationID)

  cat("\n🧹 Cleaning results...\n")
  cat("  Raw records:   ", nrow(combined), "\n")
  cat("  Final records: ", nrow(final_results), "\n")
  cat("  Unique genes:  ", length(unique(final_results$Gene)), "\n")

  final_results
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("ClinVar Gene Downloader (Structured XML)\n")
    cat("Usage:\n")
    cat("  Rscript clinvar_fixed.R <phenotype> [max_results]\n\n")
    cat("Examples:\n")
    cat("  Rscript clinvar_fixed.R migraine\n")
    cat("  Rscript clinvar_fixed.R diabetes 1000\n")
    quit(status = 1)
  }

  phenotype <- args[1]
  max_results <- if (length(args) > 1) as.numeric(args[2]) else 2000

  if (is.na(max_results) || max_results <= 0) {
    stop("max_results must be a positive number")
  }

  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  clean_phenotype <- clean_filename(phenotype)
  output_file <- file.path(output_dir, paste0(clean_phenotype, "_clinvar.csv"))
  genes_only_file <- file.path(output_dir, paste0(clean_phenotype, "_clinvar_genes.csv"))

  cat("🎯 Phenotype:        ", phenotype, "\n")
  cat("📁 Output directory: ", output_dir, "\n")
  cat("📄 Full output:      ", output_file, "\n")
  cat("🧬 Genes-only file:  ", genes_only_file, "\n")
  cat("⏰ Start time:       ", format(Sys.time()), "\n\n")

  results <- download_clinvar_genes(phenotype = phenotype, max_results = max_results, batch_size = 50)

  if (nrow(results) == 0) {
    cat("\n❌ No structured gene results found for:", phenotype, "\n")
    cat("Try a broader or alternative phenotype term.\n")
    cat("\n⏰ End time: ", format(Sys.time()), "\n")
    return(invisible(NULL))
  }

  write.csv(results, output_file, row.names = FALSE)

  genes_only <- data.frame(Gene = sort(unique(results$Gene)), stringsAsFactors = FALSE)
  write.csv(genes_only, genes_only_file, row.names = FALSE)

  cat("\n✅ SUCCESS!\n")
  cat("📊 Total records:   ", nrow(results), "\n")
  cat("🧬 Unique genes:    ", nrow(genes_only), "\n")
  cat("💾 Full results:    ", output_file, "\n")
  cat("🧬 Genes-only file: ", genes_only_file, "\n\n")

  cat("🧬 Genes found for", phenotype, ":\n")
  gene_vec <- genes_only$Gene
  for (i in seq(1, length(gene_vec), by = 8)) {
    end_idx <- min(i + 7, length(gene_vec))
    cat("   ", paste(gene_vec[i:end_idx], collapse = ", "), "\n")
  }

  cat("\n⏰ End time: ", format(Sys.time()), "\n")
}

if (!interactive()) {
  main()
}