#!/usr/bin/env Rscript
# Master Gene Downloader - Calls All Database Scripts
# Downloads genes from all available databases for a given phenotype

# Set CRAN mirror for non-interactive mode
if (!interactive() && is.null(getOption("repos"))) {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

# Required packages
required_packages <- c("dplyr", "readr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

#' Run an R script and check success based on file generation
#' @param script_path Path to the R script
#' @param phenotype The phenotype to search for
#' @return List with success status and details
run_database_script <- function(script_path, phenotype) {
  script_name <- basename(script_path)
  source_name <- gsub("\\.R$", "", script_name)

  cat("\n🔄 Running", source_name, "...\n")
  cat("   Command: Rscript", script_path, phenotype, "\n")

  clean_phenotype <- gsub("[^a-zA-Z0-9_-]", "_", phenotype)
  output_dir <- "AllPackagesGenes"

  # Correct suffix map matching all fixed scripts
  suffix_map <- list(
    "kegg"              = c("_kegg.csv",                   "_kegg_genes.csv"),
    "gene_ontology"     = c("_gene_ontology_full.csv",     "_gene_ontology_genes.csv"),
    "gwasrapidd"        = c("_gwasrapidd.csv",             "_gwasrapidd_genes.csv"),
    "reactome_pathways" = c("_reactome_pathways.csv",      "_reactome_pathways_genes.csv"),
    "gtex"              = c("_gtex_prioritized_genes.csv", "_gtex_genes.csv"),
    "pubmed_pubtator"   = c("_pubmed_pubtator.csv",        "_pubmed_genes.csv"),
    "clinvar"           = c("_clinvar.csv",                "_clinvar_genes.csv"),
    "hpo"               = c("_hpo.csv",                    "_hpo_genes.csv"),
    "omim"              = c("_omim.csv",                   "_omim_genes.csv"),
    "opentargets"       = c("_opentargets.csv",            "_opentargets_genes.csv"),
    "uniprot"           = c("_uniprot.csv",                "_uniprot_genes.csv"),
    "string_db"         = c("_string_db.csv",              "_string_db_genes.csv"),
    "disgenet"          = c("_disgenet.csv",               "_disgenet_genes.csv")
  )

  suffixes <- suffix_map[[source_name]]
  if (is.null(suffixes)) {
    suffixes <- c(paste0("_", source_name, ".csv"),
                  paste0("_", source_name, "_genes.csv"))
  }

  expected_files <- c(
    file.path(output_dir, paste0(clean_phenotype, suffixes[1])),
    file.path(output_dir, paste0(clean_phenotype, suffixes[2]))
  )

  # Check files before execution
  files_before <- sapply(expected_files, file.exists)

  # Run the script and capture output
  start_time <- Sys.time()

  tryCatch({
    result <- system2("Rscript", args = c(script_path, shQuote(phenotype)),
                  stdout = TRUE, stderr = TRUE, wait = TRUE)

    end_time <- Sys.time()
    runtime  <- as.numeric(difftime(end_time, start_time, units = "secs"))

    Sys.sleep(2)

    files_after   <- sapply(expected_files, file.exists)
    files_created <- any(files_after & !files_before)
    files_exist   <- any(files_after)
    success       <- files_created || files_exist

    exit_status <- attr(result, "status")
    if (is.null(exit_status)) exit_status <- 0

    # Initialise before conditional use — fixes scoping crash
    existing_files <- character(0)

    if (success) {
      if (files_created) {
        cat("   ✅", source_name, "completed successfully - files generated in",
            round(runtime, 1), "seconds\n")
      } else {
        cat("   ✅", source_name, "completed - output files found in",
            round(runtime, 1), "seconds\n")
      }

      existing_files <- expected_files[files_after]
      if (length(existing_files) > 0) {
        cat("   📄 Generated files:\n")
        for (file in existing_files) cat("      •", basename(file), "\n")
      }
    } else {
      cat("   ❌", source_name, "failed - no expected output files generated\n")
      cat("   ⚠️  Exit status:", exit_status, "\n")
      cat("   📁 Expected files:\n")
      for (file in expected_files) cat("      •", basename(file), "(missing)\n")
      if (length(result) > 0) {
        cat("   📝 Script output (last 3 lines):\n")
        for (line in tail(result, 3)) cat("      ", line, "\n")
      }
    }

    return(list(
      source        = source_name,
      success       = success,
      runtime       = runtime,
      exit_status   = exit_status,
      files_created = files_created,
      files_found   = existing_files,
      output        = result
    ))

  }, error = function(e) {
    end_time <- Sys.time()
    runtime  <- as.numeric(difftime(end_time, start_time, units = "secs"))

    cat("   ❌", source_name, "failed with error:", conditionMessage(e), "\n")

    Sys.sleep(1)
    files_after    <- sapply(expected_files, file.exists)
    files_created  <- any(files_after & !files_before)
    existing_files <- expected_files[files_after]

    if (files_created) {
      cat("   🔄 But files were generated despite error - marking as success\n")
      return(list(
        source        = source_name,
        success       = TRUE,
        runtime       = runtime,
        files_created = files_created,
        files_found   = existing_files,
        error         = conditionMessage(e)
      ))
    } else {
      return(list(
        source        = source_name,
        success       = FALSE,
        runtime       = runtime,
        files_created = FALSE,
        files_found   = character(0),
        error         = conditionMessage(e)
      ))
    }
  })
}

#' Count genes from output files
#' @param phenotype The phenotype searched for
#' @param source_name Name of the data source
#' @return List with gene counts and file info
count_genes_from_files <- function(phenotype, source_name) {
  clean_phenotype <- gsub("[^a-zA-Z0-9_-]", "_", phenotype)
  clean_phenotype <- gsub("_+",             "_", clean_phenotype)
  clean_phenotype <- gsub("^_|_$",          "",  clean_phenotype)
  output_dir <- "AllPackagesGenes"

  patterns <- list(
    "kegg"              = c(paste0(clean_phenotype, "_kegg.csv"),
                            paste0(clean_phenotype, "_kegg_genes.csv")),
    "gene_ontology"     = c(paste0(clean_phenotype, "_gene_ontology_full.csv"),
                            paste0(clean_phenotype, "_gene_ontology_genes.csv")),
    "clinvar"           = c(paste0(clean_phenotype, "_clinvar.csv"),
                            paste0(clean_phenotype, "_clinvar_genes.csv")),
    "hpo"               = c(paste0(clean_phenotype, "_hpo.csv"),
                            paste0(clean_phenotype, "_hpo_genes.csv")),
    "gwasrapidd"        = c(paste0(clean_phenotype, "_gwasrapidd.csv"),
                            paste0(clean_phenotype, "_gwasrapidd_genes.csv")),
    "gtex"              = c(paste0(clean_phenotype, "_gtex_prioritized_genes.csv"),
                            paste0(clean_phenotype, "_gtex_genes.csv")),
    "pubmed_pubtator"   = c(paste0(clean_phenotype, "_pubmed_pubtator.csv"),
                            paste0(clean_phenotype, "_pubmed_genes.csv")),
    "uniprot"           = c(paste0(clean_phenotype, "_uniprot.csv"),
                            paste0(clean_phenotype, "_uniprot_genes.csv")),
    "reactome_pathways" = c(paste0(clean_phenotype, "_reactome_pathways.csv"),
                            paste0(clean_phenotype, "_reactome_pathways_genes.csv")),
    "omim"              = c(paste0(clean_phenotype, "_omim.csv"),
                            paste0(clean_phenotype, "_omim_genes.csv")),
    "opentargets"       = c(paste0(clean_phenotype, "_opentargets.csv"),
                            paste0(clean_phenotype, "_opentargets_genes.csv")),
    "disgenet"          = c(paste0(clean_phenotype, "_disgenet.csv"),
                            paste0(clean_phenotype, "_disgenet_genes.csv")),
    "string_db"         = c(paste0(clean_phenotype, "_string_db.csv"),
                            paste0(clean_phenotype, "_string_db_genes.csv"))
  )

  file_patterns <- patterns[[source_name]]
  if (is.null(file_patterns)) {
    return(list(total_genes = 0, unique_genes = 0, files_found = FALSE))
  }

  total_genes  <- 0
  unique_genes <- 0
  files_found  <- FALSE

  main_file  <- file.path(output_dir, file_patterns[1])
  genes_file <- file.path(output_dir, file_patterns[2])

  if (file.exists(main_file)) {
    files_found <- TRUE
    tryCatch({
      data      <- read.csv(main_file, stringsAsFactors = FALSE)
      total_genes <- nrow(data)
      gene_cols <- c("Gene", "GeneSymbol", "Gene_Symbol", "gene_symbol", "symbol")
      gene_col  <- NULL
      for (col in gene_cols) {
        if (col %in% colnames(data)) { gene_col <- col; break }
      }
      if (!is.null(gene_col)) {
        unique_genes <- length(unique(
          data[[gene_col]][!is.na(data[[gene_col]]) & data[[gene_col]] != ""]
        ))
      }
    }, error = function(e) {
      cat("     ⚠️ Error reading", main_file, ":", conditionMessage(e), "\n")
    })
  }

  if (file.exists(genes_file)) {
    tryCatch({
      genes_data <- read.csv(genes_file, stringsAsFactors = FALSE)
      if ("Gene" %in% colnames(genes_data)) {
        unique_genes <- nrow(genes_data)
      }
    }, error = function(e) {
      cat("     ⚠️ Error reading", genes_file, ":", conditionMessage(e), "\n")
    })
  }

  return(list(total_genes = total_genes, unique_genes = unique_genes,
              files_found = files_found))
}

#' Aggregate all gene results
#' @param phenotype The phenotype searched for
#' @return Summary of all results
aggregate_all_results <- function(phenotype) {
  cat("\n📊 AGGREGATING ALL RESULTS FOR:", toupper(phenotype), "\n")
  cat("═══════════════════════════════════════════════════════\n")

  output_dir      <- "AllPackagesGenes"
  clean_phenotype <- gsub("[^a-zA-Z0-9_-]", "_", phenotype)

  all_files <- list.files(output_dir,
                          pattern    = paste0("^", clean_phenotype, "_.*\\.csv$"),
                          full.names = TRUE)

  total_unique_genes <- character()
  source_summary     <- data.frame()

  # Updated source mapping matching actual output filenames
  source_mapping <- list(
    "kegg"              = "KEGG Database",
    "gene_ontology"     = "Gene Ontology",
    "clinvar"           = "ClinVar",
    "hpo"               = "Human Phenotype Ontology",
    "gwasrapidd"        = "GWAS Catalog",
    "gtex"              = "GTEx",
    "pubmed"            = "PubMed",
    "uniprot"           = "UniProt",
    "reactome_pathways" = "Reactome",
    "omim"              = "OMIM",
    "opentargets"       = "Open Targets",
    "disgenet"          = "DisGeNET",
    "string_db"         = "STRING-DB"
  )

  for (source_key in names(source_mapping)) {
    source_name <- source_mapping[[source_key]]

    source_files <- all_files[grepl(paste0("_", source_key, "\\.csv$"), all_files)]
    genes_files  <- all_files[grepl(paste0("_", source_key, "_genes\\.csv$"), all_files)]

    # Exclude combined summary files
    source_files <- source_files[!grepl("ALL_SOURCES|SOURCES_SUMMARY", source_files)]
    genes_files  <- genes_files[!grepl("ALL_SOURCES|SOURCES_SUMMARY", genes_files)]

    total_associations <- 0
    unique_genes_count <- 0
    source_genes       <- character()

    if (length(source_files) > 0) {
      main_file <- source_files[1]
      tryCatch({
        data              <- read.csv(main_file, stringsAsFactors = FALSE)
        total_associations <- nrow(data)
        gene_cols <- c("Gene", "GeneSymbol", "Gene_Symbol", "gene_symbol", "symbol")
        gene_col  <- NULL
        for (col in gene_cols) {
          if (col %in% colnames(data)) { gene_col <- col; break }
        }
        if (!is.null(gene_col)) {
          valid_genes        <- data[[gene_col]][!is.na(data[[gene_col]]) & data[[gene_col]] != ""]
          source_genes       <- unique(valid_genes)
          unique_genes_count <- length(source_genes)
        }
      }, error = function(e) {
        cat("   ⚠️ Error reading", basename(main_file), "\n")
      })
    }

    if (length(genes_files) > 0) {
      genes_file <- genes_files[1]
      tryCatch({
        genes_data <- read.csv(genes_file, stringsAsFactors = FALSE)
        if ("Gene" %in% colnames(genes_data)) {
          source_genes       <- unique(genes_data$Gene[
            !is.na(genes_data$Gene) & genes_data$Gene != ""])
          unique_genes_count <- length(source_genes)
        }
      }, error = function(e) {
        cat("   ⚠️ Error reading", basename(genes_file), "\n")
      })
    }

    if (unique_genes_count > 0) {
      source_summary <- rbind(source_summary, data.frame(
        Source             = source_name,
        Total_Associations = total_associations,
        Unique_Genes       = unique_genes_count,
        stringsAsFactors   = FALSE
      ))
      total_unique_genes <- c(total_unique_genes, source_genes)
    }
  }

  all_unique_genes <- unique(total_unique_genes)
  total_sources    <- nrow(source_summary)

  cat("📋 RESULTS SUMMARY:\n")
  cat("═══════════════════\n")

  if (nrow(source_summary) > 0) {
    source_summary <- source_summary[order(-source_summary$Unique_Genes), ]

    for (i in 1:nrow(source_summary)) {
      row <- source_summary[i, ]
      cat(sprintf("📚 %-25s: %4d unique genes (%d total associations)\n",
                  row$Source, row$Unique_Genes, row$Total_Associations))
    }

    cat("\n🎯 OVERALL STATISTICS:\n")
    cat("═══════════════════════\n")
    cat(sprintf("🔬 Total data sources used: %d\n", total_sources))
    cat(sprintf("🧬 Total unique genes found: %d\n", length(all_unique_genes)))
    cat(sprintf("📊 Total associations: %d\n", sum(source_summary$Total_Associations)))
    cat(sprintf("📁 Results directory: %s/\n", output_dir))

    if (length(all_unique_genes) > 0) {
      combined_genes_file <- file.path(output_dir,
                                       paste0(clean_phenotype, "_ALL_SOURCES_GENES.csv"))
      write.csv(data.frame(Gene = sort(all_unique_genes)),
                combined_genes_file, row.names = FALSE)
      cat(sprintf("💾 Combined gene list saved: %s\n", combined_genes_file))
    }

    summary_file <- file.path(output_dir,
                              paste0(clean_phenotype, "_SOURCES_SUMMARY.csv"))
    write.csv(source_summary, summary_file, row.names = FALSE)
    cat(sprintf("📊 Source summary saved: %s\n", summary_file))

  } else {
    cat("❌ No genes found from any source\n")
  }

  return(list(
    total_sources      = total_sources,
    total_unique_genes = length(all_unique_genes),
    source_summary     = source_summary
  ))
}

#' Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("🧬 COMPREHENSIVE GENE DOWNLOADER\n")
    cat("═══════════════════════════════════════════\n")
    cat("Downloads genes from ALL available databases\n\n")
    cat("Usage: Rscript download_genes.R <phenotype> [--force]\n\n")
    cat("📋 EXAMPLES:\n")
    cat("   Rscript download_genes.R migraine\n")
    cat("   Rscript download_genes.R \"type 2 diabetes\"\n")
    cat("   Rscript download_genes.R \"heart disease\"\n")
    cat("   Rscript download_genes.R cancer --force\n\n")
    cat("🔧 OPTIONS:\n")
    cat("   --force    Re-run all scripts even if output files exist\n\n")
    cat("🗃️  DATABASES INCLUDED:\n")
    cat("   1.  kegg.R                 - KEGG Database\n")
    cat("       → {phenotype}_kegg.csv, {phenotype}_kegg_genes.csv\n")
    cat("   2.  gene_ontology.R        - Gene Ontology\n")
    cat("       → {phenotype}_gene_ontology_full.csv, {phenotype}_gene_ontology_genes.csv\n")
    cat("   3.  clinvar.R              - ClinVar\n")
    cat("       → {phenotype}_clinvar.csv, {phenotype}_clinvar_genes.csv\n")
    cat("   4.  hpo.R                  - Human Phenotype Ontology\n")
    cat("       → {phenotype}_hpo.csv, {phenotype}_hpo_genes.csv\n")
    cat("   5.  gwasrapidd.R           - GWAS Catalog\n")
    cat("       → {phenotype}_gwasrapidd.csv, {phenotype}_gwasrapidd_genes.csv\n")
    cat("   6.  gtex.R                 - GTEx\n")
    cat("       → {phenotype}_gtex_prioritized_genes.csv, {phenotype}_gtex_genes.csv\n")
    cat("   7.  pubmed_pubtator.R      - PubMed + PubTator3\n")
    cat("       → {phenotype}_pubmed_pubtator.csv, {phenotype}_pubmed_genes.csv\n")
    cat("   8.  uniprot.R              - UniProt\n")
    cat("       → {phenotype}_uniprot.csv, {phenotype}_uniprot_genes.csv\n")
    cat("   9.  reactome_pathways.R    - Reactome\n")
    cat("       → {phenotype}_reactome_pathways.csv, {phenotype}_reactome_pathways_genes.csv\n")
    cat("   10. omim.R                 - OMIM\n")
    cat("       → {phenotype}_omim.csv, {phenotype}_omim_genes.csv\n")
    cat("   11. opentargets.R          - Open Targets\n")
    cat("       → {phenotype}_opentargets.csv, {phenotype}_opentargets_genes.csv\n")
    cat("   12. disgenet.R             - DisGeNET\n")
    cat("       → {phenotype}_disgenet.csv, {phenotype}_disgenet_genes.csv\n")
    cat("   13. string_db.R            - STRING-DB\n")
    cat("       → {phenotype}_string_db.csv, {phenotype}_string_db_genes.csv\n\n")
    cat("📁 Output: AllPackagesGenes/ directory\n")
    quit(status = 1)
  }

  phenotype   <- args[1]
  force_rerun <- "--force" %in% args

  cat("🧬 COMPREHENSIVE GENE DOWNLOADER\n")
  cat("═══════════════════════════════════════════\n")
  cat("🎯 Target phenotype:", phenotype, "\n")
  cat("📅 Date:", format(Sys.Date()), "\n")
  cat("⏰ Start time:", format(Sys.time()), "\n")
  if (force_rerun) {
    cat("🔧 Mode: FORCE RERUN (ignoring existing files)\n")
  } else {
    cat("🔧 Mode: SMART SKIP (will skip if files exist)\n")
  }

  clean_phenotype <- gsub("[^a-zA-Z0-9_-]", "_", phenotype)

  cat("\n📁 EXPECTED OUTPUT FILES IN AllPackagesGenes/:\n")
  cat("═══════════════════════════════════════════════════\n")
  cat("Each script creates 2 files: Full results + Genes only\n\n")

  expected_files <- list(
    "KEGG Database"          = c(paste0(clean_phenotype, "_kegg.csv"),
                                  paste0(clean_phenotype, "_kegg_genes.csv")),
    "Gene Ontology"          = c(paste0(clean_phenotype, "_gene_ontology_full.csv"),
                                  paste0(clean_phenotype, "_gene_ontology_genes.csv")),
    "ClinVar"                = c(paste0(clean_phenotype, "_clinvar.csv"),
                                  paste0(clean_phenotype, "_clinvar_genes.csv")),
    "HPO"                    = c(paste0(clean_phenotype, "_hpo.csv"),
                                  paste0(clean_phenotype, "_hpo_genes.csv")),
    "GWAS Catalog"           = c(paste0(clean_phenotype, "_gwasrapidd.csv"),
                                  paste0(clean_phenotype, "_gwasrapidd_genes.csv")),
    "GTEx"                   = c(paste0(clean_phenotype, "_gtex_prioritized_genes.csv"),
                                  paste0(clean_phenotype, "_gtex_genes.csv")),
    "PubMed"                 = c(paste0(clean_phenotype, "_pubmed_pubtator.csv"),
                                  paste0(clean_phenotype, "_pubmed_genes.csv")),
    "UniProt"                = c(paste0(clean_phenotype, "_uniprot.csv"),
                                  paste0(clean_phenotype, "_uniprot_genes.csv")),
    "Reactome"               = c(paste0(clean_phenotype, "_reactome_pathways.csv"),
                                  paste0(clean_phenotype, "_reactome_pathways_genes.csv")),
    "OMIM"                   = c(paste0(clean_phenotype, "_omim.csv"),
                                  paste0(clean_phenotype, "_omim_genes.csv")),
    "Open Targets"           = c(paste0(clean_phenotype, "_opentargets.csv"),
                                  paste0(clean_phenotype, "_opentargets_genes.csv")),
    "DisGeNET"               = c(paste0(clean_phenotype, "_disgenet.csv"),
                                  paste0(clean_phenotype, "_disgenet_genes.csv")),
    "STRING-DB"              = c(paste0(clean_phenotype, "_string_db.csv"),
                                  paste0(clean_phenotype, "_string_db_genes.csv"))
  )

  for (db_name in names(expected_files)) {
    files <- expected_files[[db_name]]
    cat(sprintf("📚 %-15s: %s\n", db_name, files[1]))
    cat(sprintf("   %-15s  %s\n", "", files[2]))
  }

  cat("\n🎯 COMBINED OUTPUT FILES:\n")
  cat(sprintf("💾 %s_ALL_SOURCES_GENES.csv (all unique genes)\n", clean_phenotype))
  cat(sprintf("📊 %s_SOURCES_SUMMARY.csv (source statistics)\n", clean_phenotype))

  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("📁 Created output directory:", output_dir, "\n")
  }

  # All scripts — string.R removed, pubmed_pubtator.R is the pubmed script
  scripts_info <- list(
    "kegg.R"              = "KEGG Database - Pathways and gene associations",
    "gene_ontology.R"     = "Gene Ontology - Functional annotations",
    "clinvar.R"           = "ClinVar - Clinical variants and pathogenicity",
    "hpo.R"               = "Human Phenotype Ontology - Phenotype associations",
    "gwasrapidd.R"        = "GWAS Catalog - Genome-wide association studies",
    "gtex.R"              = "GTEx - Gene expression data",
    "pubmed_pubtator.R"   = "PubMed + PubTator3 - Literature-based gene associations",
    "uniprot.R"           = "UniProt - Protein database",
    "reactome_pathways.R" = "Reactome - Biological pathways",
    "omim.R"              = "OMIM - Mendelian inheritance database",
    "opentargets.R"       = "Open Targets - Drug target platform",
    "disgenet.R"          = "DisGeNET - Gene-disease associations",
    "string_db.R"         = "STRING-DB - Alternative protein interactions"
  )

  scripts <- names(scripts_info)

  cat("\n🚀 RUNNING", length(scripts), "DATABASE SCRIPTS\n")
  cat("═══════════════════════════════════════════\n")
  cat("📋 R SCRIPTS TO EXECUTE:\n")
  for (i in seq_along(scripts)) {
    cat(sprintf("   %2d. %-25s - %s\n", i, scripts[i], scripts_info[[scripts[i]]]))
  }
  cat("\n")

  start_time         <- Sys.time()
  results            <- list()
  successful_scripts <- 0
  failed_scripts     <- 0

  for (i in seq_along(scripts)) {
    script      <- scripts[i]
    script_path <- file.path(getwd(), script)
    description <- scripts_info[[script]]

    cat(sprintf("\n[%d/%d] 📄 CHECKING: %s\n", i, length(scripts), script))
    cat(sprintf("     📋 %s\n", description))

    script_base <- gsub("\\.R$", "", script)

    # Correct skip-check suffix map
    skip_map <- list(
      "kegg"              = "_kegg_genes.csv",
      "gene_ontology"     = "_gene_ontology_genes.csv",
      "gwasrapidd"        = "_gwasrapidd_genes.csv",
      "reactome_pathways" = "_reactome_pathways_genes.csv",
      "gtex"              = "_gtex_genes.csv",
      "pubmed_pubtator"   = "_pubmed_genes.csv"
    )
    genes_suffix   <- skip_map[[script_base]]
    if (is.null(genes_suffix)) {
      genes_suffix <- paste0("_", script_base, "_genes.csv")
    }

    # Also define full output suffix for skip check
    full_skip_map <- list(
      "kegg"              = "_kegg.csv",
      "gene_ontology"     = "_gene_ontology_full.csv",
      "gwasrapidd"        = "_gwasrapidd.csv",
      "reactome_pathways" = "_reactome_pathways.csv",
      "gtex"              = "_gtex_prioritized_genes.csv",
      "pubmed_pubtator"   = "_pubmed_pubtator.csv"
    )
    output_suffix   <- full_skip_map[[script_base]]
    if (is.null(output_suffix)) {
      output_suffix <- paste0("_", script_base, ".csv")
    }

    expected_output <- file.path(output_dir, paste0(clean_phenotype, output_suffix))
    expected_genes  <- file.path(output_dir, paste0(clean_phenotype, genes_suffix))

    files_exist        <- file.exists(expected_output) || file.exists(expected_genes)
    files_have_content <- FALSE

    if (files_exist) {
      if (file.exists(expected_output)) {
        tryCatch({
          test_data          <- read.csv(expected_output, nrows = 5)
          files_have_content <- nrow(test_data) > 0
        }, error = function(e) { files_have_content <- FALSE })
      }
      if (!files_have_content && file.exists(expected_genes)) {
        tryCatch({
          test_data          <- read.csv(expected_genes, nrows = 5)
          files_have_content <- nrow(test_data) > 0
        }, error = function(e) { files_have_content <- FALSE })
      }
    }

    if (!force_rerun && files_exist && files_have_content) {
      cat("     ✅ Output files already exist with content - SKIPPING\n")
      if (file.exists(expected_output)) cat("       📄", basename(expected_output), "\n")
      if (file.exists(expected_genes))  cat("       🧬", basename(expected_genes), "\n")
      successful_scripts          <- successful_scripts + 1
      results[[script]]           <- list(source  = script_base, success = TRUE,
                                          runtime = 0,
                                          output  = "Skipped - files exist")
    } else if (file.exists(script_path)) {
      if (force_rerun) {
        cat("     🔄 EXECUTING (FORCE mode)\n")
      } else {
        cat("     🔄 EXECUTING (output files missing or empty)\n")
        if (!file.exists(expected_output)) cat("       ❌ Missing:", basename(expected_output), "\n")
        if (!file.exists(expected_genes))  cat("       ❌ Missing:", basename(expected_genes), "\n")
        if (files_exist && !files_have_content) cat("       ⚠️  Files exist but appear empty\n")
      }

      result            <- run_database_script(script_path, phenotype)
      results[[script]] <- result

      if (result$success) {
        successful_scripts <- successful_scripts + 1
      } else {
        failed_scripts <- failed_scripts + 1
      }
    } else {
      cat("❌", script, "not found in current directory\n")
      failed_scripts <- failed_scripts + 1
    }

    Sys.sleep(1)
  }

  end_time      <- Sys.time()
  total_runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))

  cat("\n✅ SCRIPT EXECUTION COMPLETE\n")
  cat("═══════════════════════════════════════════\n")
  cat(sprintf("⏱️  Total runtime: %.1f minutes\n", total_runtime))
  cat(sprintf("✅ Successful scripts: %d/%d\n", successful_scripts, length(scripts)))
  cat(sprintf("❌ Failed scripts: %d/%d\n", failed_scripts, length(scripts)))

  if (failed_scripts > 0) {
    cat("\n❌ FAILED SCRIPTS:\n")
    for (script in names(results)) {
      if (!isTRUE(results[[script]]$success)) cat("   •", script, "\n")
    }
  }

  cat("\n⏳ Waiting for files to be written...\n")
  Sys.sleep(3)

  summary <- aggregate_all_results(phenotype)

  cat("\n🎉 COMPREHENSIVE GENE DOWNLOAD COMPLETE!\n")
  cat("═══════════════════════════════════════════\n")
  cat(sprintf("🎯 Phenotype: %s\n", phenotype))
  cat(sprintf("🔬 Sources used: %d\n", summary$total_sources))
  cat(sprintf("🧬 Total unique genes: %d\n", summary$total_unique_genes))
  cat(sprintf("📁 Results in: %s/\n", output_dir))
  cat(sprintf("⏱️  Total time: %.1f minutes\n", total_runtime))

  cat("\n💡 Next steps:\n")
  cat("   • Check AllPackagesGenes/ directory for all results\n")
  cat("   • Use *_ALL_SOURCES_GENES.csv for combined gene list\n")
  cat("   • Use *_SOURCES_SUMMARY.csv for source statistics\n")
  cat("   • Individual database results in separate CSV files\n")
}

# Run if called as script
if (!interactive()) {
  main()
}