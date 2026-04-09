#!/usr/bin/env Rscript

if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}

library(httr)
library(jsonlite)

`%||%` <- function(a, b) if (!is.null(a)) a else b

download_uniprot_genes <- function(phenotype) {
  message("🔍 Searching UniProt for phenotype: ", phenotype)
  
  base_url <- "https://rest.uniprot.org/uniprotkb/search"
  
  search_queries <- c(
    paste0("(keyword:\"", phenotype, "\") AND (organism_id:9606) AND (reviewed:true)"),
    paste0("(cc_disease:\"", phenotype, "\") AND (organism_id:9606) AND (reviewed:true)"),
    paste0("(ft_mutagen:\"", phenotype, "\") AND (organism_id:9606) AND (reviewed:true)"),
    paste0("(\"", phenotype, "\") AND (organism_id:9606) AND (reviewed:true)")
  )
  
  all_genes_list <- list()
  
  for (query in search_queries) {
    message("   Trying query: ", query)
    
    tryCatch({
      response <- GET(
        base_url,
        query = list(query = query, format = "json", size = 500),
        add_headers(`User-Agent` = "R script for UniProt query - your_email@example.com"),
        timeout(30)
      )
      
      if (status_code(response) != 200) {
        warning("Request failed for query: ", query, " Status: ", status_code(response))
        next
      }
      
      res_text <- content(response, as = "text", encoding = "UTF-8")
      data <- fromJSON(res_text, simplifyVector = FALSE)
      
      if (!is.null(data$results) && length(data$results) > 0) {
        message("   Found ", length(data$results), " entries for this query")
        
        for (entry in data$results) {
          gene_names <- c()
          
          if (!is.null(entry$genes)) {
            for (gene in entry$genes) {
              if (!is.null(gene$geneName$value)) {
                gene_names <- c(gene_names, gene$geneName$value)
              }
              if (!is.null(gene$synonyms)) {
                for (syn in gene$synonyms) {
                  if (!is.null(syn$value)) {
                    gene_names <- c(gene_names, syn$value)
                  }
                }
              }
            }
          }
          
          protein_name <- ""
          if (!is.null(entry$proteinDescription$recommendedName$fullName$value)) {
            protein_name <- entry$proteinDescription$recommendedName$fullName$value
          }
          
          unique_genes <- unique(gene_names)
          if (length(unique_genes) > 0) {
            df <- data.frame(
              Gene        = unique_genes,
              UniProt_ID  = entry$primaryAccession %||% NA,
              Protein_Name = protein_name,
              Source      = "UniProt",
              Query_Type  = query,
              stringsAsFactors = FALSE
            )
            all_genes_list[[length(all_genes_list) + 1]] <- df
          }
        }
      } else {
        message("   No entries returned for this query")
      }
      
    }, error = function(e) {
      message("   ❌ Error on query: ", query, " — ", conditionMessage(e))
    })
    
    Sys.sleep(1)
  }
  
  if (length(all_genes_list) == 0) {
    message("   No genes found for phenotype: ", phenotype)
    return(data.frame())
  }
  
  all_genes <- do.call(rbind, all_genes_list)
  all_genes <- all_genes[!duplicated(all_genes$Gene), ]
  all_genes <- all_genes[order(all_genes$Gene), ]
  
  message("   ✅ Total unique genes retrieved from UniProt: ", nrow(all_genes))
  return(all_genes)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    message("UniProt Gene Downloader")
    message("Usage:   Rscript uniprot.R <phenotype>")
    message("Example: Rscript uniprot.R migraine")
    quit(status = 1)
  }
  
  phenotype <- args[1]
  
  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  clean_phenotype <- gsub("[^a-zA-Z0-9_-]", "_", phenotype)
clean_phenotype <- gsub("_+",             "_", clean_phenotype)
clean_phenotype <- gsub("^_|_$",          "",  clean_phenotype)
  
  output_file     <- ifelse(length(args) > 1,
                            args[2],
                            file.path(output_dir, paste0(clean_phenotype, "_uniprot.csv")))
  genes_only_file <- file.path(output_dir, paste0(clean_phenotype, "_uniprot_genes.csv"))
  
  cat("🎯 Phenotype:        ", phenotype, "\n")
  cat("📁 Output directory: ", output_dir, "\n")
  cat("⏰ Start time:       ", format(Sys.time()), "\n\n")
  
  results <- download_uniprot_genes(phenotype)
  
  if (nrow(results) > 0) {
    write.csv(results, output_file, row.names = FALSE)
    
    if ("Gene" %in% colnames(results)) {
      genes_only <- data.frame(Gene = unique(results$Gene))
      write.csv(genes_only, genes_only_file, row.names = FALSE)
    }
    
    cat("\n✅ SUCCESS!\n")
    cat("📊 Total records:    ", nrow(results), "\n")
    cat("🧬 Unique genes:     ", length(unique(results$Gene)), "\n")
    cat("💾 Full results:     ", output_file, "\n")
    cat("🧬 Genes-only file:  ", genes_only_file, "\n\n")
    
    unique_genes <- sort(unique(results$Gene))
    cat("🧬 Genes found for", phenotype, ":\n")
    for (i in seq(1, length(unique_genes), by = 8)) {
      end_idx <- min(i + 7, length(unique_genes))
      cat("   ", paste(unique_genes[i:end_idx], collapse = ", "), "\n")
    }
    
  } else {
    cat("❌ No genes found for:", phenotype, "\n")
    cat("Try a different phenotype term or check your internet connection.\n")
  }
  
  cat("\n⏰ End time: ", format(Sys.time()), "\n")
}

if (!interactive()) {
  main()
}