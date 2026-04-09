#!/usr/bin/env Rscript
# Reactome pathway gene downloader using ReactomePA package
# Downloads genes from Reactome pathways related to phenotype

bioc_packages <- c("ReactomePA", "reactome.db", "org.Hs.eg.db", "AnnotationDbi")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    if (!require(BiocManager, quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

required_packages <- c("dplyr", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

download_reactome_genes <- function(phenotype) {
  cat("🔍 Searching Reactome pathways for:", phenotype, "\n")
  
  tryCatch({
    # Get all human Reactome pathways
    pathway_ids <- AnnotationDbi::keys(reactome.db, keytype = "PATHID")
    pathway_info <- AnnotationDbi::select(reactome.db,
                                          keys     = pathway_ids,
                                          columns  = c("PATHID", "PATHNAME"),
                                          keytype  = "PATHID")
    
    # Keep human pathways only
    pathway_info <- pathway_info[grepl("^Homo sapiens:", pathway_info$PATHNAME), ]
    cat("   Total human pathways available:", nrow(pathway_info), "\n")
    
    # Build search patterns from phenotype string only — no hardcoding
    phenotype_lower <- tolower(phenotype)
    phenotype_words <- unlist(strsplit(phenotype_lower, "\\s+"))
    
    search_patterns <- c(
      paste0("\\b", phenotype_lower, "\\b"),
      paste0("\\b", phenotype_lower, " pathway\\b"),
      paste0("\\b", phenotype_lower, " signaling\\b"),
      paste0("\\b", phenotype_lower, " disease\\b"),
      paste0("\\b", phenotype_lower, " disorder\\b")
    )
    
    # Add individual word patterns only for multi-word phenotypes
    if (length(phenotype_words) > 1) {
      word_patterns <- paste0("\\b", phenotype_words[nchar(phenotype_words) > 3], "\\b")
      search_patterns <- c(search_patterns, word_patterns)
    }
    
    search_patterns <- unique(search_patterns[search_patterns != "" & !is.na(search_patterns)])
    cat("   Searching with", length(search_patterns), "patterns\n")
    
    # Find matching pathways
    relevant_pathways <- character()
    
    for (pattern in search_patterns) {
      matching <- pathway_info[grepl(pattern, tolower(pathway_info$PATHNAME), perl = TRUE), ]
      if (nrow(matching) > 0) {
        relevant_pathways <- c(relevant_pathways, matching$PATHID)
        cat("   Pattern '", pattern, "' matched", nrow(matching), "pathways\n", sep = "")
      }
    }
    
    relevant_pathways <- unique(relevant_pathways)
    cat("   Found", length(relevant_pathways), "relevant Reactome pathways\n")
    
    if (length(relevant_pathways) == 0) {
      cat("   No relevant Reactome pathways found for:", phenotype, "\n")
      return(data.frame())
    }
    
    # Cap at 50 to avoid overwhelming results
    if (length(relevant_pathways) > 50) {
      cat("   Limiting to top 50 pathways\n")
      relevant_pathways <- relevant_pathways[1:50]
    }
    
    # Retrieve genes per pathway
    all_genes        <- data.frame()
    successful_pathways <- 0
    
    cat("   Retrieving genes from pathways...\n")
    
    for (pathway_id in relevant_pathways) {
      tryCatch({
        entrez_ids <- AnnotationDbi::get(pathway_id, reactome.db::reactomePATHID2EXTID)
        
        if (length(entrez_ids) > 0) {
          gene_symbols <- AnnotationDbi::select(org.Hs.eg.db,
                                                keys    = as.character(entrez_ids),
                                                columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                                                keytype = "ENTREZID")
          
          if (nrow(gene_symbols) > 0) {
            pathway_name <- pathway_info$PATHNAME[pathway_info$PATHID == pathway_id][1]
            
            gene_df <- data.frame(
              Gene         = gene_symbols$SYMBOL,
              Gene_Name    = gene_symbols$GENENAME,
              Entrez_ID    = gene_symbols$ENTREZID,
              Pathway_ID   = pathway_id,
              Pathway_Name = pathway_name,
              Source       = "Reactome_Pathways",
              Phenotype    = phenotype,
              stringsAsFactors = FALSE
            )
            
            gene_df <- gene_df[!is.na(gene_df$Gene) & gene_df$Gene != "", ]
            
            if (nrow(gene_df) > 0) {
              all_genes <- rbind(all_genes, gene_df)
              successful_pathways <- successful_pathways + 1
            }
          }
        }
      }, error = function(e) {
        cat("     ⚠️  Pathway", pathway_id, "error:", conditionMessage(e), "\n")
      })
    }
    
    cat("   Successfully processed", successful_pathways, "of",
        length(relevant_pathways), "pathways\n")
    
    if (nrow(all_genes) > 0) {
      all_genes <- all_genes[!duplicated(all_genes$Gene), ]
      cat("   ✅ Total unique genes from Reactome:", nrow(all_genes), "\n")
      return(all_genes)
    }
    
    cat("   No genes found from Reactome pathways\n")
    return(data.frame())
    
  }, error = function(e) {
    cat("❌ Reactome error:", conditionMessage(e), "\n")
    return(data.frame())
  })
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Reactome Pathway Gene Downloader\n")
    cat("Usage:   Rscript reactome_pathways.R <phenotype>\n")
    cat("Example: Rscript reactome_pathways.R migraine\n")
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
  output_file     <- file.path(output_dir, paste0(clean_phenotype, "_reactome_pathways.csv"))
  genes_only_file <- file.path(output_dir, paste0(clean_phenotype, "_reactome_pathways_genes.csv"))
  
  cat("🎯 Phenotype:        ", phenotype, "\n")
  cat("📁 Output directory: ", output_dir, "\n")
  cat("⏰ Start time:       ", format(Sys.time()), "\n\n")
  
  results <- download_reactome_genes(phenotype)
  
  if (nrow(results) > 0) {
    write.csv(results, output_file, row.names = FALSE)
    
    if ("Gene" %in% colnames(results)) {
      genes_only <- data.frame(Gene = unique(results$Gene))
      write.csv(genes_only, genes_only_file, row.names = FALSE)
    }
    
    cat("\n✅ SUCCESS!\n")
    cat("📊 Total records:   ", nrow(results), "\n")
    cat("🧬 Unique genes:    ", length(unique(results$Gene)), "\n")
    cat("💾 Full results:    ", output_file, "\n")
    cat("🧬 Genes-only file: ", genes_only_file, "\n\n")
    
    # Pathway summary
    if ("Pathway_Name" %in% colnames(results)) {
      cat("Top pathways by gene count:\n")
      pathway_counts <- sort(table(results$Pathway_Name), decreasing = TRUE)
      for (i in seq_len(min(5, length(pathway_counts)))) {
        cat(sprintf("  (%d genes) %s\n", pathway_counts[i],
                    substr(names(pathway_counts)[i], 1, 70)))
      }
    }
    
    unique_genes <- sort(unique(results$Gene))
    cat("\n🧬 Genes found for", phenotype, ":\n")
    for (i in seq(1, length(unique_genes), by = 8)) {
      end_idx <- min(i + 7, length(unique_genes))
      cat("   ", paste(unique_genes[i:end_idx], collapse = ", "), "\n")
    }
    
  } else {
    cat("❌ No genes found for:", phenotype, "\n")
    cat("Reactome may not have pathways matching this phenotype term.\n")
  }
  
  cat("\n⏰ End time: ", format(Sys.time()), "\n")
}

if (!interactive()) {
  main()
}