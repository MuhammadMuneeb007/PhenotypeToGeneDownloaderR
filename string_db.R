#!/usr/bin/env Rscript
# STRING database gene downloader
# Downloads protein-protein interaction genes from STRING database

required_packages <- c("httr", "jsonlite")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

download_string_genes <- function(phenotype) {
  cat("🔍 Searching STRING database for:", phenotype, "\n")
  
  tryCatch({
    # Step 1: Resolve phenotype string to STRING protein IDs
    cat("   Resolving phenotype to STRING protein IDs...\n")
    
    resolve_url <- "https://string-db.org/api/json/get_string_ids"
    resolve_response <- GET(
      resolve_url,
      query = list(
        identifiers     = phenotype,
        species         = 9606,
        limit           = 10,
        caller_identity = "PhenotypeToGeneDownloaderR"
      ),
      timeout(30)
    )
    
    if (status_code(resolve_response) != 200) {
      cat("   ⚠️  Could not resolve phenotype in STRING. Status:",
          status_code(resolve_response), "\n")
      return(data.frame())
    }
    
    resolved <- fromJSON(content(resolve_response, "text", encoding = "UTF-8"),
                         simplifyVector = TRUE)
    
    if (length(resolved) == 0 || nrow(resolved) == 0) {
      cat("   No STRING entries found for:", phenotype, "\n")
      return(data.frame())
    }
    
    seed_genes <- unique(resolved$preferredName)
    seed_genes <- seed_genes[!is.na(seed_genes) & seed_genes != ""]
    cat("   Resolved", length(seed_genes), "seed proteins:", paste(seed_genes, collapse = ", "), "\n")
    
    all_genes <- data.frame()
    
    # Step 2: Get interaction partners for each resolved protein
    for (seed_gene in seed_genes) {
      cat("   Fetching interaction partners for:", seed_gene, "\n")
      
      interaction_url <- "https://string-db.org/api/json/interaction_partners"
      response <- GET(
        interaction_url,
        query = list(
          identifiers     = seed_gene,
          species         = 9606,
          limit           = 50,
          required_score  = 400,
          caller_identity = "PhenotypeToGeneDownloaderR"
        ),
        timeout(30)
      )
      
      if (status_code(response) != 200) {
        cat("   ⚠️  Request failed for", seed_gene, "— skipping\n")
        Sys.sleep(1)
        next
      }
      
      interactions <- fromJSON(content(response, "text", encoding = "UTF-8"),
                               simplifyVector = TRUE)
      
      if (length(interactions) == 0 || nrow(interactions) == 0) {
        cat("   No interactions found for:", seed_gene, "\n")
        Sys.sleep(1)
        next
      }
      
      valid <- !is.na(interactions$preferredName_B) & interactions$preferredName_B != ""
      if (sum(valid) > 0) {
        batch <- data.frame(
          Gene             = interactions$preferredName_B[valid],
          Seed_Gene        = seed_gene,
          Interaction_Score = interactions$score[valid],
          STRING_ID        = interactions$stringId_B[valid],
          Source           = "STRING",
          stringsAsFactors = FALSE
        )
        all_genes <- rbind(all_genes, batch)
        cat("   Found", nrow(batch), "interaction partners\n")
      }
      
      Sys.sleep(1)
    }
    
    if (nrow(all_genes) == 0) {
      cat("   No interaction partners found\n")
      return(data.frame())
    }
    
    all_genes <- all_genes[!duplicated(all_genes$Gene), ]
    all_genes <- all_genes[order(-all_genes$Interaction_Score), ]
    cat("   ✅ Total unique genes from STRING:", nrow(all_genes), "\n")
    return(all_genes)
    
  }, error = function(e) {
    cat("❌ STRING error:", conditionMessage(e), "\n")
    return(data.frame())
  })
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("STRING Gene Downloader\n")
    cat("Usage:   Rscript string_db.R <phenotype>\n")
    cat("Example: Rscript string_db.R migraine\n")
    quit(status = 1)
  }
  
  phenotype <- args[1]
  
  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  clean_phenotype     <- gsub("[^\\w\\s-]", "", phenotype, perl = TRUE)
  output_file         <- file.path(output_dir, paste0(clean_phenotype, "_string_db.csv"))
  genes_only_file     <- file.path(output_dir, paste0(clean_phenotype, "_string_db_genes.csv"))
  
  cat("🎯 Phenotype:        ", phenotype, "\n")
  cat("📁 Output directory: ", output_dir, "\n")
  cat("⏰ Start time:       ", format(Sys.time()), "\n\n")
  
  results <- download_string_genes(phenotype)
  
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
    
    unique_genes <- sort(unique(results$Gene))
    cat("🧬 Genes found for", phenotype, ":\n")
    for (i in seq(1, length(unique_genes), by = 8)) {
      end_idx <- min(i + 7, length(unique_genes))
      cat("   ", paste(unique_genes[i:end_idx], collapse = ", "), "\n")
    }
    
  } else {
    cat("❌ No genes found for:", phenotype, "\n")
    cat("STRING may not have a direct entry for this phenotype term.\n")
  }
  
  cat("\n⏰ End time: ", format(Sys.time()), "\n")
}

if (!interactive()) {
  main()
}