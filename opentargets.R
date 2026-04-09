#!/usr/bin/env Rscript
# Open Targets Platform gene downloader
# Downloads genes from Open Targets Platform using GraphQL API

required_packages <- c("httr", "jsonlite")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# GraphQL endpoint
ot_url <- "https://api.platform.opentargets.org/api/v4/graphql"

run_query <- function(query, variables = NULL) {
  tryCatch({
    response <- POST(
      ot_url,
      body   = list(query = query, variables = variables),
      encode = "json",
      content_type_json(),
      timeout(60)
    )
    
    if (response$status_code != 200) {
      cat("⚠️  GraphQL request failed with status:", response$status_code, "\n")
      return(NULL)
    }
    
    content(response, as = "parsed", type = "application/json")
    
  }, error = function(e) {
    cat("❌ Query error:", conditionMessage(e), "\n")
    return(NULL)
  })
}

get_mondo_id <- function(phenotype_name) {
  cat("🔍 Searching Open Targets for:", phenotype_name, "\n")
  
  query <- '
    query searchDisease($term: String!) {
      search(queryString: $term, entityNames: ["disease"]) {
        hits {
          id
          name
        }
      }
    }
  '
  
  result <- run_query(query, list(term = phenotype_name))
  
  if (is.null(result)) return(NULL)
  
  hits <- result$data$search$hits
  
  if (is.null(hits) || length(hits) == 0) {
    cat("   No disease match found for:", phenotype_name, "\n")
    return(NULL)
  }
  
  mondo_id <- hits[[1]]$id
  cat("   Matched disease:", hits[[1]]$name, "(", mondo_id, ")\n")
  return(mondo_id)
}

get_genes <- function(mondo_id, page_size = 1000) {
  cat("🧬 Retrieving associated genes...\n")
  
  all_rows  <- list()
  page_index <- 0
  
  repeat {
    query <- '
      query diseaseAssociations($efoId: String!, $index: Int!, $size: Int!) {
        disease(efoId: $efoId) {
          associatedTargets(page: { index: $index, size: $size }) {
            rows {
              target {
                id
                approvedSymbol
                approvedName
              }
              score
            }
          }
        }
      }
    '
    
    result <- run_query(query, list(
      efoId = mondo_id,
      index = page_index,
      size  = page_size
    ))
    
    if (is.null(result)) break
    
    rows <- result$data$disease$associatedTargets$rows
    
    if (is.null(rows) || length(rows) == 0) break
    
    all_rows <- c(all_rows, rows)
    cat("   Page", page_index, "— retrieved", length(rows), "associations\n")
    
    # If fewer rows returned than page_size, we are on the last page
    if (length(rows) < page_size) break
    
    page_index <- page_index + 1
    Sys.sleep(0.5)
  }
  
  if (length(all_rows) == 0) {
    cat("   No genes found for:", mondo_id, "\n")
    return(data.frame())
  }
  
  data.frame(
    Gene        = sapply(all_rows, function(x) x$target$approvedSymbol),
    TargetID    = sapply(all_rows, function(x) x$target$id),
    GeneName    = sapply(all_rows, function(x) x$target$approvedName),
    Score       = sapply(all_rows, function(x) x$score),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Open Targets Gene Downloader\n")
    cat("Usage:   Rscript opentargets.R <phenotype>\n")
    cat("Example: Rscript opentargets.R migraine\n")
    quit(status = 1)
  }
  
  phenotype_name <- args[1]
  
  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  clean_phenotype <- gsub("[^\\w\\s-]", "", phenotype_name, perl = TRUE)
  output_file     <- file.path(output_dir, paste0(clean_phenotype, "_opentargets.csv"))
  genes_only_file <- file.path(output_dir, paste0(clean_phenotype, "_opentargets_genes.csv"))
  
  cat("🎯 Phenotype:        ", phenotype_name, "\n")
  cat("📁 Output directory: ", output_dir, "\n")
  cat("⏰ Start time:       ", format(Sys.time()), "\n\n")
  
  mondo_id <- get_mondo_id(phenotype_name)
  
  if (is.null(mondo_id)) {
    cat("❌ Could not resolve phenotype to a disease ID in Open Targets.\n")
    quit(status = 0)
  }
  
  genes_df <- get_genes(mondo_id)
  
  if (nrow(genes_df) > 0) {
    write.csv(genes_df, output_file, row.names = FALSE)
    
    genes_only <- data.frame(Gene = unique(genes_df$Gene))
    write.csv(genes_only, genes_only_file, row.names = FALSE)
    
    cat("\n✅ SUCCESS!\n")
    cat("📊 Total records:   ", nrow(genes_df), "\n")
    cat("🧬 Unique genes:    ", length(unique(genes_df$Gene)), "\n")
    cat("💾 Full results:    ", output_file, "\n")
    cat("🧬 Genes-only file: ", genes_only_file, "\n\n")
    
    unique_genes <- sort(unique(genes_df$Gene))
    cat("🧬 Genes found for", phenotype_name, ":\n")
    for (i in seq(1, length(unique_genes), by = 8)) {
      end_idx <- min(i + 7, length(unique_genes))
      cat("   ", paste(unique_genes[i:end_idx], collapse = ", "), "\n")
    }
    
  } else {
    cat("❌ No gene associations found for:", phenotype_name, "\n")
  }
  
  cat("\n⏰ End time: ", format(Sys.time()), "\n")
}

if (!interactive()) {
  main()
}