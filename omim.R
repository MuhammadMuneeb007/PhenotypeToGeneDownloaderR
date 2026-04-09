#!/usr/bin/env Rscript
# OMIM gene downloader
# Primary:  OMIM API
# Fallback: OMIM web search scraping
#
# Usage:
#   export OMIM_API_KEY="your_key_here"
#   Rscript omim_combined.R migraine
#
# Or:
#   Rscript omim_combined.R migraine your_key_here

required_packages <- c("httr", "jsonlite", "rvest", "stringr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

omim_base_url <- "https://api.omim.org/api"

get_api_key <- function(args) {
  if (length(args) >= 2 && nzchar(args[2])) {
    return(args[2])
  }

  env_key <- Sys.getenv("OMIM_API_KEY", unset = "")
  if (nzchar(env_key)) {
    return(env_key)
  }

  ""
}

build_search_terms <- function(phenotype) {
  p <- tolower(trimws(phenotype))

  if (p == "migraine") {
    return(unique(c(
      "migraine",
      "\"migraine with aura\"",
      "\"migraine without aura\"",
      "\"familial hemiplegic migraine\"",
      "\"hemiplegic migraine\""
    )))
  }

  unique(c(
    phenotype,
    paste(phenotype, "disease"),
    paste(phenotype, "syndrome"),
    paste("familial", phenotype)
  ))
}

safe_get_text <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  as.character(x)
}

# ─────────────────────────────────────────────
# API path
# ─────────────────────────────────────────────
omim_get <- function(path, query = list(), api_key) {
  url <- paste0(omim_base_url, path)

  response <- tryCatch(
    GET(
      url,
      query = c(query, list(format = "json")),
      add_headers(
        ApiKey = api_key,
        `Accept-Encoding` = "gzip",
        `User-Agent` = "R-OMIM-API-GeneDownloader/1.0"
      ),
      timeout(60)
    ),
    error = function(e) {
      stop("Request failed: ", conditionMessage(e))
    }
  )

  status <- status_code(response)

  if (status != 200) {
    body_text <- tryCatch(
      content(response, "text", encoding = "UTF-8"),
      error = function(e) ""
    )
    stop(
      "OMIM API request failed. HTTP ", status,
      if (nzchar(body_text)) paste0("\n", body_text) else ""
    )
  }

  content(response, as = "parsed", type = "application/json", simplifyVector = FALSE)
}

extract_api_rows_from_entry <- function(entry_wrapper, search_term, phenotype) {
  entry <- entry_wrapper$entry %||% NULL
  if (is.null(entry)) return(list())

  mim_number <- entry$mimNumber %||% NA
  title <- entry$titles$preferredTitle %||% NA
  gene_map_list <- entry$geneMapList %||% list()

  if (length(gene_map_list) == 0) return(list())

  rows <- list()

  for (gm_wrap in gene_map_list) {
    gm <- gm_wrap$geneMap %||% NULL
    if (is.null(gm)) next

    approved_symbols <- safe_get_text(gm$approvedGeneSymbols)
    fallback_symbols <- safe_get_text(gm$geneSymbols)
    gene_name <- safe_get_text(gm$geneName)
    entrez_ids <- safe_get_text(gm$geneIDs)
    ensembl_ids <- safe_get_text(gm$ensemblIDs)

    raw_symbols <- if (nzchar(approved_symbols)) approved_symbols else fallback_symbols
    if (!nzchar(raw_symbols)) next

    symbols <- unique(trimws(unlist(strsplit(raw_symbols, ","))))
    symbols <- symbols[nzchar(symbols)]

    phenotype_map_list <- gm$phenotypeMapList %||% list()
    phenotype_texts <- character()
    inheritance_texts <- character()

    if (length(phenotype_map_list) > 0) {
      for (pm_wrap in phenotype_map_list) {
        pm <- pm_wrap$phenotypeMap %||% NULL
        if (is.null(pm)) next
        phenotype_texts <- c(phenotype_texts, safe_get_text(pm$phenotype))
        inheritance_texts <- c(inheritance_texts, safe_get_text(pm$phenotypeInheritance))
      }
    }

    phenotype_texts <- unique(trimws(phenotype_texts[nzchar(phenotype_texts)]))
    inheritance_texts <- unique(trimws(inheritance_texts[nzchar(inheritance_texts)]))

    phenotype_map_text <- paste(phenotype_texts, collapse = "; ")
    inheritance_text <- paste(inheritance_texts, collapse = "; ")

    keep_entry <- TRUE
    if (nzchar(phenotype_map_text)) {
      keep_entry <- str_detect(
        tolower(phenotype_map_text),
        fixed(tolower(phenotype))
      )
    }

    if (!keep_entry) next

    for (sym in symbols) {
      rows[[length(rows) + 1]] <- data.frame(
        Gene_Symbol = sym,
        Gene_Name = gene_name,
        MIM_Number = as.character(mim_number),
        Entry_Title = as.character(title),
        Search_Term = as.character(search_term),
        Phenotype_Map = phenotype_map_text,
        Inheritance = inheritance_text,
        Entrez_GeneIDs = entrez_ids,
        Ensembl_IDs = ensembl_ids,
        Source = "OMIM_API",
        stringsAsFactors = FALSE
      )
    }
  }

  rows
}

search_omim_api <- function(term, phenotype, api_key, max_results = 200) {
  cat("🔑 API search:", term, "\n")

  batch_size <- 20
  start <- 0
  collected <- list()

  repeat {
    resp <- omim_get(
      path = "/entry/search",
      query = list(
        search = term,
        include = "geneMap",
        sort = "score desc",
        start = start,
        limit = batch_size
      ),
      api_key = api_key
    )

    entry_list <- resp$omim$searchResponse$entryList %||% list()
    total_results <- resp$omim$searchResponse$totalResults %||% 0

    if (length(entry_list) == 0) break

    cat("   Retrieved", length(entry_list), "entries (start =", start, "of", total_results, ")\n")

    for (entry_wrapper in entry_list) {
      rows <- extract_api_rows_from_entry(entry_wrapper, term, phenotype)
      if (length(rows) > 0) {
        collected <- c(collected, rows)
      }
    }

    start <- start + length(entry_list)

    if (length(entry_list) < batch_size) break
    if (start >= total_results) break
    if (start >= max_results) break

    Sys.sleep(0.25)
  }

  if (length(collected) == 0) return(data.frame())

  out <- do.call(rbind, collected)
  out <- out[!duplicated(out[, c("Gene_Symbol", "MIM_Number", "Entry_Title")]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

run_api_mode <- function(phenotype, api_key) {
  if (!nzchar(api_key)) {
    cat("ℹ️ No OMIM API key found. Skipping API mode.\n")
    return(data.frame())
  }

  search_terms <- build_search_terms(phenotype)
  all_results <- list()

  for (term in search_terms) {
    term_results <- tryCatch(
      search_omim_api(term, phenotype, api_key, max_results = 200),
      error = function(e) {
        cat("❌ API error for", shQuote(term), ":", conditionMessage(e), "\n")
        data.frame()
      }
    )

    if (nrow(term_results) > 0) {
      all_results[[length(all_results) + 1]] <- term_results
    }
  }

  if (length(all_results) == 0) return(data.frame())

  final_df <- do.call(rbind, all_results)
  final_df <- final_df[!duplicated(final_df[, c("Gene_Symbol", "MIM_Number", "Entry_Title")]), , drop = FALSE]
  final_df <- final_df[order(final_df$Gene_Symbol), , drop = FALSE]
  rownames(final_df) <- NULL
  final_df
}

# ─────────────────────────────────────────────
# Scraping fallback
# ─────────────────────────────────────────────
extract_candidate_genes_from_text <- function(text_content) {
  if (!nzchar(text_content)) return(character())

  gene_patterns <- c(
    "\\b[A-Z][A-Z0-9]{2,14}\\b",
    "\\b[A-Z][A-Z0-9-]{2,15}\\b",
    "\\b[A-Z]{2,}[0-9]+[A-Z]*\\b"
  )

  exclude_words <- c(
    "OMIM", "SEARCH", "RESULTS", "PAGE", "HOME", "GENE", "LOCUS",
    "PHENOTYPE", "SYNDROME", "DISEASE", "DISORDER", "MUTATION",
    "THE", "AND", "FOR", "WITH", "THAT", "THIS", "FROM", "WERE",
    "HAVE", "MORE", "TIME", "VERY", "CAN", "HAD", "HER", "WAS",
    "ONE", "OUR", "OUT", "DAY", "GET", "HAS", "HIM", "HOW", "ITS",
    "MAY", "NEW", "NOW", "OLD", "SEE", "TWO", "WHO", "DID", "HIS",
    "LET", "PUT", "SAY", "SHE", "TOO", "USE", "ALL", "ANY", "ARE",
    "BUT", "NOT", "YOU", "WHAT", "WHEN", "WHERE", "WHY", "ABOUT",
    "AFTER", "CLINICAL", "MOLECULAR", "GENETICS", "TYPE", "NULL",
    "TRUE", "FALSE", "LIST", "NEXT", "PREV", "LAST", "FIRST",
    "SIGN", "ALSO", "BOTH", "SOME", "SUCH", "THAN", "THEN", "THEY",
    "THUS", "WELL", "WILL", "BEEN", "EACH", "EVEN", "INTO",
    "ONLY", "SAME", "SHOW", "THEM", "THESE", "THEIR", "THERE"
  )

  matches <- character()
  for (pattern in gene_patterns) {
    matches <- c(matches, str_extract_all(text_content, pattern)[[1]])
  }

  matches <- unique(trimws(matches))
  matches <- matches[
    nzchar(matches) &
    !matches %in% exclude_words &
    nchar(matches) >= 3 &
    nchar(matches) <= 15 &
    !str_detect(matches, "^[0-9]+$") &
    str_detect(matches, "^[A-Z]") &
    str_detect(matches, "[A-Z]")
  ]

  sort(unique(matches))
}

search_omim_scrape <- function(term, phenotype) {
  cat("🌐 Scraping search:", term, "\n")

  url <- paste0(
    "https://www.omim.org/search?search=",
    URLencode(term, reserved = TRUE),
    "&sort=score+desc&limit=100"
  )

  response <- tryCatch(
    GET(
      url,
      add_headers(
        `User-Agent` = "Mozilla/5.0 (compatible; R-OMIMScraper/1.0; +research)"
      ),
      timeout(30)
    ),
    error = function(e) {
      cat("❌ Scraping request failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(response)) return(data.frame())

  if (status_code(response) != 200) {
    cat("⚠️ Scraping HTTP status:", status_code(response), "\n")
    return(data.frame())
  }

  page_content <- content(response, "text", encoding = "UTF-8")
  doc <- read_html(page_content)
  text_content <- html_text(doc)

  genes <- extract_candidate_genes_from_text(text_content)

  if (length(genes) == 0) return(data.frame())

  data.frame(
    Gene_Symbol = genes,
    Gene_Name = "",
    MIM_Number = "",
    Entry_Title = "",
    Search_Term = term,
    Phenotype_Map = phenotype,
    Inheritance = "",
    Entrez_GeneIDs = "",
    Ensembl_IDs = "",
    Source = "OMIM_SCRAPE_FALLBACK",
    stringsAsFactors = FALSE
  )
}

run_scrape_mode <- function(phenotype) {
  search_terms <- build_search_terms(phenotype)
  all_results <- list()

  for (term in search_terms) {
    res <- tryCatch(
      search_omim_scrape(term, phenotype),
      error = function(e) {
        cat("❌ Scraping error for", shQuote(term), ":", conditionMessage(e), "\n")
        data.frame()
      }
    )

    if (nrow(res) > 0) {
      all_results[[length(all_results) + 1]] <- res
    }

    Sys.sleep(1)
  }

  if (length(all_results) == 0) return(data.frame())

  final_df <- do.call(rbind, all_results)
  final_df <- final_df[!duplicated(final_df[, c("Gene_Symbol", "Search_Term", "Source")]), , drop = FALSE]
  final_df <- final_df[order(final_df$Gene_Symbol), , drop = FALSE]
  rownames(final_df) <- NULL
  final_df
}

# ─────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────
save_results <- function(df, phenotype) {
  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  clean_phenotype <- gsub("[^[:alnum:]_ -]", "", phenotype)
  clean_phenotype <- gsub("\\s+", "_", trimws(clean_phenotype))

  full_file <- file.path(output_dir, paste0(clean_phenotype, "_omim.csv"))
  genes_file <- file.path(output_dir, paste0(clean_phenotype, "_omim_genes.csv"))

  write.csv(df, full_file, row.names = FALSE)

  genes_only <- data.frame(
    Gene = sort(unique(df$Gene_Symbol)),
    stringsAsFactors = FALSE
  )
  write.csv(genes_only, genes_file, row.names = FALSE)

  list(full_file = full_file, genes_file = genes_file, genes_only = genes_only)
}

# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("OMIM Gene Downloader\n")
    cat("Usage:\n")
    cat("  export OMIM_API_KEY='your_key_here'\n")
    cat("  Rscript omim_combined.R <phenotype>\n\n")
    cat("Or:\n")
    cat("  Rscript omim_combined.R <phenotype> <api_key>\n")
    quit(status = 1)
  }

  phenotype <- args[1]
  api_key <- get_api_key(args)

  cat("🎯 Phenotype:        ", phenotype, "\n")
  cat("📁 Output directory:  AllPackagesGenes\n")
  cat("🔑 API key:          ", ifelse(nzchar(api_key), "provided", "not provided"), "\n")
  cat("⏰ Start time:       ", format(Sys.time()), "\n\n")

  final_df <- data.frame()

  # Try API first
  final_df <- run_api_mode(phenotype, api_key)

  # Fallback if needed
  if (nrow(final_df) == 0) {
    cat("\n⚠️ API mode unavailable or returned no results. Switching to OMIM scraping fallback.\n\n")
    final_df <- run_scrape_mode(phenotype)
  }

  if (nrow(final_df) == 0) {
    cat("❌ No genes found for:", phenotype, "\n")
    cat("\n⏰ End time: ", format(Sys.time()), "\n")
    quit(status = 0)
  }

  saved <- save_results(final_df, phenotype)

  cat("\n✅ SUCCESS!\n")
  cat("📊 Total rows:      ", nrow(final_df), "\n")
  cat("🧬 Unique genes:     ", nrow(saved$genes_only), "\n")
  cat("🔧 Method(s):        ", paste(unique(final_df$Source), collapse = ", "), "\n")
  cat("💾 Full results:     ", saved$full_file, "\n")
  cat("🧬 Genes-only file:  ", saved$genes_file, "\n\n")

  genes <- saved$genes_only$Gene
  cat("🧬 Genes found for", phenotype, ":\n")
  for (i in seq(1, length(genes), by = 8)) {
    end_idx <- min(i + 7, length(genes))
    cat("   ", paste(genes[i:end_idx], collapse = ", "), "\n")
  }

  cat("\n⏰ End time: ", format(Sys.time()), "\n")
}

if (!interactive()) {
  main()
}