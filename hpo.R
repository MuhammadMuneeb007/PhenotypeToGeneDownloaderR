#!/usr/bin/env Rscript
# Robust HPO gene downloader

if (!interactive() && is.null(getOption("repos"))) {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("⚠️  Missing package:", pkg, "- attempting to install\n")
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org/")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

required_packages <- c("ontologyIndex", "httr", "utils")
for (pkg in required_packages) install_if_missing(pkg)

safe_read_tsv <- function(path, comment = "#") {
  tryCatch(
    read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = comment,
      stringsAsFactors = FALSE,
      fill = TRUE,
      check.names = FALSE
    ),
    error = function(e) NULL
  )
}

download_file_if_needed <- function(url, dest, max_age_days = 7) {
  if (file.exists(dest)) {
    age <- as.numeric(difftime(Sys.time(), file.mtime(dest), units = "days"))
    if (!is.na(age) && age < max_age_days) {
      cat("   ✅ Using cached file:", dest, "(", round(age, 1), "days old)\n")
      return(dest)
    }
  }

  cat("   📥 Downloading:", url, "\n")
  ok <- tryCatch({
    r <- httr::GET(url, httr::timeout(300))
    if (httr::status_code(r) != 200) return(FALSE)
    bin <- httr::content(r, "raw")
    writeBin(bin, dest)
    TRUE
  }, error = function(e) {
    cat("   ❌ Download error:", conditionMessage(e), "\n")
    FALSE
  })

  if (!ok || !file.exists(dest) || file.info(dest)$size < 100) {
    cat("   ❌ Failed or invalid file:", dest, "\n")
    return(NULL)
  }

  cat("   ✅ Downloaded:", dest, "\n")
  dest
}

download_all_hpo_files <- function() {
  cat("📁 Checking HPO data files...\n")

  files <- list(
    genes_to_phenotype = list(
      url = "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt",
      cache = "hpo_genes_to_phenotype.txt"
    ),
    phenotype_hpoa = list(
      url = "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa",
      cache = "hpo_phenotype.hpoa"
    ),
    hpo_ontology = list(
      url = "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo",
      cache = "hpo_ontology.obo"
    )
  )

  out <- list()
  for (nm in names(files)) {
    out[[nm]] <- download_file_if_needed(files[[nm]]$url, files[[nm]]$cache, max_age_days = 7)
    Sys.sleep(0.5)
  }

  out
}

HPO_FILES <- download_all_hpo_files()

download_hpo_ontology <- function() {
  cat("🔬 Loading HPO ontology...\n")

  cache_rds <- "hpo_ontology.rds"
  if (file.exists(cache_rds)) {
    age <- as.numeric(difftime(Sys.time(), file.mtime(cache_rds), units = "days"))
    if (!is.na(age) && age < 7) {
      cat("   ✅ Loaded ontology from cache\n")
      return(readRDS(cache_rds))
    }
  }

  if (is.null(HPO_FILES$hpo_ontology) || !file.exists(HPO_FILES$hpo_ontology)) {
    cat("   ❌ Ontology file not available\n")
    return(NULL)
  }

  hpo <- tryCatch({
    ontologyIndex::get_ontology(HPO_FILES$hpo_ontology, extract_tags = "minimal")
  }, error = function(e) {
    cat("   ❌ Ontology parse error:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(hpo)) {
    saveRDS(hpo, cache_rds)
    cat("   ✅ Parsed and cached ontology\n")
  }

  hpo
}

load_gene_hpo_annotations <- function() {
  cat("📥 Loading gene-HPO annotations...\n")

  cache_rds <- "gene_hpo_annotations.rds"
  if (file.exists(cache_rds)) {
    age <- as.numeric(difftime(Sys.time(), file.mtime(cache_rds), units = "days"))
    if (!is.na(age) && age < 1) {
      x <- readRDS(cache_rds)
      cat("   ✅ Loaded", nrow(x), "annotations from cache\n")
      return(x)
    }
  }

  if (is.null(HPO_FILES$genes_to_phenotype) || !file.exists(HPO_FILES$genes_to_phenotype)) {
    cat("   ❌ genes_to_phenotype file not available\n")
    return(data.frame())
  }

  df <- safe_read_tsv(HPO_FILES$genes_to_phenotype)
  if (is.null(df) || nrow(df) == 0) {
    cat("   ❌ Could not parse genes_to_phenotype file\n")
    return(data.frame())
  }

  required <- c("ncbi_gene_id", "gene_symbol", "hpo_id", "hpo_name")
  if (!all(required %in% colnames(df))) {
    cat("   ❌ Expected columns not found. Columns were:\n")
    cat("   ", paste(colnames(df), collapse = ", "), "\n")
    return(data.frame())
  }

  df$ncbi_gene_id <- trimws(as.character(df$ncbi_gene_id))
  df$gene_symbol  <- trimws(as.character(df$gene_symbol))
  df$hpo_id       <- trimws(as.character(df$hpo_id))
  df$hpo_name     <- trimws(as.character(df$hpo_name))

  keep <- !is.na(df$gene_symbol) &
    nzchar(df$gene_symbol) &
    !grepl("^[0-9]+$", df$gene_symbol) &
    !is.na(df$hpo_id) &
    grepl("^HP:[0-9]{7}$", df$hpo_id) &
    !is.na(df$hpo_name) &
    nzchar(df$hpo_name)

  df <- df[keep, , drop = FALSE]

  out <- data.frame(
    Gene_ID = df$ncbi_gene_id,
    Gene = df$gene_symbol,
    HPO_ID = df$hpo_id,
    HPO_Term = df$hpo_name,
    stringsAsFactors = FALSE
  )

  out <- out[!duplicated(paste(out$Gene, out$HPO_ID)), , drop = FALSE]
  saveRDS(out, cache_rds)

  cat("   ✅ Loaded", nrow(out), "unique gene-HPO annotations\n")
  cat("   🧬 Unique genes:", length(unique(out$Gene)), "\n")
  cat("   🏥 Unique HPO terms:", length(unique(out$HPO_ID)), "\n")

  out
}

load_disease_annotations <- function() {
  cat("📥 Loading disease annotations...\n")

  cache_rds <- "hpo_disease_annotations.rds"
  if (file.exists(cache_rds)) {
    age <- as.numeric(difftime(Sys.time(), file.mtime(cache_rds), units = "days"))
    if (!is.na(age) && age < 7) {
      x <- readRDS(cache_rds)
      cat("   ✅ Loaded", nrow(x), "disease annotations from cache\n")
      return(x)
    }
  }

  if (is.null(HPO_FILES$phenotype_hpoa) || !file.exists(HPO_FILES$phenotype_hpoa)) {
    cat("   ❌ phenotype.hpoa file not available\n")
    return(data.frame())
  }

  lines <- readLines(HPO_FILES$phenotype_hpoa, warn = FALSE)
  lines <- lines[!grepl("^#", lines) & nzchar(trimws(lines))]

  parsed <- lapply(lines, function(x) strsplit(x, "\t", fixed = TRUE)[[1]])
  parsed <- parsed[vapply(parsed, length, integer(1)) >= 4]

  if (length(parsed) == 0) return(data.frame())

  df <- data.frame(
    Disease_ID = vapply(parsed, function(x) trimws(x[1]), character(1)),
    Disease_Name = vapply(parsed, function(x) trimws(x[2]), character(1)),
    Qualifier = vapply(parsed, function(x) if (length(x) >= 3) trimws(x[3]) else "", character(1)),
    HPO_ID = vapply(parsed, function(x) trimws(x[4]), character(1)),
    stringsAsFactors = FALSE
  )

  df <- df[
    nzchar(df$Disease_ID) &
    nzchar(df$Disease_Name) &
    grepl("^HP:[0-9]{7}$", df$HPO_ID),
    ,
    drop = FALSE
  ]

  df <- df[!duplicated(paste(df$Disease_ID, df$HPO_ID)), , drop = FALSE]
  saveRDS(df, cache_rds)

  cat("   ✅ Loaded", nrow(df), "unique disease annotations\n")
  df
}

medical_variations <- list(
  migraine = c("migraine", "cephalgia", "hemicrania"),
  headache = c("headache", "cephalgia", "cephalalgia"),
  seizure = c("seizure", "epilepsy", "convulsion"),
  autism = c("autism", "autistic", "asd"),
  intellectual_disability = c("intellectual disability", "developmental delay", "global developmental delay")
)

create_search_patterns <- function(phenotype) {
  x <- tolower(trimws(phenotype))
  words <- unlist(strsplit(x, "\\s+"))
  words <- words[nchar(words) > 2]

  patterns <- unique(c(
    x,
    paste0("\\b", x, "\\b")
  ))

  if (length(words) > 1) {
    patterns <- c(patterns, paste(words, collapse = " "))
  }

  for (nm in names(medical_variations)) {
    if (grepl(nm, x, fixed = TRUE) || grepl(x, nm, fixed = TRUE)) {
      patterns <- c(patterns, medical_variations[[nm]])
      break
    }
  }

  unique(patterns)
}

calculate_relevance_score <- function(source) {
  scores <- c(
    HPO_Term_Exact = 13,
    HPO_Term_Contains = 12,
    HPO_Disease_Exact = 11,
    HPO_Disease_Contains = 9,
    HPO_Term_Pattern = 8,
    HPO_Disease_Pattern = 7,
    KnownAssociations = 3
  )
  if (source %in% names(scores)) scores[[source]] else 1
}

get_known_phenotype_genes <- function(phenotype) {
  known_associations <- list(
    migraine = c("CACNA1A", "ATP1A2", "SCN1A", "KCNK18", "PRRT2"),
    headache = c("CACNA1A", "ATP1A2", "SCN1A", "KCNK18", "PRRT2"),
    epilepsy = c("SCN1A", "SCN2A", "SCN8A", "KCNQ2", "KCNQ3", "CDKL5", "ARX", "STXBP1"),
    seizure = c("SCN1A", "SCN2A", "SCN8A", "KCNQ2", "KCNQ3"),
    autism = c("SHANK3", "NRXN1", "NLGN3", "NLGN4", "MECP2", "FMR1", "TSC1", "TSC2"),
    diabetes = c("INS", "GCGR", "HNF1A", "HNF4A", "GCK", "PDX1"),
    cancer = c("TP53", "BRCA1", "BRCA2", "APC", "MLH1", "MSH2", "PTEN")
  )

  p <- tolower(trimws(phenotype))
  out <- c()
  for (k in names(known_associations)) {
    if (grepl(k, p, fixed = TRUE) || grepl(p, k, fixed = TRUE)) {
      out <- c(out, known_associations[[k]])
    }
  }
  unique(out)
}

download_hpo_alternative <- function(phenotype) {
  cat("🔄 Using fallback approach for:", phenotype, "\n")
  known_genes <- get_known_phenotype_genes(phenotype)

  if (length(known_genes) == 0) return(data.frame())

  data.frame(
    Gene = known_genes,
    HPO_ID = "Known_Association",
    HPO_Term = phenotype,
    Gene_ID = "",
    Source = "KnownAssociations",
    Match_Pattern = phenotype,
    stringsAsFactors = FALSE
  )
}

download_hpo_genes <- function(phenotype) {
  cat("🔍 Searching for genes associated with:", phenotype, "\n")

  annotations <- load_gene_hpo_annotations()
  disease_annotations <- load_disease_annotations()

  if (nrow(annotations) == 0) {
    cat("   ⚠️  No direct HPO annotations available\n")
    return(download_hpo_alternative(phenotype))
  }

  patterns <- create_search_patterns(phenotype)
  phenotype_lower <- tolower(trimws(phenotype))

  annotations$HPO_Term_Lower <- tolower(annotations$HPO_Term)
  all_hits <- list()

  cat("   🔎 Search patterns:", paste(head(patterns, 5), collapse = ", "), "\n")

  for (pat in patterns) {
    pat_lower <- tolower(pat)
    idx <- grepl(pat_lower, annotations$HPO_Term_Lower, fixed = TRUE)

    if (!any(idx)) next

    tmp <- annotations[idx, c("Gene", "HPO_ID", "HPO_Term", "Gene_ID"), drop = FALSE]

    exact <- tmp$HPO_Term_Lower <- tolower(tmp$HPO_Term)
    src <- ifelse(
      exact == phenotype_lower,
      "HPO_Term_Exact",
      ifelse(grepl(phenotype_lower, exact, fixed = TRUE), "HPO_Term_Contains", "HPO_Term_Pattern")
    )

    tmp$Source <- src
    tmp$Match_Pattern <- pat
    all_hits[[length(all_hits) + 1]] <- tmp
  }

  res <- if (length(all_hits) > 0) do.call(rbind, all_hits) else data.frame()

  if (nrow(disease_annotations) > 0 && nrow(res) < 50) {
    disease_annotations$Disease_Name_Lower <- tolower(disease_annotations$Disease_Name)
    didx <- grepl(phenotype_lower, disease_annotations$Disease_Name_Lower, fixed = TRUE)

    if (any(didx)) {
      matched_hpo <- unique(disease_annotations$HPO_ID[didx])
      extra <- annotations[annotations$HPO_ID %in% matched_hpo, c("Gene", "HPO_ID", "HPO_Term", "Gene_ID"), drop = FALSE]
      if (nrow(extra) > 0) {
        extra$Source <- "HPO_Disease_Pattern"
        extra$Match_Pattern <- "disease_match"
        res <- rbind(res, extra)
      }
    }
  }

  if (nrow(res) == 0) {
    cat("   ⚠️  No matches found in HPO annotations\n")
    return(download_hpo_alternative(phenotype))
  }

  res <- res[!duplicated(paste(res$Gene, res$HPO_ID)), , drop = FALSE]
  res$Relevance_Score <- vapply(res$Source, calculate_relevance_score, numeric(1))
  res <- res[order(-res$Relevance_Score, res$Gene), , drop = FALSE]
  res$Relevance_Score <- NULL

  cat("   ✅ Total associations:", nrow(res), "\n")
  cat("   🧬 Unique genes:", length(unique(res$Gene)), "\n")

  res
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("Usage: Rscript hpo.R <phenotype> [output_file]\n")
    cat("Example: Rscript hpo.R migraine\n")
    quit(status = 1)
  }

  phenotype <- args[1]

  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  clean_phenotype <- gsub("[^A-Za-z0-9_-]", "_", phenotype)
  clean_phenotype <- gsub("_{2,}", "_", clean_phenotype)
  clean_phenotype <- gsub("^_|_$", "", clean_phenotype)

  output_file <- if (length(args) > 1) args[2] else file.path(output_dir, paste0(clean_phenotype, "_hpo.csv"))
  genes_only_file <- file.path(output_dir, paste0(clean_phenotype, "_hpo_genes.csv"))

  cat("🚀 Starting HPO gene search for:", phenotype, "\n")
  cat("📁 Output directory:", output_dir, "\n")

  results <- download_hpo_genes(phenotype)

  if (nrow(results) > 0) {
    write.csv(results, output_file, row.names = FALSE)
    write.csv(data.frame(Gene = unique(results$Gene), stringsAsFactors = FALSE), genes_only_file, row.names = FALSE)

    cat("✅ Results saved to:", output_file, "\n")
    cat("🧬 Genes-only file saved to:", genes_only_file, "\n")
    cat("📊 Found", nrow(results), "gene associations\n")
    cat("🧬 Unique genes:", length(unique(results$Gene)), "\n")
  } else {
    cat("❌ No genes found for phenotype:", phenotype, "\n")
    cat("💡 Try different search terms or check HPO data availability\n")
  }
}

if (!interactive()) {
  main()
}