#!/usr/bin/env Rscript

# DisGeNET gene downloader
# Usage:
#   export DISGENET_API_KEY="your_new_api_key"
#   Rscript disgenet_fixed.R migraine
#
# Optional:
#   Rscript disgenet_fixed.R migraine ALL 0.1 1

#!/usr/bin/env Rscript

required_cran <- c("dplyr", "devtools")
# api_key <- "c1b23105-1df0-41c5-9788-d5833a6192ad"
api_key <- "c1b23105-1df0-41c5-9788-d5833a6192ad"
Sys.setenv(DISGENET_API_KEY = api_key)
#!/usr/bin/env Rscript

# Robust DisGeNET gene downloader
# Usage:
#   export DISGENET_API_KEY="your_api_key"
#   Rscript disgenet.R migraine
# Optional:
#   Rscript disgenet.R migraine ALL 0.1 1

if (!interactive() && is.null(getOption("repos"))) {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

required_cran <- c("dplyr", "devtools")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

if (!requireNamespace("disgenet2r", quietly = TRUE)) {
  devtools::install_gitlab("medbio/disgenet2r")
}
suppressPackageStartupMessages(library(disgenet2r))

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

extract_first_dataframe <- function(obj) {
  if (is.null(obj)) return(NULL)

  if (is.data.frame(obj)) return(obj)
  if (is.matrix(obj)) return(as.data.frame(obj, stringsAsFactors = FALSE))

  if (methods::is(obj, "S4")) {
    sn <- methods::slotNames(obj)
    preferred_slots <- c("qresult", "results", "result", "data", "content")
    for (nm in c(preferred_slots, sn)) {
      if (nm %in% sn) {
        val <- tryCatch(methods::slot(obj, nm), error = function(e) NULL)
        if (is.data.frame(val)) return(val)
        if (is.matrix(val)) return(as.data.frame(val, stringsAsFactors = FALSE))
        if (is.list(val)) {
          for (el in val) {
            if (is.data.frame(el)) return(el)
            if (is.matrix(el)) return(as.data.frame(el, stringsAsFactors = FALSE))
          }
        }
      }
    }
  }

  if (is.list(obj)) {
    preferred_names <- c("qresult", "results", "result", "data", "content")
    for (nm in c(preferred_names, names(obj))) {
      if (!is.null(nm) && nm %in% names(obj)) {
        val <- obj[[nm]]
        if (is.data.frame(val)) return(val)
        if (is.matrix(val)) return(as.data.frame(val, stringsAsFactors = FALSE))
      }
    }
    for (el in obj) {
      if (is.data.frame(el)) return(el)
      if (is.matrix(el)) return(as.data.frame(el, stringsAsFactors = FALSE))
    }
  }

  NULL
}

find_id_column <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)

  candidates <- c(
    "cui", "diseaseid", "disease_id", "diseaseId", "umls", "umls_id",
    "vocabulary_code", "code", "id"
  )

  cn_lower <- tolower(colnames(df))
  idx <- match(candidates, cn_lower)
  idx <- idx[!is.na(idx)]

  if (length(idx) > 0) {
    return(colnames(df)[idx[1]])
  }

  NULL
}

find_name_column <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)

  candidates <- c("name", "disease_name", "term", "description", "vocabularyname", "label")
  cn_lower <- tolower(colnames(df))
  idx <- match(candidates, cn_lower)
  idx <- idx[!is.na(idx)]

  if (length(idx) > 0) {
    return(colnames(df)[idx[1]])
  }

  NULL
}

pick_best_row <- function(df, phenotype) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)

  phenotype_lower <- tolower(trimws(phenotype))
  name_col <- find_name_column(df)

  if (!is.null(name_col)) {
    nm <- tolower(safe_trim(df[[name_col]]))

    exact_idx <- which(nm == phenotype_lower)
    if (length(exact_idx) > 0) return(df[exact_idx[1], , drop = FALSE])

    contains_idx <- which(grepl(phenotype_lower, nm, fixed = TRUE))
    if (length(contains_idx) > 0) return(df[contains_idx[1], , drop = FALSE])
  }

  df[1, , drop = FALSE]
}

build_disease_candidates <- function(raw_id, phenotype) {
  raw_id <- safe_trim(raw_id)
  phenotype <- safe_trim(phenotype)

  candidates <- c()

  if (nzchar(raw_id)) {
    candidates <- c(candidates, raw_id)

    if (grepl("^C[0-9]+$", raw_id)) {
      candidates <- c(candidates, paste0("UMLS_", raw_id))
    }

    if (grepl("^UMLS_C[0-9]+$", raw_id)) {
      candidates <- c(candidates, sub("^UMLS_", "", raw_id))
    }
  }

  if (nzchar(phenotype)) {
    candidates <- c(candidates, phenotype)
  }

  unique(candidates[nzchar(candidates)])
}

run_disease2gene <- function(disease_candidates, database, min_score, max_score) {
  for (d in disease_candidates) {
    cat("   Trying disease query:", d, "\n")

    res <- tryCatch(
      disease2gene(
        disease = d,
        database = database,
        score = c(min_score, max_score)
      ),
      error = function(e) {
        cat("   ⚠️ disease2gene failed for", d, ":", conditionMessage(e), "\n")
        NULL
      }
    )

    df <- extract_first_dataframe(res)

    if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
      cat("   ✅ Success with disease query:", d, "\n")
      return(list(query = d, data = df))
    }
  }

  NULL
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("DisGeNET Gene Downloader\n")
    cat("Usage:\n")
    cat("  Rscript disgenet.R <phenotype> [database] [min_score] [max_score]\n\n")
    cat("Examples:\n")
    cat("  Rscript disgenet.R migraine\n")
    cat("  Rscript disgenet.R migraine ALL 0.1 1\n")
    quit(status = 1)
  }

  phenotype <- args[1]
  database  <- if (length(args) > 1) toupper(args[2]) else "ALL"
  min_score <- if (length(args) > 2) as.numeric(args[3]) else 0.1
  max_score <- if (length(args) > 3) as.numeric(args[4]) else 1.0

  if (is.na(min_score) || is.na(max_score) || min_score < 0 || max_score > 1 || min_score > max_score) {
    stop("Score range must satisfy: 0 <= min_score <= max_score <= 1")
  }

  api_key <- Sys.getenv("DISGENET_API_KEY")
  if (!nzchar(safe_trim(api_key))) {
    stop("DISGENET_API_KEY is not set. Use: export DISGENET_API_KEY='your_api_key'")
  }

  Sys.setenv(DISGENET_API_KEY = api_key)

  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  clean_phenotype <- clean_filename(phenotype)
  output_file <- file.path(output_dir, paste0(clean_phenotype, "_disgenet.csv"))
  genes_only_file <- file.path(output_dir, paste0(clean_phenotype, "_disgenet_genes.csv"))

  cat("🎯 Phenotype:         ", phenotype, "\n")
  cat("🗂️ Database:          ", database, "\n")
  cat("📈 Score range:       ", min_score, "to", max_score, "\n")
  cat("📁 Output directory:  ", output_dir, "\n")
  cat("⏰ Start time:        ", format(Sys.time()), "\n\n")

  cat("🔍 Searching disease identifier for:", phenotype, "\n")

  search_obj <- tryCatch(
    get_umls_from_vocabulary(disease = phenotype, vocabulary = "NAME"),
    error = function(e) {
      cat("⚠️ Identifier search failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  search_df <- extract_first_dataframe(search_obj)

  raw_id <- ""
  matched_term <- ""

  if (!is.null(search_df) && is.data.frame(search_df) && nrow(search_df) > 0) {
    best_row <- pick_best_row(search_df, phenotype)

    id_col <- find_id_column(best_row)
    name_col <- find_name_column(best_row)

    if (!is.null(id_col)) raw_id <- safe_trim(best_row[[id_col]][1])
    if (!is.null(name_col)) matched_term <- safe_trim(best_row[[name_col]][1])
  }

  if (nzchar(raw_id)) {
    cat("✅ Selected disease ID:", raw_id, "\n")
    if (nzchar(matched_term)) {
      cat("✅ Matched term:      ", matched_term, "\n")
    }
  } else {
    cat("⚠️ No explicit disease ID extracted. Will try phenotype directly as fallback.\n")
  }

  disease_candidates <- build_disease_candidates(raw_id, phenotype)

  cat("🧬 Retrieving gene associations...\n")
  hit <- run_disease2gene(disease_candidates, database, min_score, max_score)

  if (is.null(hit)) {
    cat("❌ No results found for:", phenotype, "\n")
    cat("Tried disease queries:", paste(disease_candidates, collapse = ", "), "\n")
    cat("\n⏰ End time: ", format(Sys.time()), "\n")
    return(invisible(data.frame()))
  }

  gdas <- hit$data
  used_query <- hit$query

  write.csv(gdas, output_file, row.names = FALSE)

  genes_only <- data.frame()
  gene_col_candidates <- c("gene_symbol", "symbol", "geneSymbol")
  gene_col_candidates <- gene_col_candidates[gene_col_candidates %in% colnames(gdas)]

  if (length(gene_col_candidates) > 0) {
    gene_col <- gene_col_candidates[1]
    genes_only <- data.frame(
      Gene = sort(unique(safe_trim(gdas[[gene_col]][nzchar(safe_trim(gdas[[gene_col]]))]))),
      stringsAsFactors = FALSE
    )
    write.csv(genes_only, genes_only_file, row.names = FALSE)
  } else {
    cat("⚠️ No gene symbol column found, so genes-only file was not created.\n")
  }

  cat("\n✅ SUCCESS!\n")
  cat("🔎 Query used:         ", used_query, "\n")
  cat("📊 Total associations: ", nrow(gdas), "\n")
  if (nrow(genes_only) > 0) {
    cat("🧬 Unique genes:       ", nrow(genes_only), "\n")
  }
  cat("💾 Full results:       ", output_file, "\n")
  if (nrow(genes_only) > 0) {
    cat("🧬 Genes-only file:    ", genes_only_file, "\n")
  }

  show_cols <- intersect(c("gene_symbol", "disease_name", "score", "source"), colnames(gdas))
  if (length(show_cols) > 0) {
    cat("\n--- Top Gene-Disease Associations for", phenotype, "---\n")
    top_results <- gdas %>%
      dplyr::select(dplyr::all_of(show_cols))

    if ("score" %in% colnames(top_results)) {
      top_results <- top_results %>% dplyr::arrange(dplyr::desc(.data$score))
    }

    print(utils::head(top_results, 10))
  }

  cat("\n⏰ End time: ", format(Sys.time()), "\n")
  invisible(gdas)
}

if (!interactive()) {
  main()
}