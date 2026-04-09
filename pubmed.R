#!/usr/bin/env Rscript
# Reliable PubMed -> Gene downloader using PubMed E-utilities + PubTator3
# Usage:
#   Rscript pubmed.R <phenotype> [output_file]

required_packages <- c("httr", "jsonlite", "xml2")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_trim <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  trimws(as.character(x)[1])
}

chunk_vector <- function(x, chunk_size = 100) {
  if (length(x) == 0) return(list())
  split(x, ceiling(seq_along(x) / chunk_size))
}

get_ncbi_params <- function() {
  email <- Sys.getenv("NCBI_EMAIL", unset = "")
  api_key <- Sys.getenv("NCBI_API_KEY", unset = "")

  params <- list(tool = "phenotype_gene_downloader_r")
  if (nzchar(email)) params$email <- email
  if (nzchar(api_key)) params$api_key <- api_key
  params
}

search_pubmed_pmids <- function(phenotype, retmax = 300) {
  cat("🔍 Searching PubMed for phenotype:", phenotype, "\n")

  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

  search_terms <- unique(c(
    paste0("\"", phenotype, "\"[Title/Abstract] AND (gene OR genes OR genetic OR genomics OR mutation OR variant OR polymorphism)"),
    paste0("\"", phenotype, "\"[MeSH Terms] AND (gene OR genes OR genetic OR mutation OR variant)"),
    paste0("\"", phenotype, "\"[Title/Abstract] AND (GWAS OR genome-wide association OR transcriptome OR sequencing)")
  ))

  all_pmids <- character()

  for (term in search_terms) {
    cat("   Query:", term, "\n")

    query_params <- c(
      list(
        db = "pubmed",
        term = term,
        retmax = retmax,
        retmode = "xml",
        sort = "relevance"
      ),
      get_ncbi_params()
    )

    resp <- tryCatch(
      httr::GET(base_url, query = query_params, httr::timeout(60)),
      error = function(e) NULL
    )

    if (is.null(resp)) {
      cat("   ⚠️ Request failed for this query\n")
      next
    }

    if (httr::status_code(resp) != 200) {
      cat("   ⚠️ PubMed returned status", httr::status_code(resp), "\n")
      next
    }

    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    doc <- tryCatch(xml2::read_xml(txt), error = function(e) NULL)

    if (is.null(doc)) {
      cat("   ⚠️ Could not parse PubMed XML\n")
      next
    }

    pmids <- xml2::xml_text(xml2::xml_find_all(doc, "//IdList/Id"))
    pmids <- unique(pmids[nzchar(pmids)])

    cat("   Found", length(pmids), "PMIDs\n")
    all_pmids <- c(all_pmids, pmids)

    Sys.sleep(0.34)
  }

  all_pmids <- unique(all_pmids)
  cat("   Total unique PMIDs:", length(all_pmids), "\n")
  all_pmids
}

fetch_pubmed_metadata <- function(pmids) {
  if (length(pmids) == 0) return(data.frame())

  cat("📰 Fetching PubMed metadata for", length(pmids), "PMIDs...\n")

  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  chunks <- chunk_vector(pmids, 100)
  rows <- list()

  for (idx in seq_along(chunks)) {
    batch <- chunks[[idx]]

    query_params <- c(
      list(
        db = "pubmed",
        id = paste(batch, collapse = ","),
        retmode = "xml"
      ),
      get_ncbi_params()
    )

    resp <- tryCatch(
      httr::GET(base_url, query = query_params, httr::timeout(60)),
      error = function(e) NULL
    )

    if (is.null(resp) || httr::status_code(resp) != 200) {
      cat("   ⚠️ Metadata fetch failed for batch", idx, "\n")
      Sys.sleep(0.34)
      next
    }

    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    doc <- tryCatch(xml2::read_xml(txt), error = function(e) NULL)

    if (is.null(doc)) {
      cat("   ⚠️ Could not parse metadata XML for batch", idx, "\n")
      Sys.sleep(0.34)
      next
    }

    articles <- xml2::xml_find_all(doc, "//PubmedArticle")

    for (article in articles) {
      pmid <- safe_trim(xml2::xml_text(xml2::xml_find_first(article, ".//PMID")))
      title <- safe_trim(xml2::xml_text(xml2::xml_find_first(article, ".//ArticleTitle")))
      journal <- safe_trim(xml2::xml_text(xml2::xml_find_first(article, ".//Journal/Title")))
      year <- safe_trim(xml2::xml_text(xml2::xml_find_first(article, ".//PubDate/Year")))

      if (!nzchar(year)) {
        medline_date <- safe_trim(xml2::xml_text(xml2::xml_find_first(article, ".//PubDate/MedlineDate")))
        year <- sub("^([0-9]{4}).*", "\\1", medline_date)
        if (!grepl("^[0-9]{4}$", year)) year <- NA_character_
      }

      rows[[length(rows) + 1]] <- data.frame(
        PMID = pmid,
        Title = title,
        Journal = journal,
        Year = year,
        stringsAsFactors = FALSE
      )
    }

    Sys.sleep(0.34)
  }

  if (length(rows) == 0) return(data.frame())

  meta <- do.call(rbind, rows)
  meta <- meta[!duplicated(meta$PMID), , drop = FALSE]
  meta
}

normalize_pubtator_documents <- function(parsed) {
  if (is.null(parsed)) return(list())

  # New observed wrapper: top-level PubTator3
  if (is.list(parsed) && !is.null(parsed$PubTator3)) {
    parsed <- parsed$PubTator3
  }

  # Case 1: BioC collection with documents
  if (is.list(parsed) && !is.null(parsed$documents)) {
    return(parsed$documents)
  }

  # Case 2: nested collection object
  if (is.list(parsed) && !is.null(parsed$collection) && !is.null(parsed$collection$documents)) {
    return(parsed$collection$documents)
  }

  # Case 3: a single document object
  if (is.list(parsed) && (!is.null(parsed$id) || !is.null(parsed$passages) || !is.null(parsed$annotations))) {
    return(list(parsed))
  }

  # Case 4: already a list of document objects
  if (is.list(parsed) && length(parsed) > 0) {
    looks_like_doc <- function(x) {
      is.list(x) && (!is.null(x$id) || !is.null(x$passages) || !is.null(x$annotations))
    }
    if (all(vapply(parsed, looks_like_doc, logical(1)))) {
      return(parsed)
    }
  }

  list()
}

extract_annotations_from_doc <- function(doc) {
  out <- list()

  pmid <- safe_trim(doc$id %||% doc$pmid %||% NA_character_)

  doc_annotations <- doc$annotations %||% list()
  if (length(doc_annotations) > 0) {
    out[[length(out) + 1]] <- list(
      pmid = pmid,
      passage_text = "",
      annotations = doc_annotations
    )
  }

  passages <- doc$passages %||% list()
  if (length(passages) > 0) {
    for (passage in passages) {
      out[[length(out) + 1]] <- list(
        pmid = pmid,
        passage_text = safe_trim(passage$text %||% ""),
        annotations = passage$annotations %||% list()
      )
    }
  }

  out
}

extract_gene_rows_from_context <- function(context) {
  rows <- list()
  pmid <- context$pmid %||% ""
  passage_text <- context$passage_text %||% ""
  annotations <- context$annotations %||% list()

  if (length(annotations) == 0) return(rows)

  for (ann in annotations) {
    infons <- ann$infons %||% list()

    ann_type <- tolower(safe_trim(
      infons$type %||%
      infons$concept %||%
      infons$category %||%
      ""
    ))

    is_gene <- ann_type %in% c("gene", "genes", "geneproduct", "gene_product")
    if (!is_gene && !grepl("gene", ann_type, fixed = TRUE)) next

    gene_text <- safe_trim(ann$text %||% ann$mention %||% "")
    if (!nzchar(gene_text)) next

    gene_id <- safe_trim(
      infons$identifier %||%
      infons$Identifier %||%
      infons$geneid %||%
      infons$GeneID %||%
      NA_character_
    )

    gene_ids <- unlist(strsplit(gene_id, "[,;|]"))
    gene_ids <- unique(trimws(gene_ids))
    gene_ids <- gene_ids[nzchar(gene_ids)]
    if (length(gene_ids) == 0) gene_ids <- NA_character_

    rows[[length(rows) + 1]] <- data.frame(
      PMID = as.character(pmid),
      Gene = as.character(gene_text),
      Entrez_ID = paste(gene_ids, collapse = ";"),
      Evidence_Text = substr(as.character(passage_text), 1, 500),
      Source = "PubTator3",
      stringsAsFactors = FALSE
    )
  }

  rows
}

extract_pubtator_genes <- function(pmids) {
  if (length(pmids) == 0) return(data.frame())

  cat("🧬 Fetching PubTator gene annotations for", length(pmids), "PMIDs...\n")

  base_url <- "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson"
  chunks <- chunk_vector(pmids, 100)
  all_rows <- list()

  for (idx in seq_along(chunks)) {
    batch <- chunks[[idx]]
    cat("   Batch", idx, "of", length(chunks), "- PMIDs:", length(batch), "\n")

    resp <- tryCatch(
      httr::GET(
        base_url,
        query = list(
          pmids = paste(batch, collapse = ","),
          concepts = "gene"
        ),
        httr::timeout(120)
      ),
      error = function(e) NULL
    )

    if (is.null(resp)) {
      cat("   ⚠️ PubTator request failed for batch", idx, "\n")
      Sys.sleep(0.5)
      next
    }

    if (httr::status_code(resp) != 200) {
      cat("   ⚠️ PubTator returned status", httr::status_code(resp), "for batch", idx, "\n")
      Sys.sleep(0.5)
      next
    }

    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    if (!nzchar(txt)) {
      cat("   ⚠️ Empty PubTator response for batch", idx, "\n")
      Sys.sleep(0.5)
      next
    }

    parsed <- tryCatch(
      jsonlite::fromJSON(txt, simplifyVector = FALSE),
      error = function(e) NULL
    )

    if (is.null(parsed)) {
      cat("   ⚠️ Could not parse PubTator JSON for batch", idx, "\n")
      cat("   Debug preview:", substr(txt, 1, 200), "\n")
      Sys.sleep(0.5)
      next
    }

    docs <- normalize_pubtator_documents(parsed)

    if (length(docs) == 0) {
      cat("   ⚠️ No documents recognized in PubTator response for batch", idx, "\n")
      if (idx == 1) {
        cat("   Debug top-level names:", paste(names(parsed), collapse = ", "), "\n")
      }
      Sys.sleep(0.5)
      next
    }

    batch_rows_before <- length(all_rows)

    for (doc in docs) {
      contexts <- extract_annotations_from_doc(doc)
      if (length(contexts) == 0) next

      for (context in contexts) {
        gene_rows <- extract_gene_rows_from_context(context)
        if (length(gene_rows) > 0) {
          all_rows <- c(all_rows, gene_rows)
        }
      }
    }

    batch_rows_after <- length(all_rows)
    cat("   Extracted", batch_rows_after - batch_rows_before, "gene annotation rows in this batch\n")

    Sys.sleep(0.5)
  }

  if (length(all_rows) == 0) return(data.frame())

  out <- do.call(rbind, all_rows)
  out <- out[nzchar(out$Gene), , drop = FALSE]
  out$Gene <- trimws(out$Gene)
  out <- out[nchar(out$Gene) <= 50, , drop = FALSE]
  out <- out[!duplicated(out[, c("PMID", "Gene", "Entrez_ID")]), , drop = FALSE]

  out
}

build_gene_summary <- function(annotation_df, phenotype) {
  if (nrow(annotation_df) == 0) return(data.frame())

  annotation_df$PMID <- as.character(annotation_df$PMID)
  annotation_df$Gene <- as.character(annotation_df$Gene)
  annotation_df$Entrez_ID <- as.character(annotation_df$Entrez_ID)

  genes <- sort(unique(annotation_df$Gene))

  pmid_counts <- sapply(genes, function(g) {
    length(unique(annotation_df$PMID[annotation_df$Gene == g]))
  })

  entrez_map <- sapply(genes, function(g) {
    vals <- unique(annotation_df$Entrez_ID[annotation_df$Gene == g])
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(vals) == 0) NA_character_ else paste(vals, collapse = ";")
  })

  evidence_pmids <- sapply(genes, function(g) {
    paste(unique(annotation_df$PMID[annotation_df$Gene == g]), collapse = ";")
  })

  summary_df <- data.frame(
    Gene = genes,
    Entrez_ID = unname(entrez_map),
    PMID_Count = as.integer(unname(pmid_counts)),
    Supporting_PMIDs = unname(evidence_pmids),
    Phenotype = phenotype,
    Source = "PubMed+PubTator3",
    stringsAsFactors = FALSE
  )

  summary_df <- summary_df[order(-summary_df$PMID_Count, summary_df$Gene), , drop = FALSE]
  rownames(summary_df) <- NULL
  summary_df
}

download_pubmed_pubtator_genes <- function(phenotype) {
  pmids <- search_pubmed_pmids(phenotype)

  if (length(pmids) == 0) {
    cat("❌ No PubMed articles found for:", phenotype, "\n")
    return(list(summary = data.frame(), detailed = data.frame(), metadata = data.frame()))
  }

  metadata <- fetch_pubmed_metadata(pmids)
  annotations <- extract_pubtator_genes(pmids)

  if (nrow(annotations) == 0) {
    cat("❌ No PubTator gene annotations found for:", phenotype, "\n")
    return(list(summary = data.frame(), detailed = data.frame(), metadata = metadata))
  }

  if (nrow(metadata) > 0) {
    annotations <- merge(annotations, metadata, by = "PMID", all.x = TRUE)
  }

  summary_df <- build_gene_summary(annotations, phenotype)

  cat("✅ Total unique genes:", nrow(summary_df), "\n")
  cat("✅ Total annotated gene mentions:", nrow(annotations), "\n")

  list(summary = summary_df, detailed = annotations, metadata = metadata)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("Reliable PubMed -> Gene Downloader (PubMed + PubTator3)\n")
    cat("Usage:   Rscript pubmed.R <phenotype> [output_file]\n")
    cat("Example: Rscript pubmed.R migraine\n")
    quit(status = 1)
  }

  phenotype <- args[1]

  output_dir <- "AllPackagesGenes"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  clean_phenotype <- gsub("[^a-zA-Z0-9_-]", "_", phenotype)

  output_file <- if (length(args) > 1) {
    args[2]
  } else {
    file.path(output_dir, paste0(clean_phenotype, "_pubmed_pubtator.csv"))
  }

  genes_only_file <- file.path(output_dir, paste0(clean_phenotype, "_pubmed_genes.csv"))
  detailed_file <- file.path(output_dir, paste0(clean_phenotype, "_pubmed_pubtator_detailed.csv"))
  metadata_file <- file.path(output_dir, paste0(clean_phenotype, "_pubmed_pubtator_metadata.csv"))

  cat("🎯 Phenotype:        ", phenotype, "\n")
  cat("📁 Output directory: ", output_dir, "\n")
  cat("⏰ Start time:       ", format(Sys.time()), "\n\n")

  result <- download_pubmed_pubtator_genes(phenotype)

  summary_df <- result$summary
  detailed_df <- result$detailed
  metadata_df <- result$metadata

  if (nrow(summary_df) > 0) {
    write.csv(summary_df, output_file, row.names = FALSE)
    write.csv(detailed_df, detailed_file, row.names = FALSE)
    if (nrow(metadata_df) > 0) write.csv(metadata_df, metadata_file, row.names = FALSE)

    genes_only <- data.frame(Gene = unique(summary_df$Gene), stringsAsFactors = FALSE)
    write.csv(genes_only, genes_only_file, row.names = FALSE)

    cat("\n✅ SUCCESS!\n")
    cat("📊 Unique genes:      ", nrow(summary_df), "\n")
    cat("🧬 Gene mentions:     ", nrow(detailed_df), "\n")
    cat("📰 Articles searched: ", if (nrow(metadata_df) > 0) nrow(metadata_df) else length(unique(detailed_df$PMID)), "\n")
    cat("💾 Summary file:      ", output_file, "\n")
    cat("💾 Detailed file:     ", detailed_file, "\n")
    cat("💾 Metadata file:     ", metadata_file, "\n")
    cat("🧬 Genes-only file:   ", genes_only_file, "\n\n")

    top_n <- min(25, nrow(summary_df))
    cat("Top genes by supporting PMID count:\n")
    for (i in seq_len(top_n)) {
      cat(sprintf("  %2d. %-20s PMIDs=%d\n",
                  i, summary_df$Gene[i], summary_df$PMID_Count[i]))
    }
  } else {
    cat("❌ No genes found for:", phenotype, "\n")
    cat("This usually means PubTator returned no recognized gene annotations for the matched papers, or the response format changed again.\n")
  }

  cat("\n⏰ End time: ", format(Sys.time()), "\n")
}

if (!interactive()) {
  main()
}