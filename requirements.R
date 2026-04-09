#!/usr/bin/env Rscript
# Comprehensive R Package Installer for GeneDownloaderScriptsR
# ============================================================
# Installs all required packages for the 14 database gene downloader scripts
# 
# Usage: Rscript requirements.R
# 
# This script analyzes dependencies from all individual R scripts and installs:
# - Core CRAN packages for data manipulation and web APIs
# - Bioconductor packages for biological databases
# - Database-specific packages with API access

cat("🔧 Installing comprehensive R packages for GeneDownloaderScriptsR...\n")
cat("This may take several minutes depending on your system...\n\n")

#' Database-specific package requirements
#' Based on analysis of individual R scripts
database_requirements <- list(
  "pubmed.R" = list(
    cran = c("httr", "xml2", "jsonlite"),
    bioc = c(),
    description = "PubMed literature mining via NCBI E-utilities"
  ),
  
  "omim.R" = list(
    cran = c("httr", "rvest", "dplyr", "stringr"),
    bioc = c(),
    description = "OMIM genetic disorders web scraping"
  ),
  
  "string_db.R" = list(
    cran = c("httr", "jsonlite"),
    bioc = c(),
    description = "STRING protein interactions API"
  ),
  
  "disgenet.R" = list(
    cran = c("dplyr"),
    bioc = c("disgenet2r"),
    description = "DisGeNET gene-disease associations API"
  ),
  
  "clinvar.R" = list(
    cran = c("rentrez", "dplyr", "stringr"),
    bioc = c(),
    description = "ClinVar clinical variants via NCBI E-utilities"
  ),
  
  "reactome_pathways.R" = list(
    cran = c("dplyr", "stringr"),
    bioc = c("ReactomePA", "reactome.db", "org.Hs.eg.db", "AnnotationDbi"),
    description = "Reactome biological pathways"
  ),
  
  "kegg.R" = list(
    cran = c(),
    bioc = c("KEGGREST"),
    description = "KEGG metabolic pathways"
  ),
  
  "hpo.R" = list(
    cran = c("ontologyIndex", "dplyr", "httr", "jsonlite", "utils"),
    bioc = c(),
    description = "Human Phenotype Ontology"
  ),
  
  "gtex.R" = list(
    cran = c("httr", "jsonlite", "dplyr"),
    bioc = c(),
    description = "GTEx gene expression data"
  ),
  
  "uniprot.R" = list(
    cran = c("httr", "xml2", "dplyr"),
    bioc = c(),
    description = "UniProt protein database"
  ),
  
  "opentargets.R" = list(
    cran = c("httr", "jsonlite"),
    bioc = c(),
    description = "Open Targets drug targets GraphQL API"
  ),
  
  "gwasrapidd.R" = list(
    cran = c("dplyr"),
    bioc = c("gwasrapidd"),
    description = "GWAS Catalog associations"
  ),
  
  "gene_ontology.R" = list(
    cran = c("dplyr", "stringr"),
    bioc = c("GO.db", "org.Hs.eg.db", "AnnotationDbi"),
    description = "Gene Ontology functional annotations"
  ),
  
  "string.R" = list(
    cran = c("httr", "jsonlite", "dplyr"),
    bioc = c(),
    description = "Alternative STRING implementation"
  )
)

# Install BiocManager first (required for Bioconductor packages)
if (!require(BiocManager, quietly = TRUE)) {
  cat("📦 Installing BiocManager (required for Bioconductor packages)...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  library(BiocManager)
}

# Aggregate all required packages
all_cran_packages <- unique(unlist(sapply(database_requirements, function(x) x$cran)))
all_bioc_packages <- unique(unlist(sapply(database_requirements, function(x) x$bioc)))

# Additional utility packages for master scripts
additional_cran <- c(
  "readr",      # For CSV reading/writing (download_genes.R)
  "tools",      # For file operations
  "curl"        # For robust HTTP operations
)

# Core CRAN packages (data manipulation, web APIs)
core_cran_packages <- c(
  # Data manipulation
  "dplyr", "stringr", "readr", "tools",
  
  # Web APIs and HTTP
  "httr", "jsonlite", "xml2", "rvest", "curl",
  
  # NCBI access
  "rentrez",
  
  # Ontology analysis
  "ontologyIndex",
  
  # Utilities
  "utils"
)

# Core Bioconductor packages (biological databases)
core_bioc_packages <- c(
  # Gene annotations
  "org.Hs.eg.db", "AnnotationDbi",
  
  # Pathway databases
  "ReactomePA", "reactome.db", "GO.db", "KEGGREST",
  
  # Specialized APIs
  "gwasrapidd", "disgenet2r"
)

#' Print database requirements summary
print_database_summary <- function() {
  cat("📊 Database-Specific Package Requirements:\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  for (script_name in names(database_requirements)) {
    req <- database_requirements[[script_name]]
    cat(sprintf("📄 %s\n", script_name))
    cat(sprintf("   📝 %s\n", req$description))
    
    if (length(req$cran) > 0) {
      cat(sprintf("   📦 CRAN: %s\n", paste(req$cran, collapse = ", ")))
    }
    
    if (length(req$bioc) > 0) {
      cat(sprintf("   🧬 Bioconductor: %s\n", paste(req$bioc, collapse = ", ")))
    }
    
    cat("\n")
  }
}

#' Install CRAN packages with error handling
install_cran_packages <- function() {
  cat("📦 Installing CRAN packages...\n")
  
  // Combine core and database-specific packages
  all_cran <- unique(c(core_cran_packages, all_cran_packages, additional_cran))
  missing_cran <- all_cran[!all_cran %in% rownames(installed.packages())]
  
  if (length(missing_cran) > 0) {
    cat(sprintf("Installing %d missing CRAN packages: %s\n", 
                length(missing_cran), paste(missing_cran, collapse = ", ))
    
    for (pkg in missing_cran) {
      tryCatch({
        cat(sprintf("   Installing %s...", pkg))
        install.packages(pkg, repos = "https://cloud.r-project.org/", dependencies = TRUE)
        cat(" ✅\n")
      }, error = function(e) {
        cat(sprintf(" ❌ Error: %s\n", e$message))
      })
    }
    
    cat("✅ CRAN package installation completed\n")
  } else {
    cat("✅ All required CRAN packages already installed\n")
  }
}

#' Install Bioconductor packages with error handling
install_bioconductor_packages <- function() {
  cat("🧬 Installing Bioconductor packages...\n")
  
  // Ensure BiocManager is available
  library(BiocManager)
  
  // Combine core and database-specific packages
  all_bioc <- unique(c(core_bioc_packages, all_bioc_packages))
  missing_bioc <- all_bioc[!all_bioc %in% rownames(installed.packages())]
  
  if (length(missing_bioc) > 0) {
    cat(sprintf("Installing %d missing Bioconductor packages: %s\n", 
                length(missing_bioc), paste(missing_bioc, collapse = ", ))
    
    for (pkg in missing_bioc) {
      tryCatch({
        cat(sprintf("   Installing %s...", pkg))
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
        cat(" ✅\n")
      }, error = function(e) {
        cat(sprintf(" ❌ Error: %s\n", e$message))
        cat(sprintf("      Try manual installation: BiocManager::install('%s')\n", pkg))
      })
    }
    
    cat("✅ Bioconductor package installation completed\n")
  } else {
    cat("✅ All required Bioconductor packages already installed\n")
  }
}

#' Check installation status of all packages
check_packages <- function() {
  cat("🔍 Checking package installation status...\n")
  
  // Check CRAN packages
  all_cran <- unique(c(core_cran_packages, all_cran_packages, additional_cran))
  cat("📦 CRAN PACKAGES:\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  
  cran_missing <- c()
  for (pkg in all_cran) {
    if (pkg %in% rownames(installed.packages())) {
      cat(sprintf("✅ %s\n", pkg))
    } else {
      cat(sprintf("❌ %s (MISSING)\n", pkg))
      cran_missing <- c(cran_missing, pkg)
    }
  }
  
  // Check Bioconductor packages
  all_bioc <- unique(c(core_bioc_packages, all_bioc_packages))
  cat("\n🧬 BIOCONDUCTOR PACKAGES:\n")
  cat(paste(rep("-", 40), collapse = ""), "\n")
  
  bioc_missing <- c()
  for (pkg in all_bioc) {
    if (pkg %in% rownames(installed.packages())) {
      cat(sprintf("✅ %s\n", pkg))
    } else {
      cat(sprintf("❌ %s (MISSING)\n", pkg))
      bioc_missing <- c(bioc_missing, pkg)
    }
  }
  
  // Summary
  cat("\n📋 INSTALLATION SUMMARY:\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  cat(sprintf("📦 CRAN packages: %d installed, %d missing\n", 
              length(all_cran) - length(cran_missing), length(cran_missing)))
  cat(sprintf("🧬 Bioconductor packages: %d installed, %d missing\n", 
              length(all_bioc) - length(bioc_missing), length(bioc_missing)))
  
  if (length(cran_missing) == 0 && length(bioc_missing) == 0) {
    cat("🎉 All packages successfully installed!\n")
    cat("   Ready to run: Rscript download_genes.R <phenotype>\n")
  } else {
    cat("⚠️  Some packages are missing. Run installation functions to fix.\n")
  }
  
  return(list(
    cran_missing = cran_missing,
    bioc_missing = bioc_missing,
    all_installed = (length(cran_missing) == 0 && length(bioc_missing) == 0)
  ))
}

#' Test package loading (optional verification)
test_package_loading <- function() {
  cat("🧪 Testing package loading...\n")
  
  all_packages <- unique(c(core_cran_packages, all_cran_packages, additional_cran,
                          core_bioc_packages, all_bioc_packages))
  
  failed_loads <- c()
  for (pkg in all_packages) {
    tryCatch({
      library(pkg, character.only = TRUE, quietly = TRUE)
      cat(sprintf("✅ %s loads successfully\n", pkg))
    }, error = function(e) {
      cat(sprintf("❌ %s failed to load: %s\n", pkg, e$message))
      failed_loads <- c(failed_loads, pkg)
    })
  }
  
  if (length(failed_loads) == 0) {
    cat("🎉 All packages load successfully!\n")
  } else {
    cat(sprintf("⚠️  %d packages failed to load\n", length(failed_loads)))
  }
  
  return(failed_loads)
}

# Main execution
main <- function() {
  cat("🧬 GeneDownloaderScriptsR Package Installation\n")
  cat("=============================================\n")
  
  cat("📊 System Information:\n")
  cat(sprintf("   R version: %s\n", R.version.string))
  cat(sprintf("   Platform: %s\n", R.version$platform))
  cat(sprintf("   OS: %s\n", Sys.info()["sysname"]))
  cat("\n")
  
  // Show database summary
  print_database_summary()
  
  // Check current status
  status <- check_packages()
  
  // Install missing packages
  if (!status$all_installed) {
    cat("\n🔧 Installing missing packages...\n")
    install_cran_packages()
    install_bioconductor_packages()
    
    // Re-check status
    cat("\n🔍 Re-checking installation status...\n")
    final_status <- check_packages()
    
    if (final_status$all_installed) {
      cat("\n🎉 Installation completed successfully!\n")
      cat("   You can now run: Rscript download_genes.R <phenotype>\n")
    } else {
      cat("\n⚠️  Some packages still missing. Manual intervention may be required.\n")
    }
  }
  
  cat("\n📚 Next Steps:\n")
  cat("   1. Test installation: Rscript requirements.R --test\n")
  cat("   2. Run gene download: Rscript download_genes.R migraine\n")
  cat("   3. Run analysis: python download_genes_analysis.py migraine\n")
  cat("\n")
}

// Handle command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  if (args[1] == "--test") {
    cat("🧪 Testing package installation and loading...\n")
    check_packages()
    test_package_loading()
  } else if (args[1] == "--check") {
    check_packages()
  } else if (args[1] == "--summary") {
    print_database_summary()
  } else {
    cat("Usage: Rscript requirements.R [--test|--check|--summary]\n")
    cat("   (no args): Install all packages\n")
    cat("   --test: Test package loading\n")
    cat("   --check: Check installation status\n")
    cat("   --summary: Show database requirements\n")
  }
} else {
  // Default: run main installation
  main()
}

#' Install Bioconductor packages
install_bioconductor_packages <- function() {
  cat("🧬 Installing Bioconductor packages...\n")
  
  # Install BiocManager if not available
  if (!require(BiocManager, quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager")
  }
  
  library(BiocManager)
  
  missing_bioc <- bioconductor_packages[!bioconductor_packages %in% rownames(installed.packages())]
  
  if (length(missing_bioc) > 0) {
    cat("Installing missing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
    cat("✅ Bioconductor packages installed successfully\n")
  } else {
    cat("✅ All Bioconductor packages already installed\n")
  }
}

#' Check package installation status
check_packages <- function() {
  cat("🔍 Checking package installation status...\n")
  cat("\nCRAN PACKAGES:\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  
  for (pkg in cran_packages) {
    status <- if (pkg %in% rownames(installed.packages())) "✅ Installed" else "❌ Missing"
    cat(sprintf("%-15s: %s\n", pkg, status))
  }
  
  cat("\nBIOCONDUCTOR PACKAGES:\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  
  for (pkg in bioconductor_packages) {
    status <- if (pkg %in% rownames(installed.packages())) "✅ Installed" else "❌ Missing"
    cat(sprintf("%-15s: %s\n", pkg, status))
  }
  
  # Check if all packages are installed
  all_packages <- c(cran_packages, bioconductor_packages)
  missing_packages <- all_packages[!all_packages %in% rownames(installed.packages())]
  
  cat("\nSUMMARY:\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  cat("Total packages required:", length(all_packages), "\n")
  cat("Packages installed:", length(all_packages) - length(missing_packages), "\n")
  cat("Packages missing:", length(missing_packages), "\n")
  
  if (length(missing_packages) == 0) {
    cat("🎉 All packages are installed and ready!\n")
    return(TRUE)
  } else {
    cat("⚠️  Missing packages:", paste(missing_packages, collapse = ", "), "\n")
    return(FALSE)
  }
}

#' Install all required packages
install_all_packages <- function() {
  cat("🚀 Installing all required packages for Gene Downloader Scripts\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  install_cran_packages()
  cat("\n")
  install_bioconductor_packages()
  cat("\n")
  
  # Final check
  if (check_packages()) {
    cat("\n🎉 All packages installed successfully! You can now run the gene downloader scripts.\n")
  } else {
    cat("\n❌ Some packages failed to install. Please check the error messages above.\n")
  }
}

#' Main function for command line usage
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Gene Downloader R Package Requirements\n")
    cat("Usage: Rscript requirements.R [command]\n")
    cat("\nCommands:\n")
    cat("  install   - Install all required packages\n")
    cat("  check     - Check installation status\n")
    cat("  cran      - Install only CRAN packages\n")
    cat("  bioc      - Install only Bioconductor packages\n")
    cat("  list      - List all required packages\n")
    cat("\nExample: Rscript requirements.R install\n")
    quit(status = 0)
  }
  
  command <- args[1]
  
  switch(command,
    "install" = install_all_packages(),
    "check" = check_packages(),
    "cran" = install_cran_packages(),
    "bioc" = install_bioconductor_packages(),
    "list" = {
      cat("CRAN packages required:\n")
      cat(paste("  -", cran_packages), sep = "\n")
      cat("\nBioconductor packages required:\n")
      cat(paste("  -", bioconductor_packages), sep = "\n")
    },
    {
      cat("Unknown command:", command, "\n")
      cat("Use 'Rscript requirements.R' for help\n")
      quit(status = 1)
    }
  )
}

# Run if called as script
if (!interactive()) {
  main()
}
