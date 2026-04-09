# R Package Requirements for Gene Downloader Scripts

## CRAN Packages (install via install.packages())
- dplyr          # Data manipulation and transformation
- readr          # Reading and writing data files
- tools          # R development tools
- httr           # HTTP requests for API calls
- httr2          # Modern HTTP requests (alternative to httr)
- jsonlite       # JSON data parsing
- rvest          # Web scraping
- xml2           # XML parsing
- stringr        # String manipulation
- tidyr          # Data tidying and reshaping
- purrr          # Functional programming tools
- curl           # HTTP requests (alternative backend)
- RCurl          # Alternative HTTP client
- data.table     # Fast data manipulation
- magrittr       # Pipe operators
- tibble         # Modern data frames
- vroom          # Fast file reading
- rentrez        # NCBI E-utilities access
- easyPubMed     # PubMed data extraction
- europepmc      # Europe PMC API access
- geneLenDataBase # Gene length information
- GenomicRanges  # Genomic ranges manipulation

## Bioconductor Packages (install via BiocManager::install())
- gwasrapidd     # GWAS Catalog API access
- biomaRt        # Biomart database access
- GO.db          # Gene Ontology database
- org.Hs.eg.db  # Human gene annotations
- AnnotationDbi  # Annotation database interface
- KEGGREST       # KEGG pathway REST API
- ReactomePA     # Reactome pathway analysis
- reactome.db    # Reactome database
- GenomicFeatures # Genomic features extraction
- TxDb.Hsapiens.UCSC.hg38.knownGene # Human genome annotations
- BSgenome.Hsapiens.UCSC.hg38 # Human genome sequence
- VariantAnnotation # Variant annotation tools
- gwascat        # GWAS catalog data access
- SNPlocs.Hsapiens.dbSNP144.GRCh38 # SNP locations
- MafDb.1Kgenomes.phase3.hs37d5 # Minor allele frequencies
- PolyPhen.Hsapiens.dbSNP131 # PolyPhen predictions
- SIFT.Hsapiens.dbSNP137    # SIFT predictions
- LDlinkR        # LD Link API access
- haploR         # HaploReg API access
- rsnps          # SNP database access
- myvariant      # MyVariant.info API access

## Installation Commands

### Install all packages automatically:
```r
# Run the requirements script
Rscript requirements.R install
```

### Manual installation:
```r
# Install CRAN packages
cran_packages <- c("dplyr", "readr", "tools", "httr", "httr2", "jsonlite", "rvest", "xml2", "stringr", 
                   "tidyr", "purrr", "curl", "RCurl", "data.table", "magrittr", "tibble", "vroom",
                   "rentrez", "easyPubMed", "europepmc", "geneLenDataBase", "GenomicRanges")
install.packages(cran_packages)

# Install Bioconductor packages
if (!require(BiocManager, quietly = TRUE)) install.packages("BiocManager")
bioc_packages <- c("gwasrapidd", "biomaRt", "GO.db", "org.Hs.eg.db", "AnnotationDbi", "KEGGREST", 
                   "ReactomePA", "reactome.db", "GenomicFeatures", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                   "BSgenome.Hsapiens.UCSC.hg38", "VariantAnnotation", "gwascat", 
                   "SNPlocs.Hsapiens.dbSNP144.GRCh38", "MafDb.1Kgenomes.phase3.hs37d5",
                   "PolyPhen.Hsapiens.dbSNP131", "SIFT.Hsapiens.dbSNP137", "LDlinkR", 
                   "haploR", "rsnps", "myvariant")
BiocManager::install(bioc_packages)
```

## Package Usage by Script

### download_gene2.R (main script)
- dplyr, readr, tools, httr, jsonlite

### gwasrapidd.R
- gwasrapidd, dplyr

### gwas_rest_api.R
- httr, jsonlite

### opentargets.R
- httr, jsonlite

### disgenet.R
- httr, rvest, xml2

### pubmed.R
- httr, jsonlite

### biomart.R
- biomaRt, dplyr, stringr

### gene_ontology.R
- GO.db, org.Hs.eg.db, AnnotationDbi, dplyr, stringr

### kegg_pathways.R
- KEGGREST, org.Hs.eg.db, AnnotationDbi, dplyr, stringr

### reactome_pathways.R
- ReactomePA, reactome.db, org.Hs.eg.db, AnnotationDbi, dplyr, stringr

### curated_genes.R
- No additional packages (uses base R)
