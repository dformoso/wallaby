# Enable multithreading when possible (library dependent)
options(Ncpus = parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS=toString(parallel::detectCores()))
Sys.setenv(OMP_THREAD_LIMIT=toString(parallel::detectCores()))
Sys.setenv(OMP_NUM_THREADS=parallel::detectCores())
Sys.setenv(OMP_THREAD_LIMIT=parallel::detectCores())

# Installing Pacman as a Package Manager if not installed, then load
# Packages to install - Installing only Pacman as a Package Manager
cran_packages = c("pacman", "BiocManager", "XML")

# Load Pacman, or Install and Load if never installed
package.check <- lapply(
  cran_packages,
  FUN = function(x) {
      install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
  }
)

# Options for tricky packages a.k.a.: BSgenome.Hsapiens.UCSC.hg38
options(timeout = 300)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

# List of CRAN packages to either Load, or Install and Load
pacman::p_load(dplyr, 
               ggplot2,
               shiny,
               shinyLP,
               DT, 
               ggrepel, 
               tidyr,, 
               data.table,
               kableExtra,
               knitr,
               IRdisplay,
               shiny)

# List of Bioconductor packages to either Load, or Install and Load
pacman::p_load(BSgenome,
               BSgenome.Hsapiens.UCSC.hg38,
               GenomicFeatures,
               GenomicAlignments, 
               Rsubread, 
               Rsamtools,
               bamsignals, 
               rtracklayer,
               GenomicRanges,
               TxDb.Hsapiens.UCSC.hg38.knownGene, 
               regioneR,
               karyoploteR, 
               seqinr,
               Repitools,
               org.Hs.eg.db, 
               Organism.dplyr,
               Gviz, 
               Biostrings)