# Enable multithreading when possible (library dependent)
options(Ncpus = parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS=toString(parallel::detectCores()))
Sys.setenv(OMP_THREAD_LIMIT=toString(parallel::detectCores()))
Sys.setenv("OMP_NUM_THREADS"=parallel::detectCores())
Sys.setenv("OMP_THREAD_LIMIT"=parallel::detectCores())

# Installing Pacman as a Package Manager if not installed, then load
# Packages to install - Installing only Pacman as a Package Manager
cran_packages = c("pacman", "BiocManager", "XML")

# Load Pacman, or Install and Load if never installed
package.check <- lapply(
  cran_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)

# List of CRAN packages to either Load, or Install and Load
pacman::p_load(dplyr, 
               ggplot2,
               shiny,
               shinyLP,
               DT, 
               ggrepel, 
               tidyr,, 
               data.table)

# List of Bioconductor packages to either Load, or Install and Load
pacman::p_load(GenomicAlignments, 
               Rsubread, 
               Rsamtools,
               bamsignals, 
               rtracklayer,
               GenomicRanges, 
               BSgenome.Hsapiens.UCSC.hg38,
               TxDb.Hsapiens.UCSC.hg38.knownGene, 
               regioneR,
               karyoploteR, 
               seqinr,
               Repitools)