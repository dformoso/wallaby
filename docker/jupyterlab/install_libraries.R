# Enable multithreading when possible (library dependent)
options(Ncpus = parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS=toString(parallel::detectCores()))
Sys.setenv(OMP_THREAD_LIMIT=toString(parallel::detectCores()))
Sys.setenv(OMP_NUM_THREADS=parallel::detectCores())
Sys.setenv(OMP_THREAD_LIMIT=parallel::detectCores())

# Installing tricky packages manually
download.file("https://bioconductor.org/packages/3.12/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz", 
               destfile = "BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz", 
               method="curl")
install.packages("BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz", repos = NULL, type = "source")

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
               data.table,
               kableExtra,
               knitr,
               IRdisplay)

# List of Bioconductor packages to either Load, or Install and Load
pacman::p_load(GenomicFeatures,
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
               Repitools)