## Setup Environment
setup_enviroment <- function() {
    # Enable multithreading when possible (library dependent)
    options(Ncpus = parallel::detectCores())
    Sys.setenv(OMP_NUM_THREADS=toString(parallel::detectCores()))
    Sys.setenv(OMP_THREAD_LIMIT=toString(parallel::detectCores()))
    Sys.setenv(OMP_NUM_THREADS=parallel::detectCores())
    Sys.setenv(OMP_THREAD_LIMIT=parallel::detectCores())

    ## Load or install and load all libraries
    suppressPackageStartupMessages(library("pacman", character.only = TRUE))

    # List of CRAN packages to either Load, or Install and Load
    pacman::p_load(dplyr, ggplot2, shiny, shinyLP, DT,  ggrepel, tidyr, IRdisplay, repr,
                   data.table, kableExtra, knitr, IRdisplay, install = FALSE)

    # List of Bioconductor packages to either Load, or Install and Load
    pacman::p_load(BSgenome, BSgenome.Hsapiens.UCSC.hg38, GenomicFeatures, 
                   GenomicAlignments,  Rsubread,  Rsamtools, bamsignals,  
                   rtracklayer, GenomicRanges, org.Hs.eg.db, Organism.dplyr,
                   TxDb.Hsapiens.UCSC.hg38.knownGene, regioneR, karyoploteR,
                   seqinr, Repitools, Gviz, Biostrings, install = FALSE)
    }