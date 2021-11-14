# Installing Pacman as a Package Manager if not installed, then load
# Packages to install - Installing only Pacman as a Package Manager

## Load or install and load all libraries
install.packages("pacman", repos='http://cran.us.r-project.org')
install.packages("BiocManager", repos='http://cran.us.r-project.org')
install.packages("XML", repos='http://cran.us.r-project.org')

library("pacman")
library("BiocManager")
library("XML")

# List of CRAN packages to either Load, or Install and Load
pacman::p_load(stringr, data.table)

# List of Bioconductor packages to either Load, or Install and Load
pacman::p_load(rtracklayer, TxDb.Hsapiens.UCSC.hg38.knownGene, Repitools,
               Organism.dplyr, GenomicFeatures)