#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Variables
srr_name <- toString(args[1])
donor_name <- toString(args[2])
recipient_name <- toString(args[3])
donor_ref_genome <- toString(args[4])
recipient_ref_genome <- toString(args[5])
inputs_folder <- toString(args[6])
min_num_crossings <- as.integer(args[7])
min_num_reads <- as.integer(args[8])

# Select crossings to import and order of visualization
crossings <- c("1Md_2Mr", "1Md_1Mr", "2Md_1Mr")

## Load or install and load all libraries
library("pacman", character.only = TRUE)

# List of CRAN packages to either Load, or Install and Load
pacman::p_load(stringr, data.table, install = FALSE)

# List of Bioconductor packages to either Load, or Install and Load
pacman::p_load(rtracklayer, TxDb.Hsapiens.UCSC.hg38.knownGene, Repitools,
               Organism.dplyr, install = FALSE)

# Load all .bed files for all srrs created by the cromwell workflow
load_beds <- function(srr_names, crossings, name) {
    # Define variables to hold all srrs
    granges_all_srrs <- list()

    # Loop over each srr
    for (srr in srr_names) {

        # For each single srr
        granges <- list()
        granges_crossings <- list()

        # Load all .bed files created by the cromwell workflow
        for (cross in crossings){
            # list all recipient files
            file <- list.files(inputs_folder,
                                      pattern=paste(srr, 
                                                    '-to-',
                                                    name,
                                                    "_", cross, 
                                                    "_filtered.bed", 
                                                    sep = ""), 
                                      recursive = TRUE, 
                                      full.names = TRUE)

            # check whether the file exists and add if it does
            if (!identical(file, character(0))) {
                granges[[cross]] <- import(file)
                granges_crossings[[cross]] <- cross
            }

            # add granges for this srr to the _all object
            granges_all_srrs[[srr]] <- granges
        }
    }
    
    return (granges_all_srrs)
}

# Function to Create a Table mapping ranges of overlapping paired-end crossings 
summary_table_recipient <- function(granges, 
                                    granges_labels, 
                                    min_num_crossings,
                                    min_num_reads, 
                                    src) {

    # convert all granges to dataframes
    granges_df <- lapply(granges, annoGR2DF)
    # assign all granges labels (crossing names) as each dataframe's name
    names(granges_df) <- granges_labels
    # flatten all dataframes into a single dataframe, 
    # and use their's respective name as an identifier in a new column named 'crossing'
    merged_df <- bind_rows(granges_df, .id = "crossing")
    # convert the data.frame to a data.table
    merged_dt <- as.data.table(merged_df)

    # assign each row to a group, based on whether their ranges overlap
    merged_dt[,group := { ir <- IRanges(start, end); subjectHits(findOverlaps(ir, reduce(ir))) }, by = chr]
    # aggregate results by group, and add additional aggregated columns
    merged_final <- merged_dt[, list(start = min(start), 
                                     stop = max(end), 
                                     num_crossings = length(unique(list(crossing)[[1]])),
                                     unique_crossings = list(unique(crossing)),
                                     num_reads = length(list(name)[[1]])
                                     ), by = list(group, chr)]
    
    # filter results using a minimum number of reads and/or crossings
    merged_final <- merged_final[merged_final[, num_reads > (min_num_reads - 1)]]
    merged_final <- merged_final[merged_final[, num_crossings > (min_num_crossings - 1)]]
    
    # there might be no matches, so check for that before looking for a gene name
    if (nrow(merged_final) != 0) {
        # use the src database to look for the gene names or each crossing overlap region
        # then, add it as a new column
        merged_final$gene_name <- apply(merged_final, 1, FUN = function(x) toString(
            unique(unlist(suppressWarnings(annoGR2DF(
                                    transcripts(src, 
                                                 filter=~(GRangesFilter(
                                                     GenomicRanges::GRanges(
                                                         paste(toString(x["chr"]), ":", 
                                                               as.integer(x["start"]), "-", 
                                                               as.integer(x["stop"]), sep = "")))), 
                                                 columns=c("symbol")))$symbol)))))
                                        
        # delete the 'group' column
        merged_final <- merged_final[, !"group"]
        # add an ID to each row
        merged_final <- merged_final[, id := .I]
        setcolorder(merged_final, c("id", "chr", "start", "stop", 
                                    "num_crossings", "unique_crossings", 
                                    "num_reads", "gene_name"))
                                       
    # if there are no matches, write an NA row
    } else {
            return(data.table(id = "<NA>",
                              chr = "<NA>", 
                              start = 0, 
                              stop = 0, 
                              num_crossings = 0, 
                              unique_crossings = "<NA>", 
                              num_reads = 0, 
                              gene_name = "<NA>")
                  )
        }

    return(merged_final)
}
                                       
# Function to Create a Summary Table for many SRRs
srrs_summary_table_recipient <- function(granges_list, 
                                         min_num_crossings,
                                         min_num_reads, 
                                         src) {

    # store all tables in a single list
    srrs_list <- list()
    
    # iterate over all granges
    for (srr in names(granges_list)) {

        # create a summary table for each granges object
        granges <- unname(granges_list[srr][[1]])
        
        # check if there are valid files for the given SRRs
        if (length(granges) != 0) {
            granges_labels <- str_split(names(granges_list[srr][[1]]), " ")
            summary_table <- summary_table_recipient (granges = granges,
                                                      granges_labels = granges_labels,
                                                      min_num_crossings = min_num_crossings, 
                                                      min_num_reads = min_num_reads, 
                                                      src = src)

            # add a column for the srr
            summary_table$srr <- srr
            # add the table to the list
            srrs_list[[srr]] <- summary_table
        } else {
            srrs_list[[srr]] <- data.table(srr = srr,
                                           id = "<NA>",
                                           chr = "<NA>", 
                                           start = 0, 
                                           stop = 0, 
                                           num_crossings = 0, 
                                           unique_crossings = "<NA>", 
                                           num_reads = 0, 
                                           gene_name = "<NA>")
        }
    }

    # check if there are valid files for the given SRRs
    if (length(srrs_list) != 0) {
        # bind all the tables by row
        srrs_summary_table <- do.call(rbind, c(srrs_list, fill=TRUE))
        # reorder the table
        setcolorder(srrs_summary_table, c("srr", "id", "chr", "start", "stop", 
                                          "num_crossings", "unique_crossings", 
                                          "num_reads", "gene_name"))

        return (srrs_summary_table)
        } else {
            print("No files matching SRRs")
            return
    }
}
                                        
# Load all .bed files for all srrs created by the cromwell workflow
donor_granges_all_srrs <- load_beds(srr_name, crossings, donor_name)
recip_granges_all_srrs <- load_beds(srr_name, crossings, recipient_name)
                                        
# Create a sqlite database from TxDb and corresponding Org packages
# The database provides a convenient way to map between gene, transcript, and protein identifiers.
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
                                        
# Aggregated view of all overlapping crossings for the potential recipient for all srrs
srrs_summary_table <- srrs_summary_table_recipient (granges_list = recip_granges_all_srrs, 
                                                    min_num_crossings = min_num_crossings,
                                                    min_num_reads = min_num_reads, 
                                                    src = src)
                                        
fwrite(srrs_summary_table, paste(srr_name, "_putative_insertion_table.csv", sep = ""))