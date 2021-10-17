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
    pacman::p_load(dplyr, ggplot2, shiny, shinyLP, DT,  ggrepel,  tidyr, IRdisplay, repr, stringr,
                   data.table, kableExtra, knitr, IRdisplay, install = FALSE)

    # List of Bioconductor packages to either Load, or Install and Load
    pacman::p_load(BSgenome, BSgenome.Hsapiens.UCSC.hg38, GenomicFeatures, 
                   GenomicAlignments,  Rsubread,  Rsamtools, bamsignals, ensembldb,
                   rtracklayer, GenomicRanges, org.Hs.eg.db, Organism.dplyr,
                   TxDb.Hsapiens.UCSC.hg38.knownGene, regioneR, karyoploteR,
                   seqinr, Repitools, Gviz, Biostrings, install = FALSE)
}

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

# Load .bam files created by the cromwell workflow
load_bams <- function(srr_names, crossings, name) {
    # Define variables to hold all srrs
    bams_all_srrs <- list()

    # Loop over each srr
    for (srr in srr_names) {

        bams <- list()
        bam_crossings <- list()

        for (cross in crossings){
            # list all donor files
            file <- list.files(inputs_folder,
                                      pattern=paste(srr, 
                                                    '-to-',
                                                    name,
                                                    "_", cross, 
                                                    "_filtered.bam$", 
                                                    sep = ""), 
                                      recursive = TRUE, 
                                      full.names = TRUE)

            # check whether they exist and add if they do
            if (!identical(file, character(0))) {
                bams[[length(bams) + 1]] <- file
                bam_crossings[[length(bam_crossings) + 1]] <- cross
            }

            # add _granges for this srr to the _all object
            bams_all_srrs[[srr]] <- bams
        }
    }
    
    return (bams_all_srrs)
}

# Function to Create a Table mapping ranges of overlapping paired-end crossings 
summary_table_recipient <- function(granges, 
                                    granges_labels, 
                                    min_num_crossings,
                                    min_num_reads, 
                                    src,
                                    Hsapiens) {

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
                                        
        # use the Hsapiens databse to look up the DNA sequence for each crossing overlap region
        # then, add it as a new column
        merged_final$sequence <- apply(merged_final, 1, FUN = function(x) toString(getSeq(Hsapiens, 
                                                                           toString(x["chr"]), 
                                                                           start = as.integer(x["start"]), 
                                                                           end = as.integer(x["stop"]))))
        # delete the 'group' column
        merged_final <- merged_final[, !"group"]
        # add an ID to each row
        merged_final <- merged_final[, id := .I]
        setcolorder(merged_final, c("id", "chr", "start", "stop", 
                                    "num_crossings", "unique_crossings", 
                                    "num_reads", "gene_name", "sequence"))
                                       
    # if there are no matches, write an NA row
    } else {
            return(data.table(id = "<NA>",
                              chr = "<NA>", 
                              start = 0, 
                              stop = 0, 
                              num_crossings = 0, 
                              unique_crossings = "<NA>", 
                              num_reads = 0, 
                              gene_name = "<NA>", 
                              sequence = "<NA>")
                  )
        }

    return(merged_final)
}
                                       

# Function to Create a Table mapping ranges of overlapping paired-end crossings 
summary_table_donor <- function(granges, 
                                granges_labels, 
                                min_num_crossings, 
                                min_num_reads) {
    
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
    merged_final <- merged_dt[, list(start=min(start), 
                                     stop=max(end), 
                                     num_crossings=length(unique(list(crossing)[[1]])),
                                     unique_crossings=list(unique(crossing)),
                                     num_reads=length(list(name)[[1]])
                                     ), by=list(group,chr)]

    # filter results using a minimum number of reads and/or crossings
    merged_final <- merged_final[merged_final[, num_reads > (min_num_reads - 1)]]
    merged_final <- merged_final[merged_final[, num_crossings > (min_num_crossings - 1)]]

        # there might be no matches, so check for that before looking for a gene name
    if (nrow(merged_final) != 0) {
        # delete the 'group' column
        merged_final <- merged_final[, !"group"]
        # add an ID to each row
        merged_final <- merged_final[, id := .I]
        setcolorder(merged_final, c("id", "chr", "start", "stop", 
                                    "num_crossings", "unique_crossings", 
                                    "num_reads"))
    } else {
            return(data.table(id = "<NA>",
                              chr = "<NA>", 
                              start = 0, 
                              stop = 0, 
                              num_crossings = 0, 
                              unique_crossings = "<NA>", 
                              num_reads = 0)
                  )
        }

    return(merged_final)
}
                                       
# Function to Create a Summary Table for many SRRs
srrs_summary_table_recipient <- function(granges_list, 
                                         min_num_crossings,
                                         min_num_reads, 
                                         src,
                                         Hsapiens) {

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
                                                      min_num_crossings = 1, 
                                                      min_num_reads = min_num_reads, 
                                                      src = src, 
                                                      Hsapiens = Hsapiens)

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
                                           gene_name = "<NA>", 
                                           sequence = "<NA>")
        }
    }

    # check if there are valid files for the given SRRs
    if (length(srrs_list) != 0) {
        # bind all the tables by row
        srrs_summary_table <- do.call(rbind, c(srrs_list, fill=TRUE))
        # reorder the table
        setcolorder(srrs_summary_table, c("srr", "id", "chr", "start", "stop", 
                                          "num_crossings", "unique_crossings", 
                                          "num_reads", "gene_name", "sequence"))

        # stylize the output
        srrs_summary_table %>%
            kable("html") %>%
            kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
            column_spec(10, width = "30em", width_max = "30em") %>%
            as.character() %>%
            display_html()

        return (srrs_summary_table)
        } else {
            print("No files matching SRRs")
            return
    }
}
                                       
# Function to Create a Summary Table for many SRRs
srrs_summary_table_donor <- function(granges_list, 
                                     min_num_crossings,
                                     min_num_reads) {

    # store all tables in a single list
    srrs_list <- list()
    
    # iterate over all granges
    for (srr in names(granges_list)) {
        
        # create a summary table for each granges object
        granges <- unname(granges_list[srr][[1]])
        
        # check if there are valid files for the given SRRs
        if (length(granges) != 0) {
            granges_labels <- str_split(names(granges_list[srr][[1]]), " ")
            summary_table <- summary_table_donor (granges = granges,
                                                  granges_labels = granges_labels,
                                                  min_num_crossings = 1, 
                                                  min_num_reads = min_num_reads)

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
                                           num_reads = 0)
        }
    }
    
    # check if there are valid files for the given SRRs
    if (length(srrs_list) != 0) {
        # bind all the tables by row
        srrs_summary_table <- do.call(rbind, c(srrs_list, fill=TRUE))
        # reorder the table
        setcolorder(srrs_summary_table, c("srr", "id", "chr", "start", "stop", 
                                          "num_crossings", "unique_crossings", 
                                          "num_reads"))

        # stylize the output
        srrs_summary_table %>%
            kable("html") %>%
            kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
            as.character() %>%
            display_html()
    
    return (srrs_summary_table)
    } else {
        print("No files matching SRRs")
        return
    }
}
                                       
# Function to create a visualization for specific overlap regions
plot_reads_region <- function(srr, id = 1, crossings_table_recipient, recip_bams,
                              extend_left, extend_right, 
                              ref_genome, donor_name, recipient_name) {
    
    
    options(ucscChromosomeNames=FALSE)
    
    # extract chromosome, start, and end positions from the given overlap table
    chr <- toString(crossings_table_recipient[id,]$chr)
    start <- toString(crossings_table_recipient[id,]$start - extend_left) 
    end <- toString(crossings_table_recipient[id,]$stop + extend_right)
    crossings <- as.list(strsplit(crossings_table_recipient[id,]$unique_crossings[[1]], ","))
    gene_name <- toString(crossings_table_recipient[id,]$gene_name)
    graph_title <- paste("Donor: ", donor_name, "  -  ",
                         "Recipient: ", recipient_name, "  -  ",
                         "SRR: ", srr, "  -  ",
                         "ID: ", toString(crossings_table_recipient[id,]$id), "  -  ",
                         "Num Reads: ", toString(crossings_table_recipient[id,]$num_reads), "  -  ",
                         "Gene Name: ", if (gene_name == "") "No Match" else gene_name)
    
    # only shows bams with overlaps
    bams <- list()
    for (crossing in crossings) {
        bams <- c(bams, recip_bams[grepl(crossing, recip_bams, fixed = TRUE)])
    }
    
    # create a track which holds a schematic display of a chromosome
    i_track <- IdeogramTrack(genome = "hg38", chromosome = chr, 
                             from = as.numeric(start) - extend_left, 
                             to = as.numeric(end) + extend_right, 
                             showId = TRUE,  showBandId = TRUE, 
                             cex = 3, cex.bands = 1)
    
    # create a track which display the genomics axis
    g_track <- GenomeAxisTrack(showId = TRUE, labelPos = "alternating", cex = 2)
    
    # create a track which holds genes and exons names
    gr_track <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                genome = 'hg38', 
                                chromosome = chr, name = "Exons", 
                                background.title = "red", fill = "orange",
                                showID = TRUE, collapse = TRUE,
                                transcriptAnnotation = NULL, exonAnnotation = 'gene',
                                showExonId = TRUE, mergeGroups = TRUE,
                                showOverplotting = TRUE,
                                max.height = 10, min.height = 0.01)
    
    # create a track which holds the reads, coloring mismatches and indels
    a_tracks <- mapply(function(x, y) { AlignmentsTrack(x, name = y, isPaired = TRUE, 
                                                        stacking = 'full', max.height = 10,
                                                        chromosome = chr, min.height = 0.01, 
                                                        background.title = "blue", fill="black",
                                                        alpha = 0.90, alpha.mismatch = 1,
                                                        type = "pileup", showMismatches = TRUE, 
                                                        showIndels = TRUE, col.mates = "purple", 
                                                        from = as.numeric(start) - extend_left, 
                                                        to = as.numeric(end) + extend_right)
                                      }, bams, crossings)
    
    # create a track which holds each letter
    s_track <- SequenceTrack(readDNAStringSet(recipient_ref_genome), chromosome = chr, min.width = 0.1, cex = 0.5)
    
    # plot all tracks together
    plotTracks(c(i_track, g_track, gr_track, a_tracks, s_track), 
               chromosome = chr, col.main = "black",
               from = as.numeric(start), to = as.numeric(end), main = graph_title,
               extend.left = extend_left, extend.right = extend_right, margin=c(30,30,-10,30), 
               just.group = 'above', cex.title = 2, rotation.title = 0, 
               title.width = 8, sizes = c(0.5, 0.7, 0.75, replicate(length(crossings), 1), 0.2))
}

# Function to create a visualization of the donor crossings
create_viz_donor <- function(ref_genome = "hg38", 
                             granges, 
                             granges_labels,
                             title_prepend = "") {

    # reverse order of granges and granges_labels so that they plot in the right order
    # as plotKaryotype reverses it again
    granges <- rev(granges)
    granges_labels <- rev(granges_labels)
    
    # Set up plot parameters
    plot.type <- 4
    tracks <- length(granges)
    track_sep <- 0.05
    track_width <- 1 / (tracks) - track_sep
    window.size <- 10
    title <- paste(title_prepend, "- window size (in bases): ", window.size)
    pp <- getDefaultPlotParams(plot.type=plot.type)
    pp$leftmargin <- 0.17
    
    # Create a custom granges object for the Donor Reference Genome
    summary_fasta <- summary(read.fasta(ref_genome))
    total_genome_length <- as.integer(summary_fasta[, "Length"])
    seqname <- unique(as.character(seqnames(granges[[1]])))
    custom.genome <- toGRanges(data.frame(chr = c(seqname), 
                                          start = c(1), 
                                          end = c(total_genome_length)))
    
    # Create the object for plotting
    kp <- plotKaryotype(genome = custom.genome,
                        plot.type = plot.type, 
                        plot.params = pp, 
                        labels.plotter = NULL, 
                        main = title,
                        cex = 2)
    kpAddBaseNumbers(kp, tick.dist = window.size * 50, add.units = TRUE, cex = 2) 

    # Create the tracks in the plot, depending on how many tracks there are
    track_no <- 0
    for (grange in granges) {
        track_no <- track_no + 1
        
        r0 <- (track_no-1) * track_width + (track_no-1) * track_sep
        r1 <- track_no * track_width + (track_no-1) * track_sep
        
        kp <- suppressWarnings(kpPlotDensity(kp, data = grange, 
                                             window.size = window.size, 
                                             col = "blue", r0 = r0, r1 = r1))
        
        kpAxis(kp, ymax = kp$latest.plot$computed.values$max.density, cex = 2, r0 = r0, r1 = r1)
        kpAddLabels(kp, labels = granges_labels[track_no], 
                    r0 = r0, r1 = r1, label.margin = 0.07, cex = 2)
    }
}

# Function to create a visualization of the recipient crossings
create_viz_recipient <- function(graph_type = "recipient", 
                                 ref_genome = "hg38", 
                                 granges,
                                 granges_labels, 
                                 title_prepend = "") {

    # reverse order of granges and granges_labels so that they plot in the right order
    # as plotKaryotype reverses it again
    granges <- rev(granges)
    granges_labels <- rev(granges_labels)
    
    # Set up plot parameters
    plot.type <- 4
    tracks <- length(granges)
    track_sep <- 0.05
    track_width <- 1 / tracks - track_sep
    genome = "hg38"
    window.size <- 1e6
    title <- paste(title_prepend, "- window size (in bases): ", window.size)
    pp <- getDefaultPlotParams(plot.type=plot.type)
    pp$leftmargin <- 0.17
    pp$topmargin <- 30
    pp$bottommargin <- 30
    
    # Create the object for plotting
    kp <- plotKaryotype(genome = ref_genome,
                        plot.type = plot.type, 
                        plot.params = pp, 
                        labels.plotter = NULL, 
                        main = title,
                        cex = 2,
                        chromosomes=c("chr1", "chr2", "chr3", "chr4", "chr5", 
                                      "chr6", "chr7", "chr8", "chr9", "chr10",
                                      "chr11", "chr12", "chr13", "chr14", "chr15", 
                                      "chr16", "chr17", "chr18", "chr19", "chr20",
                                      "chr21", "chr22", "chrX", "chrY"))

    kpAddChromosomeNames(kp, srt = 90, cex = 2) 
    
    # Create the tracks in the plot, depending on how many tracks there are
    track_no <- 0
    for (grange in granges) {
        track_no <- track_no + 1

        r0 <- (track_no-1) * track_width + (track_no-1) * track_sep
        r1 <- track_no * track_width + (track_no-1) * track_sep
        
        kp <- suppressWarnings(kpPlotDensity(kp, data = grange, ymin = 0,
                                             window.size = window.size, col = "blue", 
                                             r0 = r0, r1 = r1))

        kpAxis(kp, ymax = kp$latest.plot$computed.values$max.density, 
               cex = 2, r0 = r0, r1 = r1)
        kpAddLabels(kp, labels = granges_labels[track_no], r0 = r0, 
                    r1 = r1, label.margin = 0.07, cex = 2)
    }
}
                         
# Plot all the overlapping reads for all 'ids' for all 'srr's
plot_all_srrs <- function(srr_names, srrs_summary_table, 
                          crossings, donor_granges_all_srrs, recip_granges_all_srrs, 
                          recip_bams_all_srrs, donor_ref_genome, 
                          recipient_ref_genome, donor_name, recipient_name,
                          extend_left = 20, extend_right = 20) {

    # loop over all srrs
    for (srr_name in srr_names) {
        # delete all existing plots in the "plots" folder, and create directorty if it doesn't exist
        image_folder <- paste("plots", "/", donor_name, "/", srr_name, sep = "")
        suppressWarnings(dir.create(image_folder, recursive = TRUE))
        was_deleted <- do.call(file.remove, list(list.files(image_folder, full.names = TRUE)))
        
        # extract table for srr
        crossings_table_recipient <- srrs_summary_table[srrs_summary_table[, srr == srr_name], ]
        crossings_table_recipient <- crossings_table_recipient[, !"srr"]
        # extract bams for srr
        recip_bams <- unlist(recip_bams_all_srrs[[srr_name]])
        
        # extract id list
        ids <- unlist(as.list(crossings_table_recipient[,"id"]$id))
        # skip if there's no overlaps
        if (ids[[1]] != "<NA>" || length(ids[[1]]) == 0) {
             # display main title
            display_markdown(paste("###", srr_name))
            # display graph title
            display_markdown("#### Donor reads density graph")
            # graph donor analysis
            no_tracks <- length(donor_granges_all_srrs[[srr_name]])
            png(paste(image_folder, "/donor.png", sep = ""), 
                width = 1480, height = 200 + 100*no_tracks, res = 60)
            title_prepend <- paste(srr_name, ' aligned to ', 
                                   donor_name, ', and crossed with ', 
                                   recipient_name, sep = "")
            create_viz_donor(ref_genome = donor_ref_genome, 
                             granges = donor_granges_all_srrs[[srr_name]],  
                             granges_labels = crossings,
                             title_prepend = title_prepend)
            # display plot as image
            dev.off()
            display_png(file=paste(image_folder, "/donor.png", sep=""))

            # display graph title
            display_markdown("#### Recipient reads density graph")
            # graph recipient analysis
            no_tracks <- length(recip_granges_all_srrs[[srr_name]])
            png(paste(image_folder, "/recipient.png", sep = ""), 
                width = 1480, height = 200 + 100*no_tracks, res = 60)
            title_prepend <- paste(srr_name, ' aligned to ', 
                                   recipient_name, ', and crossed with ', 
                                   donor_name, sep = "")
            create_viz_recipient(ref_genome="hg38", 
                                 granges = recip_granges_all_srrs[[srr_name]], 
                                 granges_labels = crossings, 
                                 title_prepend = title_prepend)
            # display plot as image
            dev.off()
            display_png(file=paste(image_folder, "/recipient.png", sep=""))

            # create a visualization for all 'id's
            for (idn in ids) {
                    display_markdown("#### Crossings overlap graph - Putative insertion site")
                    # calculate total graph height, making it dependent of the number of crossings
                    crossings_number <- length(as.list(strsplit(
                        crossings_table_recipient[id == idn,]$unique_crossings[[1]], ",")))
                    height <- 400 + 150 * crossings_number
                    # open image writer
                    png(paste(image_folder, "/reads_", idn, ".png", sep = ""), 
                        width = 1480, height = height, res = 60)
                        # plot graph
                        plot_reads_region(srr = srr_name,
                                          id = as.integer(idn), 
                                          crossings_table_recipient = crossings_table_recipient, 
                                          recip_bams = recip_bams,
                                          extend_left = extend_left, 
                                          extend_right = extend_right, 
                                          ref_genome = recipient_ref_genome, 
                                          donor_name, recipient_name)
                    # display plot as image
                    dev.off()
                    display_png(file=paste(image_folder, "/reads_", idn,".png", sep=""))
            }
        }
    }
}