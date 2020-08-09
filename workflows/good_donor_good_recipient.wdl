version development

import "subworkflows/align.wdl" as align
import "subworkflows/bucketize.wdl" as bucketize
import "subworkflows/cross_and_subset.wdl" as cross_and_subset
import "subworkflows/low_complex_filter.wdl" as low_complex_filter
import "subworkflows/blast.wdl" as blast
import "subworkflows/locate.wdl" as locate

workflow good_donor_good_recipient {

    input {
        Directory blastdb
        File donor_ref_genome
        File recipient_ref_genome
        File fastq_1
        File fastq_2
    }

    # Donor Reference Genome
    # Align 
    call align.main as donor_align {
        input:
            ref_genome = donor_ref_genome,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
    }

    # Bucketize
    call bucketize.main as donor_bucketize {
        input:
            index = donor_align.index,
            bam = donor_align.bam,
            bai = donor_align.bai
    }

    # Recipient Reference Genome
    # Align 
    call align.main as recipient_align {
        input:
            ref_genome = recipient_ref_genome,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
    }

    # Bucketize
    call bucketize.main as recipient_bucketize {
        input:
            index = recipient_align.index,
            bam = recipient_align.bam,
            bai = recipient_align.bai
    }

    # Cross BAM files buckets
    call cross_and_subset.main as cross_subset {
        input:
            donor_bams = donor_bucketize.bams,
            donor_qnames = donor_bucketize.qnames,
            recipient_bams = recipient_bucketize.bams,
            recipient_qnames = recipient_bucketize.qnames
    }

    # Filter out low complexity sequences
    call low_complex_filter.main as low_complex_filter {
        input:
            bams = cross_subset.bams
    }

    # Locate sequences in reference genome
    call locate.main as donor_locator {
        input:
            bams = cross_subset.bams,
            ref_genome = donor_ref_genome
    }

    call locate.main as recipient_locator {
        input:
            bams = cross_subset.bams,
            ref_genome = recipient_ref_genome
    }

    # Blastn search over all crossings of interest
    #call blast.main as blaster {
    #    input:
    #        fastas = low_complex_filter.fastas,
    #        blastdb = blastdb
    #}

    output {
        Array[File] fastas = low_complex_filter.fastas
        #Array[File] blastns = blaster.blastns
    }

}