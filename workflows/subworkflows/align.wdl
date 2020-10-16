version 1.0

import "../tasks/bwa.wdl" as bwa
import "../tasks/samtools.wdl" as samtools
import "../tasks/structs/bwa.wdl" as bwa_struct
import "../tasks/structs/compute.wdl"

workflow main {

    input {
        File ref_genome
        File fastq_1
        File fastq_2
        String base_filename
        Int ignore_matches_shorted_than
        Int ignore_gaps_longer_than
        Int discard_if_repeated_in_ref_genome_more_than
        Int matching_score
        Int mismatch_penalty
        Int gap_open_penalty
        Int gap_extension_penalty
        Boolean output_all_found_alignments
        Resources resources
    }

    # Aligning FASTQ files to Reference Genomes
    call bwa.index as bwa_indexer { 
        input: 
            fasta = ref_genome,
            resources = resources

    }

    call bwa.align as bwa_aligner { 
        input: 
            bwa_index = bwa_indexer.index_object, 
            fastq_1 = fastq_1, 
            fastq_2 = fastq_2, 
            out_file = "~{base_filename}.sam",
            ignore_matches_shorted_than = ignore_matches_shorted_than,
            ignore_gaps_longer_than = ignore_gaps_longer_than,
            discard_if_repeated_in_ref_genome_more_than = discard_if_repeated_in_ref_genome_more_than,
            matching_score = matching_score,
            mismatch_penalty = mismatch_penalty,
            gap_open_penalty = gap_open_penalty,
            gap_extension_penalty = gap_extension_penalty,
            output_all_found_alignments = output_all_found_alignments,
            resources = resources
    }

    # Creating a BAM indexed file
    call samtools.view as samtools_sam_to_bam { 
        input: 
            file = bwa_aligner.out, 
            out_file = "~{base_filename}.bam",
            resources = resources
    }
    
    call samtools.sort as samtools_sorter { 
        input: 
            file = samtools_sam_to_bam.out, 
            out_file = "~{base_filename}.sorted.bam",
            resources = resources
    }

    call samtools.index as samtools_indexer { 
        input: 
            file = samtools_sorter.out, 
            out_file = "~{base_filename}.sorted.bam.bai",
            resources = resources
    }

    output {
        BWAIndex index = bwa_indexer.index_object
        File bam = samtools_sorter.out
        File bai = samtools_indexer.out
    }

}