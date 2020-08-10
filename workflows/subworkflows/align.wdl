version development

import "../tasks/bwa.wdl" as bwa
import "../tasks/samtools.wdl" as samtools
import "../tasks/structs/structures.wdl"

workflow main {

    input {
        File ref_genome
        File fastq_1
        File fastq_2
        String base_filename
    }

    # Aligning FASTQ files to Reference Genomes
    call bwa.index as bwa_indexer { 
        input: 
            fasta = ref_genome 
    }

    call bwa.align as bwa_aligner { 
        input: 
            bwa_index = bwa_indexer.index, 
            fastq_1 = fastq_1, 
            fastq_2 = fastq_2, 
            out_file = "~{base_filename}.sam"
    }

    # Creating a BAM indexed file
    call samtools.view as samtools_sam_to_bam { 
        input: 
            file = bwa_aligner.out, 
            out_file = "~{base_filename}.bam" 
    }
    
    call samtools.sort as samtools_sorter { 
        input: 
            file = samtools_sam_to_bam.out, 
            out_file = "~{base_filename}.sorted.bam" 
    }

    call samtools.index as samtools_indexer { 
        input: 
            file = samtools_sorter.out, 
            out_file = "~{base_filename}.sorted.bam.bai" 
    }

    output {
        BWAIndex index = bwa_indexer.index
        File bam = samtools_sorter.out
    }

}
