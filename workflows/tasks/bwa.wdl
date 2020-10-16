version 1.0

import "structs/bwa.wdl" as bwa_structs
import "structs/compute.wdl"

task index {
    
    input {
        File fasta 
        String basename_fasta = basename(fasta) 
        Resources resources
    }

    command <<<
        ln ~{fasta} .
        bwa index ~{basename_fasta} 
    >>>

    output {
        BWAIndex index_object = {
            "fasta": "~{basename_fasta}",
            "amb": "~{basename_fasta}.amb",
            "ann": "~{basename_fasta}.ann",
            "bwt": "~{basename_fasta}.bwt",
            "pac": "~{basename_fasta}.pac",
            "sa": "~{basename_fasta}.sa"
        }
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/bwa:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task align {

    input {
        File fastq_1
        File fastq_2
        BWAIndex bwa_index
        String out_file = "reads-to-ref-genome.sam"
        Int ignore_matches_shorted_than = 19
        Int ignore_gaps_longer_than = 100
        Int discard_if_repeated_in_ref_genome_more_than = 10000
        Int matching_score = 1
        Int mismatch_penalty = 4
        Int gap_open_penalty = 6
        Int gap_extension_penalty = 1
        Boolean output_all_found_alignments = true
        Resources resources
    }

    command <<<
        bwa mem \
            ~{"-k " + ignore_matches_shorted_than} \
            ~{"-w " + ignore_gaps_longer_than} \
            ~{"-c " + discard_if_repeated_in_ref_genome_more_than} \
            ~{"-A " + matching_score} \
            ~{"-B " + mismatch_penalty} \
            ~{"-O " + gap_open_penalty} \
            ~{"-E " + gap_extension_penalty} \
            ~{"-t " + resources.cpu} \
            ~{true="-a" false="" output_all_found_alignments} \
            ~{bwa_index.fasta} \
            ~{fastq_1} \
            ~{fastq_2} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/bwa:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}
