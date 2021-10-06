version 1.0

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
        Array[File] index_object = glob("*")
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/aligners:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task align {

    input {
        File fastq_1
        File fastq_2
        Array[File] index_object
        String out_file = "reads-to-ref-genome.sam"

        Int bwa_ignore_matches_shorter_than = 19
        Int bwa_ignore_gaps_longer_than = 100
        Int bwa_discard_if_repeated_in_ref_genome_more_than = 10000
        Int bwa_matching_score = 1
        Int bwa_mismatch_penalty = 4
        Int bwa_gap_open_penalty = 6
        Int bwa_gap_extension_penalty = 1
        Boolean bwa_output_all_found_alignments = true

        Resources resources
    }

    command <<<
        for object in ~{sep=" " index_object}
        do
            ln $object .
        done

        bwa mem \
            ~{"-k " + bwa_ignore_matches_shorter_than} \
            ~{"-w " + bwa_ignore_gaps_longer_than} \
            ~{"-c " + bwa_discard_if_repeated_in_ref_genome_more_than} \
            ~{"-A " + bwa_matching_score} \
            ~{"-B " + bwa_mismatch_penalty} \
            ~{"-O " + bwa_gap_open_penalty} \
            ~{"-E " + bwa_gap_extension_penalty} \
            ~{"-t " + resources.cpu} \
            ~{true="-a" false="" bwa_output_all_found_alignments} \
            "$(ls -- *.fa | head -n 1)" \
            ~{fastq_1} \
            ~{fastq_2} \
            > ~{out_file}

            mv Aligned.out.bam ~{out_file}
        fi
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/aligners:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task align_convert_index {

    input {
        Array[File] index_object
        String base_filename
        File fastq_1
        File fastq_2

        Int bwa_ignore_matches_shorter_than = 19
        Int bwa_ignore_gaps_longer_than = 100
        Int bwa_discard_if_repeated_in_ref_genome_more_than = 10000
        Int bwa_matching_score = 1
        Int bwa_mismatch_penalty = 4
        Int bwa_gap_open_penalty = 6
        Int bwa_gap_extension_penalty = 1
        Boolean bwa_output_all_found_alignments = true

        Boolean include_header = false
        Boolean output_count = false
        Boolean output_file = true
        String include = ""
        String exclude = ""

        Resources resources
    }

    command <<<
        for object in ~{sep=" " index_object}
            do
                ln $object .
            done

            bwa mem \
                ~{"-k " + bwa_ignore_matches_shorter_than} \
                ~{"-w " + bwa_ignore_gaps_longer_than} \
                ~{"-c " + bwa_discard_if_repeated_in_ref_genome_more_than} \
                ~{"-A " + bwa_matching_score} \
                ~{"-B " + bwa_mismatch_penalty} \
                ~{"-O " + bwa_gap_open_penalty} \
                ~{"-E " + bwa_gap_extension_penalty} \
                ~{"-t " + resources.cpu} \
                ~{true="-a" false="" bwa_output_all_found_alignments} \
                "$(ls -- *.fa | head -n 1)" \
                ~{fastq_1} \
                ~{fastq_2} \
                > "~{base_filename}.sam"

            samtools view \
                ~{"-f " + include} \
                ~{"-F " + exclude} \
                ~{true="-h" false="" include_header} \
                ~{true="-b" false="" output_file} \
                ~{true="-c" false="" output_count} \
                ~{"-@ " + resources.cpu} \
                "~{base_filename}.sam" \
                > "~{base_filename}_raw.bam"

            rm -rf "~{base_filename}.sam"

            samtools sort \
                ~{"-@ " + resources.cpu} \
                "~{base_filename}_raw.bam" \
                > "~{base_filename}.bam"

            rm -rf "~{base_filename}_raw.bam"

            samtools index \
                ~{"-@ " + resources.cpu} \
                "~{base_filename}.bam" \
                "~{base_filename}.bai"
    >>>

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/bwa_samtools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }

    output {
        File bam = "~{base_filename}.bam"
        File bai = "~{base_filename}.bai"
    }

}