version 1.0

import "structs/compute.wdl"

task index {
    
    input {
        File fasta 
        String aligner_type
        String basename_fasta = basename(fasta) 
        Resources resources
    }

    command <<<
        if [ "~{aligner_type}" = 'bwa' ]; then
            ln ~{fasta} .
            bwa index ~{basename_fasta}
        elif [ "~{aligner_type}" = 'star' ]; then
            STAR \
                ~{"--runThreadN " + resources.cpu} \
                --runMode genomeGenerate \
                ~{"--genomeDir ./"} \
                ~{"--limitGenomeGenerateRAM 40000000000"} \
                ~{"--genomeFastaFiles " + fasta}
        fi
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
        String aligner_type
        Array[File] index_object
        String out_file = "reads-to-ref-genome.sam"

        Int bwa_ignore_matches_shorted_than = 19
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
        if [ "~{aligner_type}" = "bwa" ]; then

            for object in ~{sep=" " index_object}
            do
                ln $object .
            done

            bwa mem \
                ~{"-k " + bwa_ignore_matches_shorted_than} \
                ~{"-w " + bwa_ignore_gaps_longer_than} \
                ~{"-c " + bwa_discard_if_repeated_in_ref_genome_more_than} \
                ~{"-A " + bwa_matching_score} \
                ~{"-B " + bwa_mismatch_penalty} \
                ~{"-O " + bwa_gap_open_penalty} \
                ~{"-E " + bwa_gap_extension_penalty} \
                ~{"-t " + resources.cpu} \
                ~{true="-a" false="" bwa_output_all_found_alignments} \
                "$(ls -- *.fasta | head -n 1)" \
                ~{fastq_1} \
                ~{fastq_2} \
                > ~{out_file}

        elif [ "~{aligner_type}" = "star" ]; then

            mkdir index_dir
            for object in ~{sep=" " index_object}
            do
                ln $object index_dir/
            done

            STAR \
                ~{"--genomeDir index_dir"} \
                ~{"--runThreadN " + resources.cpu} \
                ~{"--limitGenomeGenerateRAM 40000000000"} \
                --readFilesIn ~{fastq_1} ~{fastq_2} \
                --outStd SAM \
                --outSAMtype BAM Unsorted \
                --outSAMunmapped Within KeepPairs \
                --outSAMattributes All \
                --outFilterMultimapNmax 100 \
                --winAnchorMultimapNmax 100 \
                --chimOutType WithinBAM SoftClip

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
