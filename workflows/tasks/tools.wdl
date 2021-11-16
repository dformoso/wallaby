version 1.0

import "structs/compute.wdl"

task inner_join {
    
    input {
        File file_1
        File file_2
        String out_file
        Resources resources
    }

    command <<<
        comm -12 \
            ~{file_1} \
            ~{file_2} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task seq_ids_from_fastq {
    
    input {
        File file
        String out_file
        Resources resources
    }

    command <<<
        cat \
            ~{file} \
            | grep @ | sed 's/^@//' | sed 's/\/.$//' \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task concat_text {
    
    input {
        File file_1
        File file_2
        String out_file
        Resources resources
    }

    command <<<
        cat \
            ~{file_1} \
            ~{file_2} \
            | sort ~{"--parallel=" + resources.cpu} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task merge_csvs {
    
    input {
        Array[File] csvs
        Resources resources
    }

    command <<<
        awk '(NR == 1) || (FNR > 1)' ~{sep=' ' csvs} > putative_insertion_table.csv
    >>>

    output {
        File out = "putative_insertion_table.csv"
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

task summary_and_inputs {
    
    input {
        String donor_name
        String recipient_name
        File donor_ref_genome
        File donor_ref_genome_fai
        File donor_ref_genome_gff
        Resources resources
    }

    command <<<
        echo "~{donor_name}, ~{recipient_name}" > donor_and_recipient.csv
        ln ~{donor_ref_genome} .
        ln ~{donor_ref_genome_fai} .
        ln ~{donor_ref_genome_gff} .
    >>>

    output {
        File summary = "donor_and_recipient.csv"
        File out_donor_ref_genome = basename(donor_ref_genome)
        File out_donor_ref_genome_fai = basename(donor_ref_genome_fai)
        File out_donor_ref_genome_gff = basename(donor_ref_genome_gff)
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/samtools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}