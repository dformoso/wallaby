version development

task inner_join {
    
    input {
        File file_1
        File file_2
        String out_file
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
        cpu: "24"
        memory: "16GB"
        docker: "dformoso/samtools:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}

task seq_ids_from_fastq {
    
    input {
        File file
        String out_file
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
        cpu: "24"
        memory: "16GB"
        docker: "dformoso/samtools:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}

task concat_text {
    
    input {
        File file_1
        File file_2
        String out_file

        Int threads = 24
    }

    command <<<
        cat \
            ~{file_1} \
            ~{file_2} \
            | sort ~{"--parallel=" + threads} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: "24"
        memory: "16GB"
        docker: "dformoso/samtools:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}