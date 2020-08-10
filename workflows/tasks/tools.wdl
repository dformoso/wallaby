version development

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
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
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
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
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
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}