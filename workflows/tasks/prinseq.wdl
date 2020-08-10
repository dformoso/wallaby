version development

import "structs/compute.wdl"

task matching {
    input {
        File file
        String out_file
        String lc_method
        String lc_threshold
        Resources resources
    }

    command <<<
        prinseq \
            ~{"-fastq=" + file} \
            ~{"-out_good=stdout"} \
            ~{"-out_bad=null"} \
            ~{"-lc_method=" + lc_method} \
            ~{"-lc_threshold=" + lc_threshold} \
            > ~{out_file}
    >>>

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/prinseq:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}

