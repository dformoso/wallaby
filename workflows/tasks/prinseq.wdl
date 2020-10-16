version 1.0

import "structs/compute.wdl"

task matching {
    input {
        File file
        String out_file
        Int filter_shorter_than = 5
        Int filter_longer_than = 10000
        Int filter_if_gc_content_lower_than = 10
        Int filter_if_gc_content_higher_than = 90
        Int filter_if_avg_quality_below = 20
        String low_complexity_method = 'dust'
        String low_complexity_threshold = '7'
        Resources resources
    }

    command <<<
        prinseq \
            ~{"-fastq=" + file} \
            ~{"-out_good=stdout"} \
            ~{"-out_bad=null"} \
            ~{"-min_len="+ filter_shorter_than} \
            ~{"-max_len="+ filter_longer_than} \
            ~{"-min_gc="+ filter_if_gc_content_lower_than} \
            ~{"-max_gc="+ filter_if_gc_content_higher_than} \
            ~{"-min_qual_mean="+ filter_if_avg_quality_below} \
            ~{"-lc_method=" + low_complexity_method} \
            ~{"-lc_threshold=" + low_complexity_threshold} \
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

