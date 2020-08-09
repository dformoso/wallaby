version development

task matching {
    input {
        File file
        String out_file
        String lc_method = "dust"
        String lc_threshold = "7"
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
        cpu: "24"
        memory: "16GB"
        docker: "dformoso/prinseq:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}

