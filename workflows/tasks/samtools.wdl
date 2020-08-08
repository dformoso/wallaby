##################################
########     SAMTOOLS     ########
######## TASK DEFINITIONS ########
##################################
version development

task view {
    input {
        File file
        String out_file
        
        Boolean include_header = false
        Boolean output_count = false
        Boolean output_file = true
        String include = ""
        String exclude = ""
        Int threads = 24
    }

    command <<<
        samtools view \
            ~{"-f " + include} \
            ~{"-F " + exclude} \
            ~{true="-h" false="" include_header} \
            ~{true="-b" false="" output_file} \
            ~{true="-c" false="" output_count} \
            ~{"-@ " + threads} \
            ~{file} \
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


task sort {
    input {
        File file
        String out_file
        
        Int threads = 24
    }

    command <<<
        samtools sort \
            ~{"-@ " + threads} \
            ~{file} \
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


task index {
    input {
        File file
        String out_file
        
        Int threads = 24
    }

    command <<<
        samtools index \
            ~{"-@ " + threads} \
            ~{file} \
            ~{out_file}
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


task merge {
    input {
        Array[File] files
        String out_file

        Int threads = 24
    }

    command <<<
        samtools merge \
            ~{"-@ " + threads} \
            ~{out_file} \
            ~{sep=" " files}
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


task stats {
    input {
        File file
        String out_file

        Int threads = 24
    }

    command <<<
        set -e
        if [ -s "~{file}" ]
        then
        samtools stats \
            ~{"-@ " + threads} \
            ~{file} \
            > ~{out_file}
        fi
    >>>

    output {
        File? stats = out_file
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

task flagstats {
    input {
        File file
        String out_file

        Int threads = 24
    }

    command <<<
        samtools flagstat \
            ~{"-@ " + threads} \
            ~{file} \
            > ~{out_file}
    >>>

    output {
        File? flagstats = out_file
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

# Create files containing a unique, ordered, list of QNAMES (the identifiers for the reads, or Read IDs)
task extract_qnames {
    input {
        File file
        String out_file

        Int threads = 24
    }

    command <<<
        samtools view \
            ~{"-@ " + threads} \
            ~{file} \
            | cut -f1 | sort | uniq \
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

task count {
    input {
        File file
        String out_file
        
        Boolean include_header = false
        Boolean output_count = true
        Boolean output_file = false
        String include = ""
        String exclude = ""
        Int threads = 24
    }

    command <<<
        set -e
        if [ -s "~{file}" ]
        then
        samtools view \
            ~{"-f " + include} \
            ~{"-F " + exclude} \
            ~{true="-h" false="" include_header} \
            ~{true="-b" false="" output_file} \
            ~{true="-c" false="" output_count} \
            ~{"-@ " + threads} \
            ~{file} \
            > ~{out_file}
        fi
    >>>

    output {
        File? bam = out_file
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

# Extract list of reads from BAM file
# Output paired reads to separate files, 
# discarding singletons, supplementary and secondary reads
task bam_to_fastas {
    input {
        File file
        String out_file_1
        String out_file_2

    }

    command <<<
        samtools fasta \
            ~{"-1 " + out_file_1} \
            ~{"-2 " + out_file_2} \
            ~{"-0 /dev/null"} \
            ~{"-s /dev/null"} \
            ~{"-n " + file}
    >>>

    output {
        File out_1 = out_file_1
        File out_2 = out_file_2
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