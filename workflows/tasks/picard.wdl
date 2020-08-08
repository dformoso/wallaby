##################################
########      PICARD      ########
######## TASK DEFINITIONS ########
##################################
version development

task filter_reads {
    input {
        File file
        File txt
        String out_file

        String filter = "includeReadList"
    }

    command {
        set -e
        if [ -s "~{txt}" ]
        then
            if [ -s "~{file}" ] 
            then 
            java -jar /usr/picard/picard.jar FilterSamReads \
                ~{"INPUT=" + file} \
                ~{"OUTPUT=" + out_file} \
                ~{"READ_LIST_FILE=" + txt} \
                ~{"FILTER=" + filter}
            fi
        fi
    }

    output {
        File? out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: "24"
        memory: "16GB"
        docker: "dformoso/picard-tools:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}

task filter_valid_reads {
    input {
        File file
        File txt
        String out_file

        String filter = "includeReadList"
    }

    command {
        set -e
        if [ -s "~{txt}" ]
        then
            if [ -s "~{file}" ] 
            then 
            java -jar /usr/picard/picard.jar FilterSamReads \
                ~{"INPUT=" + file} \
                ~{"OUTPUT=" + out_file} \
                ~{"READ_LIST_FILE=" + txt} \
                ~{"FILTER=" + filter}
            fi
        fi
    }

    output {
        File out = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: "24"
        memory: "16GB"
        docker: "dformoso/picard-tools:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}

task mark_duplicates {
    input {
        File file
        String out_file

        String validation_stringency = "SILENT"
        String assume_sorted = "false"
        String minimum_distance = "-1"
        String validation_stringency = "SILENT"
    }

    command {
        set -e
        java -jar /usr/picard/picard.jar MarkDuplicates \
            ~{"INPUT=" + file} \
            ~{"OUTPUT=" + out_file} \
            ~{"ASSUME_SORTED=" + assume_sorted} \
            ~{"MINIMUM_DISTANCE=" + minimum_distance} \
            ~{"VALIDATION_STRINGENCY=" + validation_stringency}
    }

    output {
        File? deduped = out_file
    }

    runtime {
        continueOnReturnCode: false
        cpu: "24"
        memory: "16GB"
        docker: "dformoso/picard-tools:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}


task insert_size_metrics {
    input {
        File file
        String out_file
        String pdf_filename

        String validation_stringency = "SILENT"
        String assume_sorted = "false"
    }

    command {
        set -e
        java -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
            ~{"INPUT=" + file} \
            ~{"OUTPUT=" + out_file} \
            ~{"HISTOGRAM_FILE=" + pdf_filename} \
            ~{"ASSUME_SORTED=" + assume_sorted} \
            ~{"VALIDATION_STRINGENCY=" + validation_stringency} 
    }

    output {
        File? insert_size_metrics = out_file
        File? pdf = pdf_filename
    }

    runtime {
        continueOnReturnCode: false
        cpu: "24"
        memory: "16GB"
        docker: "dformoso/picard-tools:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}

task sam_to_fastq {
    input {
        File file
        String out_fastq_1
        String out_fastq_2

        String validation_stringency = "SILENT"
    }

    command {
        set -e
        if [ -s "~{file}" ] 
        then 
        java -jar /usr/picard/picard.jar SamToFastq \
            ~{"INPUT=" + file} \
            ~{"FASTQ=" + out_fastq_1} \
            ~{"SECOND_END_FASTQ=" + out_fastq_2} \
            ~{"VALIDATION_STRINGENCY=" + validation_stringency} 
        fi
    }

    output {
        File fastq_1 = out_fastq_1
        File fastq_2 = out_fastq_2
    }

    runtime {
        continueOnReturnCode: false
        cpu: "24"
        memory: "16GB"
        docker: "dformoso/picard-tools:latest"
        disks: "local-disk 100GB HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 0
        zones: "us-central1-c"
    }
}