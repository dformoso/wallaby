version development

import "structs/compute.wdl"

task trim {
    input {
        File fastq_1
        File fastq_2
        String out_fastq_1_paired = "~{basename(fastq_1)}.trimmed.fastq"
        String out_fastq_1_unpaired = "~{basename(fastq_1)}.unpaired.trimmed.fastq"
        String out_fastq_2_paired = "~{basename(fastq_2)}.trimmed.fastq"
        String out_fastq_2_unpaired = "~{basename(fastq_2)}.unpaired.trimmed.fastq"
        String out_trimlog = "out_trimlog.txt"
        String? adapter = "TruSeq3-PE.fa"
        Int? seed_mismatches = 2
        Int? paired_clip_threshold = 30
        Int? unpaired_clip_threshold = 10
        Int? leading = 3
        Int? trailing = 3
        Int? sliding_window_quality = 4
        Int? sliding_window_length = 15
        Int? min_length = 36
        Boolean? is_phred33 = true
        Resources resources
    }

    command <<<
        set -e
        trimmomatic \
            PE \
            ~{true="-phred33" false="-phred64" is_phred33} \
            ~{"-threads " + resources.cpu} \
            ~{"-trimlog " + out_trimlog} \
            ~{fastq_1} \
            ~{fastq_2} \
            ~{out_fastq_1_paired} \
            ~{out_fastq_1_unpaired} \
            ~{out_fastq_2_paired} \
            ~{out_fastq_2_unpaired} \
            ILLUMINACLIP:~{adapter}:~{seed_mismatches}:~{paired_clip_threshold}:~{unpaired_clip_threshold} \
            ~{"LEADING:" + leading} \
            ~{"TRAILING:" + trailing} \
            SLIDINGWINDOW:~{sliding_window_length}:~{sliding_window_quality} \
            ~{"MINLEN:" + min_length}
    >>>

    output {
        File fastq_1_paired = out_fastq_1_paired 
        File fastq_1_unpaired = out_fastq_1_unpaired 
        File fastq_2_paired = out_fastq_2_paired 
        File fastq_2_unpaired = out_fastq_2_unpaired 
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/trimmomatic:latest"
        disks: resources.disks
        gpuType: resources.gpuType
        gpuCount: resources.gpuCount
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}
