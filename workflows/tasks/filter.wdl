version 1.0

import "structs/compute.wdl"

task bam_reads {
    input {
        File bam
        String filter = "true"
        Resources resources

        String validation_stringency = "SILENT"

        Int filter_shorter_than = 5
        Int filter_longer_than = 10000
        Int filter_if_gc_content_lower_than = 10
        Int filter_if_gc_content_higher_than = 90
        Int filter_if_avg_quality_below = 20
        String low_complexity_method = 'dust'
        String low_complexity_threshold = '7'

        String filter_type = "includeReadList"
    }

    command <<<

        if [ "~{filter}" = "true" ]; then

            # Convert Bam/Sam to FastQ
            set -e
            if [ -s "~{bam}" ] 
            then 
            java -Xmx40G -jar /usr/picard/picard.jar SamToFastq \
                ~{"INPUT=" + bam} \
                ~{"FASTQ=" + "raw.1.fastq"} \
                ~{"SECOND_END_FASTQ=" + "raw.2.fastq"} \
                ~{"VALIDATION_STRINGENCY=" + validation_stringency} 
            fi

            # Filter FASTQ files to match given parameters
            prinseq \
                ~{"-fastq=" + "raw.1.fastq"} \
                ~{"-out_good=stdout"} \
                ~{"-out_bad=null"} \
                ~{"-min_len="+ filter_shorter_than} \
                ~{"-max_len="+ filter_longer_than} \
                ~{"-min_gc="+ filter_if_gc_content_lower_than} \
                ~{"-max_gc="+ filter_if_gc_content_higher_than} \
                ~{"-min_qual_mean="+ filter_if_avg_quality_below} \
                ~{"-lc_method=" + low_complexity_method} \
                ~{"-lc_threshold=" + low_complexity_threshold} \
                > "filtered_reads.1.fastq"

            prinseq \
                ~{"-fastq=" + "raw.2.fastq"} \
                ~{"-out_good=stdout"} \
                ~{"-out_bad=null"} \
                ~{"-min_len="+ filter_shorter_than} \
                ~{"-max_len="+ filter_longer_than} \
                ~{"-min_gc="+ filter_if_gc_content_lower_than} \
                ~{"-max_gc="+ filter_if_gc_content_higher_than} \
                ~{"-min_qual_mean="+ filter_if_avg_quality_below} \
                ~{"-lc_method=" + low_complexity_method} \
                ~{"-lc_threshold=" + low_complexity_threshold} \
                > "filtered_reads.2.fastq"

            # Extract the SEQ_IDs from the FASTQ files into .txt files
            cat \
                "filtered_reads.1.fastq" \
                | grep @ | sed 's/^@//' | sed 's/\/.$//' \
                > "filtered_reads.1.txt"

            cat \
                "filtered_reads.2.fastq" \
                | grep @ | sed 's/^@//' | sed 's/\/.$//' \
                > "filtered_reads.2.txt"

            # Concatanate the text files for fastq1 and fastq2 into one
            cat \
                "filtered_reads.1.txt" \
                "filtered_reads.2.txt" \
                | sort ~{"--parallel=" + resources.cpu} \
                > "filtered_reads.merged.txt"

            # Filter the BAM file so that they only have our selected SEQ_IDs
            if [ -s "filtered_reads.merged.txt" ]
            then
                if [ -s "~{bam}" ] 
                then 
                java -jar /usr/picard/picard.jar FilterSamReads \
                    ~{"INPUT=" + bam} \
                    ~{"OUTPUT=" + "~{basename(bam, ".bam")}_filtered.bam"} \
                    ~{"READ_LIST_FILE=" + "filtered_reads.merged.txt"} \
                    ~{"FILTER=" + filter_type}
                fi
            fi

            # Remove all intermediate files
            #rf -rf "filtered_reads.1.fastq"
            #rf -rf "filtered_reads.2.fastq"
            #rm -rf "filtered_reads.1.txt"
            #rm -rf "filtered_reads.2.txt"
            #rm -rf "filtered_reads.merged.txt"
        fi
    >>>

    output {
        File bams = "~{basename(bam, ".bam")}_filtered.bam"
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/filter-tools:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}