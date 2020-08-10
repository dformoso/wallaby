version development

import "../tasks/samtools.wdl" as samtools
import "../tasks/structs/compute.wdl"

workflow main {

    input {
        Array[File] bams
        Resources resources
    }

    scatter (bam in bams) {
        
        call samtools.stats as samtools_stats { 
            input: 
                file = bam, 
                out_file = "~{basename(bam)}_stats.txt",
                resources = resources
        }

        call samtools.flagstats as samtools_flagstats { 
            input: 
                file = bam, 
                out_file = "~{basename(bam)}_flagstats.txt",
                resources = resources
        }

        call samtools.count as samtools_count { 
            input: 
                file = bam, 
                out_file = "~{basename(bam)}_count.txt",
                resources = resources
        }
    }

    output {
        Array[File] stats = samtools_stats.out
        Array[File] flagstats = samtools_flagstats.out
        Array[File] count = samtools_count.out
    }
}