version development

import "../tasks/samtools.wdl" as samtools
import "../tasks/structs/compute.wdl"

workflow main {

    input {
        Array[File] bams
        File ref_genome
        Resources resources
    }

    scatter (bam in bams) {

        call samtools.mpileup as mpileup { 
            input: 
                fasta = ref_genome, 
                bam = bam,
                out_file = "~{basename(bam)}_~{basename(ref_genome)}.mpileup",
                resources = resources
        }
    }

    output {
        Array[File] mpileups = mpileup.out
    }
}