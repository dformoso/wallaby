version development

import "../tasks/samtools.wdl" as samtools

workflow main {

    input {
        Array[Pair[String,File]] bams
        File ref_genome
    }

    scatter (pair in bams) {

        String filename = pair.left
        File file = pair.right

        call samtools.mpileup as mpileup { 
            input: 
                fasta = ref_genome, 
                bam = file,
                out_file = "~{filename}_~{basename(ref_genome)}.mpileup"
        }
    }

    output {
        Array[File] mpileups = mpileup.out
    }
}