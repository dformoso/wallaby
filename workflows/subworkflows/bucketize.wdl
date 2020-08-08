#####################################
########       BUCKET        ########
######## PIPELINE DEFINITION ########
#####################################
version development

import "../tasks/bwa.wdl" as bwa
import "../tasks/samtools.wdl" as samtools
import "../tasks/structs/structures.wdl"

workflow main {

    input {
        BWAIndex index
        File bam
        File bai
    }

    String base_filename = "reads-to-~{basename(index.fasta, ".fasta")}"

    # Read 1 is Mapped and Read 2 is Mapped in/to Ref Genome - MM
    call samtools.view as MM { input: file = bam, out_file = "~{base_filename}_MM.bam", include = "3", exclude = "2048" }

    # Read 1 is Mapped and Read 2 is Unmapped in/to Ref Genome - MU
    call samtools.view as MU_R1 { input: file = bam, out_file = "~{base_filename}_MU_R1.bam", include = "73", exclude = "4022" }
    call samtools.view as MU_R2 { input: file = bam, out_file = "~{base_filename}_MU_R2.bam", include = "133", exclude = "3962" }
    call samtools.merge as MU { input: files = [MU_R1.out, MU_R2.out], out_file = "~{base_filename}_MU.bam" }

    # Read 1 is Unmapped and Read 2 is Mapped in/to Ref Genome - UM
    call samtools.view as UM_R1 { input: file = bam, out_file = "~{base_filename}_UM_R1.bam", include = "69", exclude = "4026" }
    call samtools.view as UM_R2 { input: file = bam, out_file = "~{base_filename}_UM_R2.bam", include = "137", exclude = "3958" }
    call samtools.merge as UM { input: files = [UM_R1.out, UM_R2.out], out_file = "~{base_filename}_UM.bam" }

    # Read 1 is Unmapped and Read 2 is Unmapped in/to Ref Genome - UU
    call samtools.view as UU_R1 { input: file = bam, out_file = "~{base_filename}_UU_R1.bam", include = "77" }
    call samtools.view as UU_R2 { input: file = bam, out_file = "~{base_filename}_UU_R2.bam", include = "141" }
    call samtools.merge as UU { input: files = [UU_R1.out, UU_R2.out], out_file = "~{base_filename}_UU.bam" }

    # Create files containing a unique and ordered, list of QNAMES (Read IDs)
    call samtools.extract_qnames as MM_qnames { input: file = MM.out, out_file = "~{base_filename}_MM_qnames.txt" }
    call samtools.extract_qnames as MU_qnames { input: file = MU.out, out_file = "~{base_filename}_MU_qnames.txt"}
    call samtools.extract_qnames as UM_qnames { input: file = UM.out, out_file = "~{base_filename}_UM_qnames.txt" }
    call samtools.extract_qnames as UU_qnames { input: file = UU.out, out_file = "~{base_filename}_UU_qnames.txt"}

    # Create a statistics, flag statistics, and count, file for each BAM file    
    Array[File] bam_files = [MM.out, MU.out, UM.out, UU.out]

    scatter (bam_file in bam_files) {
        call samtools.stats { input: file = bam_file, out_file = "~{basename(bam_file)}_stats.txt" }
        call samtools.flagstats { input: file = bam_file, out_file = "~{basename(bam_file)}_flagstats.txt" }
        call samtools.count as MM_bam_count { input: file = bam_file, out_file = "~{basename(bam_file)}_count.txt" }
    }

    output {
        SplitBAMs bams = {
            "MM" : MM.out,
            "MU" : MU.out,
            "UM" : UM.out,
            "UU" : UU.out
        }

        SplitQNAMEs qnames = {
            "MM" : MM_qnames.out,
            "MU" : MU_qnames.out,
            "UM" : UM_qnames.out,
            "UU" : UU_qnames.out
        }
    }

}


