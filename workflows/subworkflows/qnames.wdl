version 1.0

import "../tasks/samtools.wdl" as samtools

workflow main {

    input {
        File MM_bam
        File MU_bam
        File UM_bam
        File UU_bam
        String base_filename

        Resources resources
    }
        
    call samtools.extract_qnames as MM_qnames { 
        input: 
            file = MM_bam, 
            out_file = "~{base_filename}_MM_qnames.txt",
            resources = resources
    }

    call samtools.extract_qnames as MU_qnames { 
        input: 
            file = MU_bam, 
            out_file = "~{base_filename}_MU_qnames.txt",
            resources = resources
    }

    call samtools.extract_qnames as UM_qnames { 
        input: 
            file = UM_bam, 
            out_file = "~{base_filename}_UM_qnames.txt",
            resources = resources
    }

    call samtools.extract_qnames as UU_qnames { 
        input: 
            file = UU_bam, 
            out_file = "~{base_filename}_UU_qnames.txt",
            resources = resources
    }

    output {
        File MM = MM_qnames.out
        File MU = MU_qnames.out
        File UM = UM_qnames.out
        File UU = UU_qnames.out
    }

}