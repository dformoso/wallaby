version development

import "../tasks/blast.wdl" as blast
import "../tasks/structs/compute.wdl"

workflow main {

    input {
        Array[File] fastas
        Directory blastdb
        String evalue
        Resources resources
    }

    call blast.n as blaster {
        input: 
            fastas = fastas,
            blastdb = blastdb,
            resources = resources,
            evalue = evalue
    }
    
    output {
        Array[File] blastns = blaster.out
    }
}