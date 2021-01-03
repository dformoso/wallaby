version 1.0

import "../tasks/blast.wdl" as blast
import "../tasks/structs/compute.wdl"

workflow main {

    input {
        Array[File] fastas
        File blastdb
        Int evalue
        Resources resources
    }

    call blast.n as blaster {
        input: 
            fastas = fastas,
            blastdb = blastdb,
            evalue = evalue,
            resources = resources
    }
    
    output {
        Array[File] blasted = blaster.blasted
    }
}