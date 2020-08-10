version development

import "../tasks/blast.wdl" as blast

workflow main {

    input {
        Array[File] fastas
        Directory blastdb
    }

    call blast.n as blaster {
        input: 
            fastas = fastas,
            blastdb = blastdb
    }
    
    output {
        Array[File] blastns = blaster.out
    }
}