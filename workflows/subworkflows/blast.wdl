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
        File donor_MMd_MUr = blaster.donor_MMd_MUr
        File donor_MUd_UMr = blaster.donor_MUd_UMr
        File donor_UMd_MUr = blaster.donor_UMd_MUr
        File recipient_MMd_MUr = blaster.recipient_MMd_MUr
        File recipient_MUd_UMr = blaster.recipient_MUd_UMr
        File recipient_UMd_MUr = blaster.recipient_UMd_MUr
    }
}