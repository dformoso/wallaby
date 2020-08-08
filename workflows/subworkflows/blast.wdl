#####################################
########        BLAST       #########
######## PIPELINE DEFINITION ########
#####################################

version development

import "../tasks/blast.wdl" as blast

workflow main {

    input {
        Array[File] fastas
        Directory blastdb
    }

    # Keeping only the following files:
    ## donor_MMd_MUr
    ## donor_MMd_UMr
    ## donor_MUd_MUr
    ## donor_MUd_UMr
    ## donor_UMd_MUr
    ## donor_UMd_UMr
    ## recipient_MMd_MUr
    ## recipient_MMd_UMr
    ## recipient_MUd_MUr
    ## recipient_MUd_UMr
    ## recipient_UMd_MUr
    ## recipient_UMd_UMr

    # Select fasta files of interest for analysis
    scatter (fasta in fastas) {

        if (basename(fasta, ".bam.1.fasta")=="donor_MMd_MUr" || basename(fasta, ".bam.1.fasta")=="donor_MMd_UMr" || basename(fasta, ".bam.1.fasta")=="donor_MUd_MUr" || basename(fasta, ".bam.1.fasta")=="donor_MUd_UMr" || basename(fasta, ".bam.1.fasta")=="donor_UMd_MUr" || basename(fasta, ".bam.1.fasta")=="donor_UMd_UMr" || basename(fasta, ".bam.1.fasta")=="recipient_MMd_MUr" || basename(fasta, ".bam.1.fasta")=="recipient_MMd_UMr" || basename(fasta, ".bam.1.fasta")=="recipient_MUd_MUr" || basename(fasta, ".bam.1.fasta")=="recipient_MUd_UMr" || basename(fasta, ".bam.1.fasta")=="recipient_UMd_MUr" || basename(fasta, ".bam.1.fasta")=="recipient_UMd_UMr" || basename(fasta, ".bam.2.fasta")=="donor_MMd_MUr" || basename(fasta, ".bam.2.fasta")=="donor_MMd_UMr" || basename(fasta, ".bam.2.fasta")=="donor_MUd_MUr" || basename(fasta, ".bam.2.fasta")=="donor_MUd_UMr" || basename(fasta, ".bam.2.fasta")=="donor_UMd_MUr" || basename(fasta, ".bam.2.fasta")=="donor_UMd_UMr" || basename(fasta, ".bam.2.fasta")=="recipient_MMd_MUr" || basename(fasta, ".bam.2.fasta")=="recipient_MMd_UMr" || basename(fasta, ".bam.2.fasta")=="recipient_MUd_MUr" || basename(fasta, ".bam.2.fasta")=="recipient_MUd_UMr" || basename(fasta, ".bam.2.fasta")=="recipient_UMd_MUr" || basename(fasta, ".bam.2.fasta")=="recipient_UMd_UMr") {
            call blast.n as blaster {
                input: 
                    fasta = fasta,
                    blastdb = blastdb
            }
            File blast_file = blaster.out
        }

    }

    output {
        Array[File] blastns = select_all(blast_file)
    }

}


