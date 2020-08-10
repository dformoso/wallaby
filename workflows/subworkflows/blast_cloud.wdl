version development

import "../tasks/blast_cloud.wdl" as blast

workflow main {

    input {
        Array[File] fastas
        Directory blastdb
    }

    # Select fasta files of interest for analysis
    scatter (fasta in fastas) {

        if (
            basename(fasta, ".bam.1.fasta")=="donor_MMd_MUr" || 
            basename(fasta, ".bam.1.fasta")=="donor_MMd_UMr" || 
            basename(fasta, ".bam.1.fasta")=="donor_MUd_MUr" || 
            basename(fasta, ".bam.1.fasta")=="donor_MUd_UMr" || 
            basename(fasta, ".bam.1.fasta")=="donor_UMd_MUr" || 
            basename(fasta, ".bam.1.fasta")=="donor_UMd_UMr" || 
            basename(fasta, ".bam.1.fasta")=="recipient_MMd_MUr" || 
            basename(fasta, ".bam.1.fasta")=="recipient_MMd_UMr" || 
            basename(fasta, ".bam.1.fasta")=="recipient_MUd_MUr" || 
            basename(fasta, ".bam.1.fasta")=="recipient_MUd_UMr" || 
            basename(fasta, ".bam.1.fasta")=="recipient_UMd_MUr" || 
            basename(fasta, ".bam.1.fasta")=="recipient_UMd_UMr" || 
            basename(fasta, ".bam.2.fasta")=="donor_MMd_MUr" || 
            basename(fasta, ".bam.2.fasta")=="donor_MMd_UMr" || 
            basename(fasta, ".bam.2.fasta")=="donor_MUd_MUr" || 
            basename(fasta, ".bam.2.fasta")=="donor_MUd_UMr" || 
            basename(fasta, ".bam.2.fasta")=="donor_UMd_MUr" || 
            basename(fasta, ".bam.2.fasta")=="donor_UMd_UMr" || 
            basename(fasta, ".bam.2.fasta")=="recipient_MMd_MUr" || 
            basename(fasta, ".bam.2.fasta")=="recipient_MMd_UMr" || 
            basename(fasta, ".bam.2.fasta")=="recipient_MUd_MUr" || 
            basename(fasta, ".bam.2.fasta")=="recipient_MUd_UMr" || 
            basename(fasta, ".bam.2.fasta")=="recipient_UMd_MUr" || 
            basename(fasta, ".bam.2.fasta")=="recipient_UMd_UMr"
            ) {

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

