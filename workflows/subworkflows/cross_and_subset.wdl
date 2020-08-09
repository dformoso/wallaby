version development

import "../tasks/samtools.wdl" as samtools
import "../tasks/tools.wdl" as tools
import "../tasks/picard.wdl" as picard
import "../tasks/structs/structures.wdl"

workflow main {

    input {
        SplitBAMs donor_bams
        SplitQNAMEs donor_qnames
        SplitBAMs recipient_bams
        SplitQNAMEs recipient_qnames
    }

    # Performs inner joins between donor (left), and recipient (right), for the following pairs:     
    # MM_MM, MM_MU, MM_UM, MM_UU, 
    # MU_MM, MU_MU, MU_UM, MU_UU, 
    # UM_MM, UM_MU, UM_UM, UM_UU, 
    # UU_MM, UU_MU, UU_UM, UU_UU

    call tools.inner_join as MMd_MMr_txt { input: file_1 = donor_qnames.MM, file_2 = recipient_qnames.MM, out_file = "MMd_MMr.txt" }
    call tools.inner_join as MMd_MUr_txt { input: file_1 = donor_qnames.MM, file_2 = recipient_qnames.MU, out_file = "MMd_MUr.txt" }
    call tools.inner_join as MMd_UMr_txt { input: file_1 = donor_qnames.MM, file_2 = recipient_qnames.UM, out_file = "MMd_UMr.txt" }
    call tools.inner_join as MMd_UUr_txt { input: file_1 = donor_qnames.MM, file_2 = recipient_qnames.UU, out_file = "MMd_UUr.txt" }
    call tools.inner_join as MUd_MMr_txt { input: file_1 = donor_qnames.MU, file_2 = recipient_qnames.MM, out_file = "MUd_MMr.txt" }
    call tools.inner_join as MUd_MUr_txt { input: file_1 = donor_qnames.MU, file_2 = recipient_qnames.MU, out_file = "MUd_MUr.txt" }
    call tools.inner_join as MUd_UMr_txt { input: file_1 = donor_qnames.MU, file_2 = recipient_qnames.UM, out_file = "MUd_UMr.txt" }
    call tools.inner_join as MUd_UUr_txt { input: file_1 = donor_qnames.MU, file_2 = recipient_qnames.UU, out_file = "MUd_UUr.txt" }
    call tools.inner_join as UMd_MMr_txt { input: file_1 = donor_qnames.UM, file_2 = recipient_qnames.MM, out_file = "UMd_MMr.txt" }
    call tools.inner_join as UMd_MUr_txt { input: file_1 = donor_qnames.UM, file_2 = recipient_qnames.MU, out_file = "UMd_MUr.txt" }
    call tools.inner_join as UMd_UMr_txt { input: file_1 = donor_qnames.UM, file_2 = recipient_qnames.UM, out_file = "UMd_UMr.txt" }
    call tools.inner_join as UMd_UUr_txt { input: file_1 = donor_qnames.UM, file_2 = recipient_qnames.UU, out_file = "UMd_UUr.txt" }
    call tools.inner_join as UUd_MMr_txt { input: file_1 = donor_qnames.UU, file_2 = recipient_qnames.MM, out_file = "UUd_MMr.txt" }
    call tools.inner_join as UUd_MUr_txt { input: file_1 = donor_qnames.UU, file_2 = recipient_qnames.MU, out_file = "UUd_MUr.txt" }
    call tools.inner_join as UUd_UMr_txt { input: file_1 = donor_qnames.UU, file_2 = recipient_qnames.UM, out_file = "UUd_UMr.txt" }
    call tools.inner_join as UUd_UUr_txt { input: file_1 = donor_qnames.UU, file_2 = recipient_qnames.UU, out_file = "UUd_UUr.txt" }

    # Subset BAM files on the inner joins between donor/recipient
    # Donor
    call picard.filter_reads as donor_MMd_MMr { input: file = donor_bams.MM, txt = MMd_MMr_txt.out, out_file = "donor_MMd_MMr.bam" }
    call picard.filter_reads as donor_MMd_MUr { input: file = donor_bams.MM, txt = MMd_MUr_txt.out, out_file = "donor_MMd_MUr.bam" }
    call picard.filter_reads as donor_MMd_UMr { input: file = donor_bams.MM, txt = MMd_UMr_txt.out, out_file = "donor_MMd_UMr.bam" }
    call picard.filter_reads as donor_MMd_UUr { input: file = donor_bams.MM, txt = MMd_UUr_txt.out, out_file = "donor_MMd_UUr.bam" }
    call picard.filter_reads as donor_MUd_MMr { input: file = donor_bams.MU, txt = MUd_MMr_txt.out, out_file = "donor_MUd_MMr.bam" }
    call picard.filter_reads as donor_MUd_MUr { input: file = donor_bams.MU, txt = MUd_MUr_txt.out, out_file = "donor_MUd_MUr.bam" }
    call picard.filter_reads as donor_MUd_UMr { input: file = donor_bams.MU, txt = MUd_UMr_txt.out, out_file = "donor_MUd_UMr.bam" }
    call picard.filter_reads as donor_MUd_UUr { input: file = donor_bams.MU, txt = MUd_UUr_txt.out, out_file = "donor_MUd_UUr.bam" }
    call picard.filter_reads as donor_UMd_MMr { input: file = donor_bams.UM, txt = UMd_MMr_txt.out, out_file = "donor_UMd_MMr.bam" }
    call picard.filter_reads as donor_UMd_MUr { input: file = donor_bams.UM, txt = UMd_MUr_txt.out, out_file = "donor_UMd_MUr.bam" }
    call picard.filter_reads as donor_UMd_UMr { input: file = donor_bams.UM, txt = UMd_UMr_txt.out, out_file = "donor_UMd_UMr.bam" }
    call picard.filter_reads as donor_UMd_UUr { input: file = donor_bams.UM, txt = UMd_UUr_txt.out, out_file = "donor_UMd_UUr.bam" }
    call picard.filter_reads as donor_UUd_MMr { input: file = donor_bams.UU, txt = UUd_MMr_txt.out, out_file = "donor_UUd_MMr.bam" }
    call picard.filter_reads as donor_UUd_MUr { input: file = donor_bams.UU, txt = UUd_MUr_txt.out, out_file = "donor_UUd_MUr.bam" }
    call picard.filter_reads as donor_UUd_UMr { input: file = donor_bams.UU, txt = UUd_UMr_txt.out, out_file = "donor_UUd_UMr.bam" }
    call picard.filter_reads as donor_UUd_UUr { input: file = donor_bams.UU, txt = UUd_UUr_txt.out, out_file = "donor_UUd_UUr.bam" }

    # Recipient
    call picard.filter_reads as recipient_MMd_MMr { input: file = recipient_bams.MM, txt = MMd_MMr_txt.out, out_file = "recipient_MMd_MMr.bam" }
    call picard.filter_reads as recipient_MUd_MMr { input: file = recipient_bams.MM, txt = MUd_MMr_txt.out, out_file = "recipient_MUd_MMr.bam" }
    call picard.filter_reads as recipient_UMd_MMr { input: file = recipient_bams.MM, txt = UMd_MMr_txt.out, out_file = "recipient_UMd_MMr.bam" }
    call picard.filter_reads as recipient_UUd_MMr { input: file = recipient_bams.MM, txt = UUd_MMr_txt.out, out_file = "recipient_UUd_MMr.bam" }
    call picard.filter_reads as recipient_MMd_MUr { input: file = recipient_bams.MU, txt = MMd_MUr_txt.out, out_file = "recipient_MMd_MUr.bam" }
    call picard.filter_reads as recipient_MUd_MUr { input: file = recipient_bams.MU, txt = MUd_MUr_txt.out, out_file = "recipient_MUd_MUr.bam" }
    call picard.filter_reads as recipient_UMd_MUr { input: file = recipient_bams.MU, txt = UMd_MUr_txt.out, out_file = "recipient_UMd_MUr.bam" }
    call picard.filter_reads as recipient_UUd_MUr { input: file = recipient_bams.MU, txt = UUd_MUr_txt.out, out_file = "recipient_UUd_MUr.bam" }
    call picard.filter_reads as recipient_MMd_UMr { input: file = recipient_bams.UM, txt = MMd_UMr_txt.out, out_file = "recipient_MMd_UMr.bam" }
    call picard.filter_reads as recipient_MUd_UMr { input: file = recipient_bams.UM, txt = MUd_UMr_txt.out, out_file = "recipient_MUd_UMr.bam" }
    call picard.filter_reads as recipient_UMd_UMr { input: file = recipient_bams.UM, txt = UMd_UMr_txt.out, out_file = "recipient_UMd_UMr.bam" }
    call picard.filter_reads as recipient_UUd_UMr { input: file = recipient_bams.UM, txt = UUd_UMr_txt.out, out_file = "recipient_UUd_UMr.bam" }
    call picard.filter_reads as recipient_MMd_UUr { input: file = recipient_bams.UU, txt = MMd_UUr_txt.out, out_file = "recipient_MMd_UUr.bam" }
    call picard.filter_reads as recipient_MUd_UUr { input: file = recipient_bams.UU, txt = MUd_UUr_txt.out, out_file = "recipient_MUd_UUr.bam" }
    call picard.filter_reads as recipient_UMd_UUr { input: file = recipient_bams.UU, txt = UMd_UUr_txt.out, out_file = "recipient_UMd_UUr.bam" }
    call picard.filter_reads as recipient_UUd_UUr { input: file = recipient_bams.UU, txt = UUd_UUr_txt.out, out_file = "recipient_UUd_UUr.bam" }

    # Create a statistics, flag statistics, and count, file for each BAM file    
    Array[File?] all_bams = [
        donor_MMd_MMr.out, donor_MMd_MUr.out, donor_MMd_UMr.out, donor_MMd_UUr.out, 
        donor_MUd_MMr.out, donor_MUd_MUr.out, donor_MUd_UMr.out, donor_MUd_UUr.out, 
        donor_UMd_MMr.out, donor_UMd_MUr.out, donor_UMd_UMr.out, donor_UMd_UUr.out, 
        donor_UUd_MMr.out, donor_UUd_MUr.out, donor_UUd_UMr.out, donor_UUd_UUr.out, 
        recipient_MMd_MMr.out, recipient_MUd_MMr.out, recipient_UMd_MMr.out, recipient_UUd_MMr.out, 
        recipient_MMd_MUr.out, recipient_MUd_MUr.out, recipient_UMd_MUr.out, recipient_UUd_MUr.out, 
        recipient_MMd_UMr.out, recipient_MUd_UMr.out, recipient_UMd_UMr.out, recipient_UUd_UMr.out, 
        recipient_MMd_UUr.out, recipient_MUd_UUr.out, recipient_UMd_UUr.out, recipient_UUd_UUr.out]

    scatter(bam in select_all(all_bams)) {
        call samtools.stats { input: file = bam, out_file = "~{basename(bam)}_stats.txt" }
        call samtools.flagstats { input: file = bam, out_file = "~{basename(bam)}_flagstats.txt" }
        call samtools.count { input: file = bam, out_file = "~{basename(bam)}_count.txt" }
    }
    
    scatter (bam_name in select_all(all_bams)) { 
        String bam_names = "~{basename(bam_name)}" 
    }

    output {
        Array[Pair[String,File]] bams = zip(bam_names, select_all(all_bams))
    }

}