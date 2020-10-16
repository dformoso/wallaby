version 1.0

import "qnames.wdl" as qnames
import "../tasks/tools.wdl" as tools
import "../tasks/picard.wdl" as picard
import "../tasks/structs/compute.wdl"

workflow main {

    input {
        File donor_MM
        File donor_MU
        File donor_UM
        File donor_UU
        File recipient_MM
        File recipient_MU
        File recipient_UM
        File recipient_UU

        Resources resources
    }

    ## Extract QNAMEs
    call qnames.main as donor_qnames {
        input:
            MM_bam = donor_MM,
            MU_bam = donor_MU,
            UM_bam = donor_UM,
            UU_bam = donor_UU,
            base_filename = "reads-to-donor",
            resources = resources
    }

    call qnames.main as recipient_qnames {
        input:
            MM_bam = recipient_MM,
            MU_bam = recipient_MU,
            UM_bam = recipient_UM,
            UU_bam = recipient_UU,
            base_filename = "reads-to-recipient",
            resources = resources
    }

    # Performs inner joins between donor (left), 
    # and recipient (right), for the following pairs:     
    # MM_MM, MM_MU, MM_UM, MM_UU, 
    # MU_MM, MU_MU, MU_UM, MU_UU, 
    # UM_MM, UM_MU, UM_UM, UM_UU, 
    # UU_MM, UU_MU, UU_UM, UU_UU

    call tools.inner_join as MMd_MMr_txt { 
        input: 
            file_1 = donor_qnames.MM, 
            file_2 = recipient_qnames.MM, 
            out_file = "MMd_MMr.txt",
            resources = resources
        }

    call tools.inner_join as MMd_MUr_txt { 
        input: 
            file_1 = donor_qnames.MM, 
            file_2 = recipient_qnames.MU, 
            out_file = "MMd_MUr.txt",
            resources = resources
        }

    call tools.inner_join as MMd_UMr_txt { 
        input: 
            file_1 = donor_qnames.MM, 
            file_2 = recipient_qnames.UM, 
            out_file = "MMd_UMr.txt",
            resources = resources
        }

    call tools.inner_join as MMd_UUr_txt { 
        input: 
            file_1 = donor_qnames.MM, 
            file_2 = recipient_qnames.UU, 
            out_file = "MMd_UUr.txt",
            resources = resources
        }

    call tools.inner_join as MUd_MMr_txt { 
        input: 
            file_1 = donor_qnames.MU, 
            file_2 = recipient_qnames.MM, 
            out_file = "MUd_MMr.txt",
            resources = resources
        }

    call tools.inner_join as MUd_MUr_txt { 
        input: 
            file_1 = donor_qnames.MU, 
            file_2 = recipient_qnames.MU, 
            out_file = "MUd_MUr.txt",
            resources = resources
        }

    call tools.inner_join as MUd_UMr_txt { 
        input: 
            file_1 = donor_qnames.MU, 
            file_2 = recipient_qnames.UM, 
            out_file = "MUd_UMr.txt",
            resources = resources
        }

    call tools.inner_join as MUd_UUr_txt { 
        input: 
            file_1 = donor_qnames.MU, 
            file_2 = recipient_qnames.UU, 
            out_file = "MUd_UUr.txt",
            resources = resources
        }

    call tools.inner_join as UMd_MMr_txt { 
        input: 
            file_1 = donor_qnames.UM, 
            file_2 = recipient_qnames.MM, 
            out_file = "UMd_MMr.txt",
            resources = resources
        }

    call tools.inner_join as UMd_MUr_txt { 
        input: 
            file_1 = donor_qnames.UM, 
            file_2 = recipient_qnames.MU, 
            out_file = "UMd_MUr.txt",
            resources = resources
        }

    call tools.inner_join as UMd_UMr_txt { 
        input: 
            file_1 = donor_qnames.UM, 
            file_2 = recipient_qnames.UM, 
            out_file = "UMd_UMr.txt",
            resources = resources
        }

    call tools.inner_join as UMd_UUr_txt { 
        input: 
            file_1 = donor_qnames.UM, 
            file_2 = recipient_qnames.UU, 
            out_file = "UMd_UUr.txt",
            resources = resources
        }

    call tools.inner_join as UUd_MMr_txt { 
        input: 
            file_1 = donor_qnames.UU, 
            file_2 = recipient_qnames.MM, 
            out_file = "UUd_MMr.txt",
            resources = resources
        }

    call tools.inner_join as UUd_MUr_txt { 
        input: 
            file_1 = donor_qnames.UU, 
            file_2 = recipient_qnames.MU, 
            out_file = "UUd_MUr.txt",
            resources = resources
        }

    call tools.inner_join as UUd_UMr_txt { 
        input: 
            file_1 = donor_qnames.UU, 
            file_2 = recipient_qnames.UM, 
            out_file = "UUd_UMr.txt",
            resources = resources
        }

    call tools.inner_join as UUd_UUr_txt { 
        input: 
            file_1 = donor_qnames.UU, 
            file_2 = recipient_qnames.UU, 
            out_file = "UUd_UUr.txt",
            resources = resources
        }


    # Subset BAM files on the inner joins 
    # between donor/recipient
    # for donor/recipient BAM files
    ## Donor
    call picard.filter_reads as donor_MMd_MMr { 
        input: 
            file = donor_MM, 
            txt = MMd_MMr_txt.out, 
            out_file = "reads-to-donor_MMd_MMr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_MMd_MUr { 
        input: 
            file = donor_MM, 
            txt = MMd_MUr_txt.out, 
            out_file = "reads-to-donor_MMd_MUr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_MMd_UMr { 
        input: 
            file = donor_MM, 
            txt = MMd_UMr_txt.out, 
            out_file = "reads-to-donor_MMd_UMr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_MMd_UUr { 
        input: 
            file = donor_MM, 
            txt = MMd_UUr_txt.out, 
            out_file = "reads-to-donor_MMd_UUr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_MUd_MMr { 
        input: 
            file = donor_MU, 
            txt = MUd_MMr_txt.out, 
            out_file = "reads-to-donor_MUd_MMr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_MUd_MUr { 
        input: 
            file = donor_MU, 
            txt = MUd_MUr_txt.out, 
            out_file = "reads-to-donor_MUd_MUr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_MUd_UMr { 
        input: 
            file = donor_MU, 
            txt = MUd_UMr_txt.out, 
            out_file = "reads-to-donor_MUd_UMr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_MUd_UUr { 
        input: 
            file = donor_MU, 
            txt = MUd_UUr_txt.out, 
            out_file = "reads-to-donor_MUd_UUr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_UMd_MMr { 
        input: 
            file = donor_UM, 
            txt = UMd_MMr_txt.out, 
            out_file = "reads-to-donor_UMd_MMr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_UMd_MUr { 
        input: 
            file = donor_UM, 
            txt = UMd_MUr_txt.out, 
            out_file = "reads-to-donor_UMd_MUr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_UMd_UMr { 
        input: 
            file = donor_UM, 
            txt = UMd_UMr_txt.out, 
            out_file = "reads-to-donor_UMd_UMr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_UMd_UUr { 
        input: 
            file = donor_UM, 
            txt = UMd_UUr_txt.out, 
            out_file = "reads-to-donor_UMd_UUr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_UUd_MMr { 
        input: 
            file = donor_UU, 
            txt = UUd_MMr_txt.out, 
            out_file = "reads-to-donor_UUd_MMr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_UUd_MUr { 
        input: 
            file = donor_UU, 
            txt = UUd_MUr_txt.out, 
            out_file = "reads-to-donor_UUd_MUr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_UUd_UMr { 
        input: 
            file = donor_UU, 
            txt = UUd_UMr_txt.out, 
            out_file = "reads-to-donor_UUd_UMr.bam",
            resources = resources
        }

    call picard.filter_reads as donor_UUd_UUr { 
        input: 
            file = donor_UU, 
            txt = UUd_UUr_txt.out, 
            out_file = "reads-to-donor_UUd_UUr.bam",
            resources = resources
        }


    ## Recipient
    call picard.filter_reads as recipient_MMd_MMr { 
        input: 
            file = recipient_MM, 
            txt = MMd_MMr_txt.out, 
            out_file = "reads-to-recipient_MMd_MMr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_MUd_MMr { 
        input: 
            file = recipient_MM, 
            txt = MUd_MMr_txt.out, 
            out_file = "reads-to-recipient_MUd_MMr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_UMd_MMr { 
        input: 
            file = recipient_MM, 
            txt = UMd_MMr_txt.out, 
            out_file = "reads-to-recipient_UMd_MMr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_UUd_MMr { 
        input: 
            file = recipient_MM, 
            txt = UUd_MMr_txt.out, 
            out_file = "reads-to-recipient_UUd_MMr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_MMd_MUr { 
        input: 
            file = recipient_MU, 
            txt = MMd_MUr_txt.out, 
            out_file = "reads-to-recipient_MMd_MUr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_MUd_MUr { 
        input: 
            file = recipient_MU, 
            txt = MUd_MUr_txt.out, 
            out_file = "reads-to-recipient_MUd_MUr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_UMd_MUr { 
        input: 
            file = recipient_MU, 
            txt = UMd_MUr_txt.out, 
            out_file = "reads-to-recipient_UMd_MUr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_UUd_MUr { 
        input: 
            file = recipient_MU, 
            txt = UUd_MUr_txt.out, 
            out_file = "reads-to-recipient_UUd_MUr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_MMd_UMr { 
        input: 
            file = recipient_UM, 
            txt = MMd_UMr_txt.out, 
            out_file = "reads-to-recipient_MMd_UMr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_MUd_UMr { 
        input: 
            file = recipient_UM, 
            txt = MUd_UMr_txt.out, 
            out_file = "reads-to-recipient_MUd_UMr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_UMd_UMr { 
        input: 
            file = recipient_UM, 
            txt = UMd_UMr_txt.out, 
            out_file = "reads-to-recipient_UMd_UMr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_UUd_UMr { 
        input: 
            file = recipient_UM, 
            txt = UUd_UMr_txt.out, 
            out_file = "reads-to-recipient_UUd_UMr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_MMd_UUr { 
        input: 
            file = recipient_UU, 
            txt = MMd_UUr_txt.out, 
            out_file = "reads-to-recipient_MMd_UUr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_MUd_UUr { 
        input: 
            file = recipient_UU, 
            txt = MUd_UUr_txt.out, 
            out_file = "reads-to-recipient_MUd_UUr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_UMd_UUr { 
        input: 
            file = recipient_UU, 
            txt = UMd_UUr_txt.out, 
            out_file = "reads-to-recipient_UMd_UUr.bam",
            resources = resources
        }

    call picard.filter_reads as recipient_UUd_UUr { 
        input: 
            file = recipient_UU, 
            txt = UUd_UUr_txt.out, 
            out_file = "reads-to-recipient_UUd_UUr.bam",
            resources = resources
        }

    Array[File?] all_bams = [
        donor_MMd_MMr.out, donor_MMd_MUr.out, 
        donor_MMd_UMr.out, donor_MMd_UUr.out, 
        donor_MUd_MMr.out, donor_MUd_MUr.out, 
        donor_MUd_UMr.out, donor_MUd_UUr.out, 
        donor_UMd_MMr.out, donor_UMd_MUr.out, 
        donor_UMd_UMr.out, donor_UMd_UUr.out, 
        donor_UUd_MMr.out, donor_UUd_MUr.out, 
        donor_UUd_UMr.out, donor_UUd_UUr.out, 
        recipient_MMd_MMr.out, recipient_MUd_MMr.out, 
        recipient_UMd_MMr.out, recipient_UUd_MMr.out, 
        recipient_MMd_MUr.out, recipient_MUd_MUr.out, 
        recipient_UMd_MUr.out, recipient_UUd_MUr.out, 
        recipient_MMd_UMr.out, recipient_MUd_UMr.out, 
        recipient_UMd_UMr.out, recipient_UUd_UMr.out, 
        recipient_MMd_UUr.out, recipient_MUd_UUr.out, 
        recipient_UMd_UUr.out, recipient_UUd_UUr.out]

    output {
        Array[File] bams = select_all(all_bams)
    }

}




