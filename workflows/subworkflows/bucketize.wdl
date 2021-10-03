version 1.0

import "../tasks/samtools.wdl" as samtools

workflow main {

    input {
        File bam
        String base_filename
        Resources resources
    }

    # Read 1 is Mapped and Read 2 is Mapped in/to Ref Genome - MM
    call samtools.view as MM_bam { 
        input: 
            file = bam, 
            out_file = "~{base_filename}_MM.bam", 
            include = "3", exclude = "2048",
            resources = resources
    }

    call samtools.index as MM_bai { 
    input: 
        file = MM_bam.out, 
        out_file = "~{base_filename}_MM.bai", 
        resources = resources
    }

    # Read 1 is Mapped and Read 2 is Unmapped in/to Ref Genome - MU
    call samtools.view as MU_R1 { 
        input: 
            file = bam, 
            out_file = "~{base_filename}_MU_R1.bam", 
            include = "73", exclude = "4022",
            resources = resources
    }

    call samtools.view as MU_R2 { 
        input: 
            file = bam, 
            out_file = "~{base_filename}_MU_R2.bam", 
            include = "133", exclude = "3962",
            resources = resources
    }

    call samtools.merge as MU_bam { 
        input: 
            files = [MU_R1.out, MU_R2.out], 
            out_file = "~{base_filename}_MU.bam",
            resources = resources
    }

    call samtools.index as MU_bai { 
    input: 
        file = MU_bam.out, 
        out_file = "~{base_filename}_MU.bai", 
        resources = resources
    }

    # Read 1 is Unmapped and Read 2 is Mapped in/to Ref Genome - UM
    call samtools.view as UM_R1 { 
        input: 
            file = bam, 
            out_file = "~{base_filename}_UM_R1.bam", 
            include = "69", exclude = "4026",
            resources = resources
    }

    call samtools.view as UM_R2 { 
        input: 
            file = bam, 
            out_file = "~{base_filename}_UM_R2.bam", 
            include = "137", exclude = "3958",
            resources = resources
    }

    call samtools.merge as UM_bam { 
        input: 
            files = [UM_R1.out, UM_R2.out], 
            out_file = "~{base_filename}_UM.bam",
            resources = resources
    }

    call samtools.index as UM_bai { 
    input: 
        file = UM_bam.out, 
        out_file = "~{base_filename}_UM.bai", 
        resources = resources
    }

    # Read 1 is Unmapped and Read 2 is Unmapped in/to Ref Genome - UU
    #call samtools.view as UU_R1 { 
    #    input: 
    #        file = bam, 
    #        out_file = "~{base_filename}_UU_R1.bam", 
    #        include = "77",
    #        resources = resources
    #}

    #call samtools.view as UU_R2 { 
    #    input: 
    #        file = bam, 
    #        out_file = "~{base_filename}_UU_R2.bam", 
    #        include = "141",
    #        resources = resources
    #}

    #call samtools.merge as UU_bam { 
    #    input: 
    #        files = [UU_R1.out, UU_R2.out], 
    #        out_file = "~{base_filename}_UU.bam",
    #        resources = resources
    #}

    #call samtools.index as UU_bai { 
    #input: 
    #    file = UU_bam.out, 
    #    out_file = "~{base_filename}_UU.bai", 
    #    resources = resources
    #}

    output {
        File MM = MM_bam.out
        File MU = MU_bam.out
        File UM = UM_bam.out
        #File UU = UU_bam.out
        Array[File] bams = [
            MM_bam.out,
            MU_bam.out,
            UM_bam.out
            #UU_bam.out
        ]
        Array[File] bais = [
            MM_bai.out,
            MU_bai.out,
            UM_bai.out
            #UU_bai.out
        ]
    }

}


