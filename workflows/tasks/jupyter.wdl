version 1.0

import "structs/compute.wdl"

task run_notebook {
    input {
        Array[File] inputs
        String srr_name
        String donor_name
        File donor_ref_genome
        String recipient_name
        File recipient_ref_genome
        File notebook
        Resources resources
    }

    command <<<
        # workaround on bug: https://github.com/broadinstitute/cromwell/issues/4361
        # copies all files from the 'inputs' folder into the 'files' folder
        mkdir files
        export FILES="`pwd`/files"
        cp ~{sep=" " inputs} `echo $FILES`
        cp ~{donor_ref_genome} `echo $FILES`
        cp ~{recipient_ref_genome} `echo $FILES`
        papermill --prepare-only \
            ~{notebook} \
            ~{srr_name}-donor-to-recipient.ipynb \
            -p srr_names ~{srr_name} \
            -p donor_name ~{donor_name} \
            -p recipient_name ~{recipient_name} \
            -p inputs_folder `echo $FILES` \
            -p donor_ref_genome ~{donor_ref_genome} \
            -p recipient_ref_genome ~{recipient_ref_genome} \
            && jupyter nbconvert --execute --to html --output-dir \
                . ~{srr_name}-donor-to-recipient.ipynb \
                --TemplateExporter.exclude_input=True --no-prompt --no-input
        rm -rf `echo $FILES`
        
    >>>

    output {
        File out = "~{srr_name}-donor-to-recipient.html"
    }

    runtime {
        continueOnReturnCode: false
        cpu: resources.cpu
        memory: resources.memory_gb
        docker: "dformoso/jupyterlab:latest"
        disks: resources.disks
        zones: resources.zones
        preemptible: resources.preemptible
        maxRetries: resources.maxRetries
    }
}
