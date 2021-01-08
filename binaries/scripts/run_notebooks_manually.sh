for SRR_NAME in SRR5090597 SRR5090599
do
	sudo papermill --prepare-only \
    tertiary-donor-to-recipient.ipynb ${SRR_NAME}-donor-to-recipient.ipynb \
    -p srr_name ${SRR_NAME} \
    -p donor_name "hpv16" \
    -p recipient_name "USCShg38" \
    -p inputs_folder "notebook_test_data/" \
    -p donor_ref_genome "../../wallaby/data/ref_genomes/hpv/HPV16.fasta" \
    -p recipient_ref_genome "../../wallaby/data/ref_genomes/human/USCS.hg38.fasta" \
    && jupyter nbconvert --execute --to html --output-dir hpv16_summary ${SRR_NAME}-donor-to-recipient.ipynb \
    --TemplateExporter.exclude_input=True --no-prompt --no-input
done





for SRR_NAME in SRR12091998 SRR12091997 SRR12091996 SRR12091995 SRR12091994 SRR12091993 SRR12091992 SRR12091991 SRR12091990

do
	sudo papermill --prepare-only \
    tertiary-donor-to-recipient.ipynb ${SRR_NAME}-donor-to-recipient.ipynb \
    -p srr_name ${SRR_NAME} \
    -p donor_name "sars-cov-2" \
    -p recipient_name "USCShg38" \
    -p inputs_folder "GSE153277/" \
    -p donor_ref_genome "../../wallaby/data/ref_genomes/sars-cov-2/wuhan-hu-1.fasta" \
    -p recipient_ref_genome "../../wallaby/data/ref_genomes/human/USCS.hg38.fasta" \
    && jupyter nbconvert --execute --to html --output-dir GSE153277_summary ${SRR_NAME}-donor-to-recipient.ipynb \
    --TemplateExporter.exclude_input=True --no-prompt --no-input
done






for SRR_NAME in SRR12507513 SRR12507514 SRR12507515 SRR12507516 SRR12507517 SRR12507518 SRR12507519 SRR12507520 SRR12507521 SRR12507522 SRR12507523 SRR12507524 SRR12507525 SRR12507526 SRR12507527 SRR12507528 SRR12507529 SRR12507530 SRR12507531 SRR12507532 SRR12507533 SRR12507534 SRR12507535 SRR12507536 SRR12507537 SRR12507538 SRR12507539 SRR12507540 SRR12507541 SRR12507542 SRR12507543 SRR12507544 SRR12507545 SRR12507546 SRR12507547 SRR12507548 SRR12507549 SRR12507550 SRR12507551 SRR12507552

do
	sudo papermill --prepare-only \
    tertiary-donor-to-recipient.ipynb ${SRR_NAME}-donor-to-recipient.ipynb \
    -p srr_name ${SRR_NAME} \
    -p donor_name "sars-cov-2" \
    -p recipient_name "USCShg38" \
    -p inputs_folder "GSE156754/" \
    -p donor_ref_genome "../../wallaby/data/ref_genomes/sars-cov-2/wuhan-hu-1.fasta" \
    -p recipient_ref_genome "../../wallaby/data/ref_genomes/human/USCS.hg38.fasta" \
    && jupyter nbconvert --execute --to html --output-dir GSE156754_summary ${SRR_NAME}-donor-to-recipient.ipynb \
    --TemplateExporter.exclude_input=True --no-prompt --no-input
done






for SRR_NAME in SRR12134528 SRR12134529 SRR12134530 SRR12134531 SRR12134532 SRR12134533 SRR12134534 SRR12134535 SRR12134536 SRR12134537 SRR12134538 SRR12134539 SRR12134540 SRR12134541 SRR12134542 SRR12134543 SRR12134544 SRR12134545

do
	sudo papermill --prepare-only \
    tertiary-donor-to-recipient.ipynb ${SRR_NAME}-donor-to-recipient.ipynb \
    -p srr_name ${SRR_NAME} \
    -p donor_name "sars-cov-2" \
    -p recipient_name "USCShg38" \
    -p inputs_folder "GSE153684/" \
    -p donor_ref_genome "../../wallaby/data/ref_genomes/sars-cov-2/wuhan-hu-1.fasta" \
    -p recipient_ref_genome "../../wallaby/data/ref_genomes/human/USCS.hg38.fasta" \
    && jupyter nbconvert --execute --to html --output-dir GSE153684_summary ${SRR_NAME}-donor-to-recipient.ipynb \
    --TemplateExporter.exclude_input=True --no-prompt --no-input
done

