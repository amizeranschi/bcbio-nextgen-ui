#!/bin/bash

variant_annotation_dir="${bcbio_runs_input}/${action_name}/variant_annotation"
echo ""
echo "--- [$(date +"%F %R")] Starting the Variant Calling Workflow"
echo "--- [$(date +"%F %R")] Using Python version: $(python --version)"
# echo "--- [$(date +"%F %R")] Current PATH variable: ${PATH}"
echo "--- [$(date +"%F %R")] Using configuration from directory: " ${path_to_scripts}

##########################################################################################################################################################################################
                                                                            # VARIANT CALLING WORKFLOW #
##########################################################################################################################################################################################
# echo "${bcbio_install_path%?} AICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
# echo "${bcbio_runs_final} AICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
# echo "${bcbio_runs_input} AICIIIIIIIIIIIIIIIIIIIIIIIIIIIIII!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
# TODO restore/resume a previous analysis

# Perform variant calling
## Perform variant calling if the output (final) VCF file does not already exist:

if [ ! -f ${bcbio_runs_final}/*_${action_name}/*-gatk-haplotype*.vcf.gz ]; then
    echo "--- [$(date +"%F %R")] Running variant calling"
    cd ${bcbio_runs_input}
    ## If there's no variant_regions.bed file provided, create one with the full chromosomes of the bcbio genome
    if [[ ! -f full_chr_${bcbio_genome}.bed ]]; then
        if [[ ! -f ${genome_dir}/seq/${bcbio_genome%?}.fa.fai ]]; then
            echo "--- [$(date +"%F %R")] Creating an index for the genome ${bcbio_genome%?}"
            faidx ${genome_dir}/seq/${bcbio_genome%?}.fa
        fi
        
        echo "--- [$(date +"%F %R")] Creating a BED file for the full chromosomes of genome ${bcbio_genome%?}"
        # touch full_chr_${bcbio_genome%?}.bed
        awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${genome_dir}/seq/${bcbio_genome%?}.fa.fai > full_chr_${bcbio_genome%?}.bed
    fi

    if [[ ! -f variant_regions.bed ]]; then
        cp full_chr_${bcbio_genome%?}.bed variant_regions.bed
    fi
    # cd ${bcbio_runs_input}
    bcbio_nextgen.py -w template gatk-variant.yaml ${action_name}.csv *.gz

    # TODO exclusion of low complexity regions

    # copy the variant_regiond.bed file to the config directory
    cp ${bcbio_runs_input}/variant_regions.bed ${bcbio_workflow_config}

    # go to work dir
    cd ${bcbio_workflow_work}
    # TODO cluster usage or not

    echo "--- [$(date +"%F %R")] Running bcbio-nextgen locally, using ${bcbio_main_cores%?} CPU cores"
    bcbio_nextgen.py ${bcbio_workflow_config}/${action_name}.yaml -n ${bcbio_main_cores%?}

else
    echo "--- [$(date +"%F %R")] Skipping variant calling because the VCF file exists: "${bcbio_runs_final}/*_${action_name}/*-gatk-haplotype*.vcf.gz

fi

echo "--- [$(date +"%F %R")] Finished variant calling"

##########################################################################################################################################################################################
                                                                             # VARIANT ANNOTATION #
##########################################################################################################################################################################################

echo "--- [$(date +"%F %R")] STARTING VARIANT ANNOTATION"

if [[ ${bcbio_variant_annotation%?} == "yes" ]]; then
    # remove previous annotation and create variant annotation directory
    rm -rf ${variant_annotation_dir}
    mkdir ${variant_annotation_dir}
    cd ${variant_annotation_dir}

    ## copy the joint VCF file with a proper name for variant annotation
	if [ -f ${bcbio_runs_final}/*_${action_name}/*-gatk-haplotype*.vcf.gz ]; then
		cp -f  ${bcbio_runs_final}/*_${action_name}/*-gatk-haplotype*.vcf.gz ${action_name}-small-var.vcf.gz
	fi
    ## run small variants annotation if there is the case
    bash ${path_to_scripts}/small_variants_annotation.sh $1

    ## run structural variants annotation if there is the case
    bash ${path_to_scripts}/structural_variants_annotation.sh $1

fi


## print message for workflow completed
echo "--- [$(date +"%F %R")] Variant calling workflow is finished."
