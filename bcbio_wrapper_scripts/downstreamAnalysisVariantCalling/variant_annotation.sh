#!/bin/bash

##########################################################################################################################################################################################
                                                                             # VARIANT ANNOTATION #
##########################################################################################################################################################################################

echo " --- [$(date +"%F %R")] STARTING VARIANT ANNOTATION"

if [[ ${bcbio_variant_annotation} == "yes" ]]; then
    # remove previous annotation and create variant annotation directory
    rm -rf ${variant_annotation_dir}
    mkdir ${variant_annotation_dir}
    cd ${variant_annotation_dir}

    ## copy the joint VCF file with a proper name for variant annotation
	if [ -f ${bcbio_runs_final}/*_${action_name}/*-gatk-haplotype*.vcf.gz ]; then
		cp -f  ${bcbio_runs_final}/*_${action_name}/*-gatk-haplotype*.vcf.gz ${action_name}-small-var.vcf.gz
	fi
    
    ## run small variants annotation if there is the case
    bash ${path_downstream_analysis}/small_variants_annotation.sh

    ## run structural variants annotation if there is the case
    bash ${path_downstream_analysis}/structural_variants_annotation.sh
fi