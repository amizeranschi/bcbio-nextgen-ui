#!/bin/bash

##########################################################################################################################################################################################
                                                                             # VARIANT ANNOTATION #
##########################################################################################################################################################################################

echo " --- [$(date +"%F %R")] Running variant annotation in the following directory: ${variant_annotation_dir}"

# remove previous annotation and create variant annotation directory
if [ -d ${variant_annotation_dir} ]; then
   rm -rf ${variant_annotation_dir}
fi
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
