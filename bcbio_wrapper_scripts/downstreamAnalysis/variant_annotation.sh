#!/bin/bash

##########################################################################################################################################################################################
                                                                             # VARIANT ANNOTATION #
##########################################################################################################################################################################################

echo " --- [$(date +"%F %R")] Running variant annotation in the following directory: ${path_downstream_analysis}"

# remove previous annotation and create variant annotation directory
if [ -d ${path_downstream_analysis} ]; then
   rm -rf ${path_downstream_analysis}
fi
mkdir ${path_downstream_analysis}
cd ${path_downstream_analysis}

## copy the joint VCF file with a proper name for variant annotation
if [ -f ${bcbio_runs_final}/*_${action_name}/*-gatk-haplotype*.vcf.gz ]; then
   cp -f  ${bcbio_runs_final}/*_${action_name}/*-gatk-haplotype*.vcf.gz ${action_name}-small-var.vcf.gz
fi
    
## run small variants annotation if there is the case
bash ${path_to_scripts}/downstreamAnalysis/small_variants_annotation.sh

## run structural variants annotation if there is the case
bash ${path_to_scripts}/downstreamAnalysis/structural_variants_annotation.sh
