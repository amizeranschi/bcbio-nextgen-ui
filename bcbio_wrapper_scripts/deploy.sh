#!/bin/bash

## Give excecution rights to every script
chmod +x environment/*.sh
chmod +x workflows/*.sh
#chmod +x downstreamAnalysis*/*.sh
chmod +x *.sh

##########################################################################################################################################################################################
                                                                            # Environment #
##########################################################################################################################################################################################

# CREATE MAIN WORKSPACE, VARIABLES AND SET PATHS
## Parse YAML configuration file to get all user's analsis data
source environment/parse_yaml.sh $1
## Set environment variables, links for installers and paths
source environment/set_environment_variables.sh

## Call the module for the setup of  bcbio_nextgen and the environment
## either install tool or create symlinks to already installed version
## Create Python environments to install the following tools
##                      * bioconda 
##                      * r-base 
##                      * r-xml 
##                      * bedops 
##                      * vcftools
##                      * sra-tools 
##                      * perl-net-ssleay 
##                      * gitentrez-direct
##                      * tassel
##                      * beagle
##                      * plink 
##                      * openssl=1.1.1l 
##                      * ensembl-vep=105
##                      * bcftools>=1.13
##                      * bcbio-monitor
##                      * pytz
##                      * python-dateutil
##                      * faststructure

bash environment/setup_environment_module.sh

export PATH="${bcbio_install_path}/extra3/bin:${bcbio_install_path}/anaconda/bin:${bcbio_install_path}/tools/bin:${bcbio_install_path}/extra2/bin${PATH:+:${PATH}}"
echo " --- [$(date +"%F %R")] The PATH IS: ${PATH}"

##########################################################################################################################################################################################
                                                                            # Workflow #
##########################################################################################################################################################################################

## GO TO SAMPLES MODULE
## Handle the samples for the upcoming analysis

bash ${path_to_scripts}/workflows/samples_module.sh

## Setup configuration files for the workflow
bash ${path_to_scripts}/workflows/config_module.sh 

## RUN WORKFLOW
if [[ ${bcbio_workflow} == "variant_calling" ]]; then
   bash  ${path_to_scripts}/workflows/run_variant_calling.sh
fi

if [[ ${bcbio_workflow} == "atac_seq" ]]; then
   bash  ${path_to_scripts}/workflows/run_atac_seq.sh
fi

if [[ ${bcbio_workflow} == "bulk_rna_seq" ]]; then
    bash  ${path_to_scripts}/workflows/run_bulk_rna_seq.sh
fi

##########################################################################################################################################################################################
                                                                            # Downstream Analysis #
##########################################################################################################################################################################################

## Run downstream analysis for variant calling: variant annotation and gene annotation
if [[ ${bcbio_workflow} == "variant_calling" ]]; then
    
    if [[ ${bcbio_variant_annotation} == "yes" ]]; then
        echo " --- [$(date +"%F %R")] Starting DOWNSTREAM ANALYSIS for VARIANT CALLING WORKFLOW."

        # run variant annotation if requested
        bash ${path_downstream_analysis}/variant_annotation.sh
        # perform gene annotation on the results from variant annotation
        echo " --- [$(date +"%F %R")] Starting gene annotation for Variant Calling Workflow, see output in ${path_to_scripts}/downstreamAnalysisVariantCalling."
        Rscript --vanilla ${path_downstream_analysis}/gene_annotation_variant_calling.R ${path_downstream_analysis} ${variant_annotation_dir}/${vcf_file_name}-vep.table ${bcbio_vep_species} ${gtf_file_location}
        

    fi

fi

## Run downstream analysis for atac_seq
if [[ ${bcbio_workflow} == "atac_seq" ]]; then
    # to do add path to peaks file
    echo " --- [$(date +"%F %R")] Starting DOWNSTREAM ANALYSIS for ATAC-seq/ChIP-seq WORKFLOW."
    Rscript --vanilla ${path_downstream_analysis}/chIP_seq-downstreamAnalysis.R ${path_downstream_analysis_bulk}  ${bcbio_vep_species} ${gtf_file_location}

fi

## Run downstream analysis for bulk_rna_seq
if [[ ${bcbio_workflow} == "bulk_rna_seq" ]]; then
    echo " --- [$(date +"%F %R")] Starting DOWNSTREAM ANALYSIS for BULK RNA-seq WORKFLOW."
    Rscript --vanilla ${path_downstream_analysis}/bulk_rna_seq-downstream_analysis.R ${path_downstream_analysis} ${counts_file} ${metadata_file} ${bcbio_vep_species} ${gtf_file_location}
    python3 create_json_downstream_page.py ${bcbio_workflow} ${path_downstream_analysis_bulk}
fi

mv ${path_downstream_analysis}/*.txt ${dowstreamResultsWorkflow}
mv ${path_downstream_analysis}/*.png ${dowstreamResultsWorkflow}
