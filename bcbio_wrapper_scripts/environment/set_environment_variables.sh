#!/bin/bash

#####################################################################################################################################################
                                            # VARIABLES USED INSIDE EACH WORKFLOW FOR EASE OF ACCESS #
#####################################################################################################################################################

if [[ ${bcbio_existing_version} = "yes" ]]; then
    export bcbio_install_path="${bcbio_path_to_existing}"
fi

# path to all workflow analyses on the system
export bcbio_runs="${HOME}/bcbio_runs/"
# path to the current workflow directory
export workflow_name="workflow_${bcbio_workflow}"
export bcbio_workflow_dir="${bcbio_runs}${workflow_name}"
# path to the input files and anlysis directory
export bcbio_runs_input="${bcbio_runs}${workflow_name}/input"
# extraction of the name of the csv file to keep for the analysis flow name
export action_name=$(echo ${bcbio_csv_file_path##*/} | cut -f 1 -d '.')


# create the directories
if [ ! -d ${bcbio_runs} ]; then
    mkdir ${bcbio_runs}
fi
if [ ! -d ${bcbio_workflow_dir} ]; then
    mkdir ${bcbio_workflow_dir}
fi
if [ ! -d ${bcbio_runs_input} ]; then
    mkdir ${bcbio_runs_input}
fi


# directory where bcbio genome is stored
# seq_dir_genome="${bcbio_install_path}/genomes/${bcbio_species}/${bcbio_genome}/seq"
export genome_dir="${bcbio_install_path}/genomes/${bcbio_species}/${bcbio_genome}"
export gtf_file_location=" ${genome_dir}/rnaseq/ref-transcripts.gtf"

# directories for the analysis directory tree of bcbio
export bcbio_runs_final="${bcbio_runs_input}/${action_name}/final"
export bcbio_workflow_config="${bcbio_runs_input}/${action_name}/config"
export bcbio_workflow_work="${bcbio_runs_input}/${action_name}/work"

## Store current path to the scripts
export path_to_scripts=$PWD
export path_to_web="${path_to_scripts}/web"
export path_downstream_analysis="${bcbio_runs_input}/${action_name}/downstreamAnalysis"
mkdir ${path_downstream_analysis}

# set variables for variant annotation and gene annotation in downstream analysis
if [[ ${bcbio_workflow} == "variant_calling" ]]; then
    vcf_file="${action_name}-small-var.vcf.gz"
    vcf_file_name=$(echo "${vcf_file}" | cut -f 1 -d '.')
fi

export counts_file="${bcbio_runs_final}/*${action_name}/counts/*.csv"
export metadata_file="${bcbio_runs_final}/*${action_name}/metadata.csv"

#####################################################################################################################################################
                                                     # URLS AND USEFUL TEMPLATES #
#####################################################################################################################################################
## variables to store yaml template urls
export variant_calling_yaml="https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/templates/gatk-variant.yaml"
export atac_seq_yaml="http://s3.amazonaws.com/bcbio-nextgen/atac_userstory_data/atac-example.yaml"
export bulk_rna_seq_yaml="https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/rnaseq-seqc.yaml"

## upgrade the extra3 environment from python3.7 to python3.9
#export miniconda3="https://repo.anaconda.com/miniconda/Miniconda3-py37_4.12.0-Linux-x86_64.sh"
#export miniconda3="https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh"
export miniconda3="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

export miniconda2="https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh"
export install_script_bcbio="https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py"

#  GET UTILITARIES AND SCRIPTS
## check if wget is available. if not, use curl instead
echo " --- [$(date +"%F %R")] Downloading utilitaries for bcbio_nexgen and miniconda installation"

if [ -x "$(command -v wget)" ]; then
    wget ${miniconda3}
    wget ${miniconda2}
    wget ${install_script_bcbio}
else
    curl -L -C - -O ${miniconda3}
    curl -L -C - -O ${miniconda2}
    curl -L -C - -O ${install_script_bcbio}
fi
