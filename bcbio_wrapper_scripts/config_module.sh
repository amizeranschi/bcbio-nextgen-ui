#!/bin/bash
## temorary TO DO CHANGE 
curr_dir="/export/home/acs/stud/m/maria.nastase0912/bcbio_nextgen_usability_improvements/bcbio_wrapper_scripts"


# variables to store yaml template urls
variant_calling_yaml="https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/templates/gatk-variant.yaml"
# PARSE THE INPUT YAML CONFIG FILE

## Function for reading from a simple YAML file, adapted from: https://stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script/21189044#21189044
## Modified the printf line in the awk script code below to add the word "export ", such that it creates environment variables accessible from the rest of the analysis scripts
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("export %s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

echo "using configuration file: " $1


## read the config file and create the variables using the "bcbio_" prefix
eval $(parse_yaml $1 "bcbio_")

echo "" 
echo "## printing the whole contents of the file:"

parse_yaml $1 "bcbio_"


bcbio_runs="${HOME}/bcbio_runs/"
workflow_name="workflow_${bcbio_workflow}"
bcbio_workflow_dir="${bcbio_runs}${workflow_name%?}"
bcbio_runs_input="${bcbio_runs}${workflow_name%?}/input"

# Go to config directory of the workflow

cd ${bcbio_runs_input}

# download a template yaml file, describing the analysis

if [[ ${bcbio_workflow%?} == "variant_calling" ]]; then
    # if [[ ! -f *.yaml ]]; then
    if [ -x "$(command -v wget)" ]; then
        wget ${variant_calling_yaml}
     else
        curl -L -C - -O ${variant_calling_yaml}
        # fi
    fi
   ## edit the settings in the illumina-rnaseq.yaml file, to make the analysis work for our files:

    sed -i 's/genome_build: hg38/genome_build: '${bcbio_genome}/ gatk-variant.yaml

    sed -i 's/aligner: bwa/aligner: bowtie2/' gatk-variant.yaml

    sed -i 's/recalibrate: gatk/# recalibrate: gatk/' gatk-variant.yaml

    # copy csv file from the location given in input to the config directory

    cp ${bcbio_csv_file_path%?} ${bcbio_runs_input}

    bash ${curr_dir}/run_variant_calling.sh ${curr_dir}/infos.yaml

fi