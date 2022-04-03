#!/bin/bash
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
bcbio_workflow="${bcbio_runs}${workflow_name%?}"
bcbio_runs_input="${bcbio_runs}${workflow_name%?}/input"
bcbio_runs_config="${bcbio_runs}${workflow_name%?}/config"
bcbio_runs_work="${bcbio_runs}${workflow_name%?}/work"
bcbio_runs_final="${bcbio_runs}${workflow_name%?}/final"

cd ${bcbio_workflow}

# SAMPLES MODULE
## if user wants to download data 
if [[ ${bcbio_download_samples%?} == "yes" ]]; then

    ## create a list of sample IDs
  IFS=', ' read -r -a sample_list <<< ${bcbio_samples%?}
  ## download samples in input directory
  cd ${bcbio_runs_input}
  
  for sample in ${sample_list[@]}
  do
        echo "$sample"
       # fasterq-dump --split-files -O . -t . ${sample}
        ## TODO figure out a way to automatically change names or determine the types of the samples
        ## TODO rename and gzip the samples, then remove the .fastq files
  done

  ## for now renaming manually
  mv SRR6059150.fastq 500-F-Rep3.fastq

    gzip -c 500-F-Rep3.fastq > 500-F-Rep3.fastq.gz

    mv SRR6059151.fastq 500-I-Rep3.fastq

    gzip -c 500-I-Rep3.fastq > 500-I-Rep3.fastq.gz

    mv SRR6783014.fastq 500-I-Rep1.fastq

    gzip -c 500-I-Rep1.fastq > 500-I-Rep1.fastq.gz

    mv SRR6783015.fastq 500-F-Rep2.fastq

    gzip -c 500-F-Rep2.fastq > 500-F-Rep2.fastq.gz

    mv SRR6783016.fastq 500-I-Rep2.fastq

    gzip -c 500-I-Rep2.fastq > 500-I-Rep2.fastq.gz

    mv SRR6784354.fastq 500-F-Rep1.fastq

    gzip -c 500-F-Rep1.fastq > 500-F-Rep1.fastq.gz

    rm *.fastq

    echo " --- [$(date +"%F %R")] Preparing configuration files for bcbio_nextgen in ${bcbio_workflow%?}/config"

    # TODO Go to config directory and get template yaml file
    # TODO copy .csv file provided by user to config
    # TODO config yaml file for the exisiting set of samples
    # TODO check for .bed file, if none, create one using bcbio genome
    # TODO exclusion of low compatibility regions
fi