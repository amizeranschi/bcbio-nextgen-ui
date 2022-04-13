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
bcbio_workflow_dir="${bcbio_runs}${workflow_name%?}"
bcbio_runs_input="${bcbio_runs}${workflow_name%?}/input"


echo ""
echo " --- [$(date +"%F %R")] Starting the variant calling workflow"
echo " --- [$(date +"%F %R")] Using Python version: $(python --version)"
echo " --- [$(date +"%F %R")] Current PATH variable: ${PATH}"
echo " --- [$(date +"%F %R")] Using configuration file: " $1


# TODO restore/resume a previous analysis


# Perform variant calling
## Perform variant calling if the output (final) VCF file does not already exist:
if [ ! -f ${bcbio_runs_final}-gatk-haplotype*.vcf.gz ]; then
    echo " --- [$(date +"%F %R")] Running variant calling"
    cd ${bcbio_runs_input}
    ## If there's no variant_regions.bed file provided, create one with the full chromosomes of the bcbio genome
    if [[ ! -f full_chr_${bcbio_genome}.bed ]]; then
        if [[ ! -f ${bcbio_install_path?%}/genomes/${bcbio_species%?}/${bcbio_genome%?}/seq/${bcbio_genome%?}.fa.fai ]]; then
            echo " --- [$(date +"%F %R")] Creating an index for the genome ${bcbio_genome%?}"
            faidx ${bcbio_install_path%?}/genomes/${bcbio_species%?}/${bcbio_genome%?}/seq/${bcbio_genome%?}.fa
        fi
        
        echo " --- [$(date +"%F %R")] Creating a BED file for the full chromosomes of genome ${bcbio_genome%?}"
        # touch full_chr_${bcbio_genome%?}.bed
        awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${bcbio_install_path%?}/genomes/${bcbio_species%?}/${bcbio_genome%?}/seq/${bcbio_genome%?}.fa.fai > full_chr_${bcbio_genome%?}.bed
    fi

    if [[ ! -f variant_regions.bed ]]; then
        cp full_chr_${bcbio_genome%?}.bed variant_regions.bed
    fi
    # cd ${bcbio_runs_input}
    bcbio_nextgen.py -w template gatk-variant.yaml varcall-analysis.csv *.gz

    # TODO exclusion of low complexity regions

    # copy the variant_regiond.bed file to the config directory
    cp ${bcbio_runs_input}/variant_regions.bed ${bcbio_runs_input}/varcall-analysis/config

    # go to work dir
    cd ${bcbio_runs_input}/varcall-analysis/work
    # TODO cluster usage or not
    echo " --- [$(date +"%F %R")] Running bcbio-nextgen locally, using ${bcbio_main_cores%?} CPU cores"
    bcbio_nextgen.py ${bcbio_runs_input}/varcall-analysis/config/varcall-analysis.yaml -n ${bcbio_main_cores%?}

else
    echo " --- [$(date +"%F %R")] Skipping variant calling because the VCF file exists: "${bcbio_runs_final}-gatk-haplotype*.vcf.gz

fi
## TODO perform variant annotation
