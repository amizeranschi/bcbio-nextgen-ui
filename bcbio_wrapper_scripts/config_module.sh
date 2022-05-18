#!/bin/bash

## download a template yaml file, describing the analysis
cd ${bcbio_runs_input}

if [[ ${bcbio_workflow%?} == "variant_calling" ]]; then
   if [ -x "$(command -v wget)" ]; then
      wget ${variant_calling_yaml}
   else
      curl -L -C - -O ${variant_calling_yaml}
   fi
   ## edit the settings in the illumina-rnaseq.yaml file, to make the analysis work for our files:
   echo "--- [$(date +"%F %R")] Configuring yaml template file for the variant calling workflow"

   sed -i 's/genome_build: hg38/genome_build: '${bcbio_genome}/ gatk-variant.yaml
   sed -i 's/aligner: bwa/aligner: bowtie2/' gatk-variant.yaml
   sed -i 's/recalibrate: gatk/# recalibrate: gatk/' gatk-variant.yaml
   sed -i "/algorithm:$/a\      effects: false" gatk-variant.yaml

   ## copy csv file from the location given in input to the config directory
   
   echo "--- [$(date +"%F %R")] Copying csv  file provided by the user for the variant calling workflow"
   cp ${bcbio_csv_file_path%?} ${bcbio_runs_input}

   bash  ${path_to_scripts}/run_variant_calling.sh

fi

if [[ ${bcbio_workflow%?} == "atac_seq" ]]; then
   if [ -x "$(command -v wget)" ]; then
      wget --no-check-certificate ${atac_seq_yaml} -O atac-example.yaml
   else
      curl -L -C - -O ${atac_seq_yaml}
   fi

   ## edit the settings in the illumina-rnaseq.yaml file, to make the analysis work for our files:
   echo "--- [$(date +"%F %R")] Configuring yaml template file for the ATAC-seq workflow"

   sed -i 's/genome_build: hg38/genome_build: '${bcbio_genome}/ atac-example.yaml
   sed -i 's/aligner: bwaa/aligner: bowtie2/' atac-example.yaml
   # pip3 install pyyaml
   # python3 ${path_to_scripts}/add_to_yaml.py ${bcbio_runs_input}/atac-example.yaml

   # copy csv file from the location given in input to the config directory
   echo "--- [$(date +"%F %R")] Copying csv  file provided by the user for the ATAC-seq workflow"
   echo "--- [$(date +"%F %R")] For ChIP-seq, bcbio_nextgen requires batch and phenotype in the csv columns for the ATAC-seq workflow"

   cp ${bcbio_csv_file_path%?} ${bcbio_runs_input}

   bash  ${path_to_scripts}/run_atac_seq.sh

fi