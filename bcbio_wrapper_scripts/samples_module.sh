#!/bin/bash

##########################################################################################################################################################################################
                                                                            # SAMPLES MODULE #
##########################################################################################################################################################################################

## create a list of sample IDs
IFS=', ' read -r -a sample_list <<< ${bcbio_samples%?}
## create a list of name samples
IFS=', ' read -r -a sample_name_list <<< ${bcbio_samples_fastq%?}

## get list lengths to determine from the user if the sample is multiple
count_samples="${#sample_list[@]}"
count_user_name_samples="${#sample_name_list[@]}"

## get the number of samples from fastq
number_of_samples=$((${count_user_name_samples}/${count_samples}))

## if user wants to download data 
if [[ ${bcbio_download_samples%?} = "yes" ]]; then

   ## download samples in input directory  
   cd ${bcbio_runs_input}
   ## counter to store the index
   cnt=0
   ## download, rename and gzip
   echo " --- [$(date +"%F %R")] Downloading samples using sra-tools in ${bcbio_runs_input}."
   echo " --- [$(date +"%F %R")] Renaming and bgzipping files."

  for sample in ${sample_list[@]}
  do
      fasterq-dump --split-files -O . -t . ${sample}
      if [[ ${number_of_samples} = 1 ]]; then
         ## rename samples as user input
         mv ${sample}.fastq ${sample_name_list[${cnt}]}.fastq
         ## bgzip the samples
         bgzip -c ${sample_name_list[cnt]}.fastq > ${sample_name_list[${cnt}]}.fastq.gz
      fi
      if [[ ${number_of_samples} = 2 ]]; then
         ## rename samples as user input
         mv ${sample}_1.fastq ${sample_name_list[$((${cnt}*2))]}.fastq
         mv ${sample}_2.fastq ${sample_name_list[$((${cnt}*2+1))]}.fastq
         ## bgzip the samples
         bgzip -c ${sample_name_list[$((${cnt}*2))]}.fastq > ${sample_name_list[$((${cnt}*2))]}.fastq.gz
         bgzip -c ${sample_name_list[$((${cnt}*2+1))]}.fastq > ${sample_name_list[$((${cnt}*2+1))]}.fastq.gz
      fi

      ## TODO compute for 3 samples

      cnt=$((${cnt} + 1))
   done
   
   ## cleanup
   rm *.fastq

fi

# for data already on the disk
if [[ ${bcbio_download_samples%?} = "no" ]]; then
   ## go to source path where the samples are stored on the system
   cd ${bcbio_path_to_samples_on_sys%?}
   ## create symlinks in the input directory for the samples
   echo "--- [$(date +"%F %R")] Create symlinks in the input directory for the samples. "
   for FILE in *
   do
      file_name=$(echo "${FILE}" | cut -f 1 -d '.')
      for val in ${sample_name_list[@]}
      do
         if [[ ${file_name} = ${val} ]]; then
            echo "${FILE}"
            cp ${FILE} ${bcbio_runs_input}
         fi
      done
   done
fi

echo "--- [$(date +"%F %R")] Preparing configuration files for bcbio_nextgen in ${bcbio_workflow_dir}"

if [[ ${bcbio_workflow%?} = "variant_calling" ]]; then
   bash ${path_to_scripts}/config_module.sh 
fi

if [[ ${bcbio_workflow%?} = "atac_seq" ]]; then
   bash ${path_to_scripts}/config_module.sh
fi

if [[ ${bcbio_workflow%?} = "bulk_rna_seq" ]]; then
   bash ${path_to_scripts}/config_module.sh
fi
