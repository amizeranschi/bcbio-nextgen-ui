#!/bin/bash

##########################################################################################################################################################################################
                                                                         # ALREADY INSTALLED VERSION OF BCBIO #
##########################################################################################################################################################################################

## IF BCBIO IS ALREADY INSTALLED
if [[ ${bcbio_existing_version%?} = "yes" ]]; then

   echo " --- [$(date +"%F %R")] bcbio_nextgen already installedy installed on the system."
   echo " --- [$(date +"%F %R")] Use bcbio_nextgen already installed version."

   cd ${bcbio_install_path%?}
   bash ${path_to_scripts}/setup_python2_env.sh
   bash ${path_to_scripts}/setup_python3_env.sh
   
   ## set the path with all the utils
   echo " --- [$(date +"%F %R")] Setting the PATH for Python 3 and Python2 environment installation."
   export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:${bcbio_install_path%?}/extra3/bin:${bcbio_install_path%?}/extra2/bin:${PATH}
   echo " --- [$(date +"%F %R")] The PATH IS: ${PATH}"
   ln -s ${bcbio_install_path%?}/anaconda/envs/python2/bin/python ${bcbio_install_path%?}/tools/bin/python2

   ## check if a genome is installed on the system
   echo " --- [$(date +"%F %R")] Check if there is any genome installed on the system"
   if [[ ! -d ${bcbio_install_path%?}/genomes/ ]]; then
      echo " --- [$(date +"%F %R")] No genome directory found."
      echo " --- [$(date +"%F %R")] Installing genome ${bcbio_genome%?} in the ${bcbio_install_path%?}/genomes path"
   fi
   ## check if the requested genome is installed on the system
   ## if not, install it
   if [[ ! -d ${bcbio_install_path%?}/genomes/${bcbio_species%?}/${bcbio_genome%?} ]]; then
      echo " --- [$(date +"%F %R")] The requested genome ${bcbio_genome%?} was not found."
      if [[ ${bcbio_annotated_species%?} = "yes" ]]; then
         # TODO check the arguments
         echo " --- [$(date +"%F %R")] Installing genome ${bcbio_genome%?} in the ${bcbio_install_path%?}/genomes path"
         bcbio_nextgen.py upgrade -u skip --genomes ${bcbio_genome%?} --datatarget variation --datatarget rnaseq --datatarget smallrna \
            --aligners bwa --aligners bowtie2 --aligners star --isolate --cores ${bcbio_total_cores%?}
      else
         ## install from custom genome, bcbio needs to have an already installed genome in order to install a custom genome
         echo " --- [$(date +"%F %R")] Installing Bcbio-nextgen with genome and transcriptome annotations for the sacCer3 reference"
         echo " --- [$(date +"%F %R")] This uses a small amount of storage and is required for creating a baseline Bcbio-nextgen configuration"
         
         ## install custom genome
         bcbio_nexgen.py upgrade -u skip --genomes sacCer3 --datatarget variation --datatarget rnaseq --datatarget smallrna \
            --aligners bwa --aligners bowtie2 --aligners hisat2 --aligners star --isolate --cores ${bcbio_total_cores} --mamba
         
         echo " --- [$(date +"%F %R")] Installing custom genome and transcriptome annotations from the user-provided FASTA and GTF files for the ${bcbio_genome} reference"
         
         bcbio_setup_genome.py -f ${bcbio_genome_fasta%?} -g ${bcbio_transcriptome_gtf%?} -i bwa bowtie2 star seq \
            -n ${bcbio_species%?} -b ${bcbio_genome%?} -c ${bcbio_total_cores%?} --buildversion ${bcbio_genome%?}
      fi

   fi

   echo " --- [$(date +"%F %R")] The requested genome ${bcbio_genome%?} was found at ${bcbio_install_path%?}."
   bcbio_conda install -y biobambam=2.0.87 -c bioconda
fi

##########################################################################################################################################################################################
                                                                         # BCBIO INSTALLATION #
##########################################################################################################################################################################################

if [[ ${bcbio_existing_version%?} = "no" ]]; then
   echo " --- [$(date +"%F %R")] bcbio_nextgen not on the system."
   echo " --- [$(date +"%F %R")] bcbio_nextgen will start installation."
   mkdir ${bcbio_install_path%?}
   cd ${bcbio_install_path%?}
   bash ${path_to_scripts}/setup_python2_env.sh
   bash ${path_to_scripts}/setup_python3_env.sh
   bash ${path_to_scripts}/install_bcbio_nextgen.sh

   ## set the path with all the utils
   echo " --- [$(date +"%F %R")] Setting the PATH for Python 3 and Python2 environment installation."
   export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:${bcbio_install_path%?}/extra3/bin:${bcbio_install_path%?}/extra_conda2/bin$:{PATH}
   echo " --- [$(date +"%F %R")] The PATH IS: ${PATH}"
   
   ln -s ${bcbio_install_path%?}/anaconda/envs/python2/bin/python ${bcbio_install_path%?}/tools/bin/python2

   ## install genomic data as described in the config file
   if [[ ${bcbio_annotated_species%?} = "yes" ]]; then
      # TODO check the arguments
      echo " --- [$(date +"%F %R")] Installing genome ${bcbio_genome%?} in the ${bcbio_install_path%?}/genomes path"
      bcbio_nextgen.py upgrade -u skip --genomes ${bcbio_genome%?} --datatarget variation --datatarget rnaseq --datatarget smallrna \
               --aligners bwa --aligners bowtie2 --aligners star --isolate --cores ${bcbio_total_cores%?}
   else
      ## install from custom genome, bcbio needs to have an already installed genome in order to install a custom genome
      echo " --- [$(date +"%F %R")] Installing Bcbio-nextgen with genome and transcriptome annotations for the sacCer3 reference"
      echo " --- [$(date +"%F %R")] This uses a small amount of storage and is required for creating a baseline Bcbio-nextgen configuration"
      
      ## install custom genome
      bcbio_nexgen.py upgrade -u skip --genomes sacCer3 --datatarget variation --datatarget rnaseq --datatarget smallrna \
               --aligners bwa --aligners bowtie2 --aligners hisat2 --aligners star --isolate --cores ${bcbio_total_cores} --mamba
      
      echo " --- [$(date +"%F %R")] Installing custom genome and transcriptome annotations from the user-provided FASTA and GTF files for the ${bcbio_genome} reference"
      
      bcbio_setup_genome.py -f ${bcbio_genome_fasta%?} -g ${bcbio_transcriptome_gtf%?} -i bwa bowtie2 star seq \
               -n ${bcbio_species%?} -b ${bcbio_genome%?} -c ${bcbio_total_cores%?} --buildversion ${bcbio_genome%?}
   fi
fi

##########################################################################################################################################################################################
                                                                         # VEP GENOME INSTALLATION #
##########################################################################################################################################################################################

## Install VEP cache data, if required
if [[ ${bcbio_annotated_species%?} = "yes" ]]; then
   if [[ ! -d ${bcbio_install_path%?}/genomes/${bcbio_species%?}/${bcbio_genome%?}/vep/${bcbio_vep_species%?}/${bcbio_ensembl_ver%?}_${bcbio_vep_assembly%?} ]]; then
      ## Configure VEP --- download cached data for the relevant species
      vep_install -s ${bcbio_vep_species%?} --NO_HTSLIB -a c -c ${bcbio_install_path%?}/genomes/${bcbio_species%?}/${bcbio_genome%?}/vep \
         --NO_UPDATE --VERSION ${bcbio_ensembl_ver%?} --ASSEMBLY ${bcbio_vep_assembly%?}
   fi
fi


echo " --- [$(date +"%F %R")] Creating analysis environment for ${bcbio_workflow%?} workflow analysis"

## Cleanup install scripts
rm ${path_to_scripts}/bcbio_nextgen_install.py*
rm ${path_to_scripts}/Miniconda3-*.sh*
rm ${path_to_scripts}/Miniconda2-*.sh*


## GO TO SAMPLES MODULE
## Handle the samples for the upcoming analysis

bash ${path_to_scripts}/samples_module.sh
