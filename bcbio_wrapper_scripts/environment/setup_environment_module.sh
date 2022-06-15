#!/bin/bash

##########################################################################################################################################################################################
                                                                         # PYTHON ENVIRONMENTS #
##########################################################################################################################################################################################

if [[ ! -d ${bcbio_install_path%?} ]]; then
   mkdir ${bcbio_install_path%?}
fi

cd ${bcbio_install_path%?}

## install or check the Python2 and Python3 environments
bash ${path_to_scripts}/environment/setup_python2_env.sh
bash ${path_to_scripts}/environment/setup_python3_env.sh

##########################################################################################################################################################################################
                                                                         # ALREADY INSTALLED VERSION OF BCBIO #
##########################################################################################################################################################################################

## IF BCBIO IS ALREADY INSTALLED
if [[ ${bcbio_existing_version%?} = "yes" ]]; then

   echo " --- [$(date +"%F %R")] bcbio_nextgen already installedy installed on the system."
   echo " --- [$(date +"%F %R")] Use bcbio_nextgen already installed version."
   
   ## set the path with all the utils
   echo " --- [$(date +"%F %R")] Setting the PATH for Python 3 and Python2 environment installation."
   export PATH="${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:${bcbio_install_path%?}/extra3/bin:${bcbio_install_path%?}/extra2/bin${PATH:+:${PATH}}"
   echo " --- [$(date +"%F %R")] The PATH IS: ${PATH}"

   ## check if a genome is installed on the system
   ## if not, install it
   echo " --- [$(date +"%F %R")] Check if there is any genome installed on the system"
   if [[ ! -d ${bcbio_install_path%?}/genomes/${bcbio_species%?}/${bcbio_genome%?} ]]; then
      echo " --- [$(date +"%F %R")] The requested genome ${bcbio_genome%?} was not found."
      bash ${path_to_scripts}/install_genome.sh

   fi

   echo " --- [$(date +"%F %R")] The requested genome ${bcbio_genome%?} can be found at ${bcbio_install_path%?}."

   ## set biobambam version to be stable
   # ${bcbio_install_path%?}/anaconda/bin/bcbio_conda install -y biobambam=2.0.87 -c bioconda
fi

##########################################################################################################################################################################################
                                                                         # BCBIO INSTALLATION #
##########################################################################################################################################################################################

if [[ ${bcbio_existing_version%?} = "no" ]]; then
   echo " --- [$(date +"%F %R")] bcbio_nextgen not on the system."
   echo " --- [$(date +"%F %R")] bcbio_nextgen will start installation."

   ## install bcbio and the genome
   bash ${path_to_scripts}/environment/install_bcbio_nextgen.sh

   ## set the path with all the utils
   echo " --- [$(date +"%F %R")] Setting the PATH for Python 3 and Python2 environment installation."
   export PATH="${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:${bcbio_install_path%?}/extra3/bin:${bcbio_install_path%?}/extra2/bin${PATH:+:${PATH}}"
   echo " --- [$(date +"%F %R")] The PATH IS: ${PATH}"
   
   bash ${path_to_scripts}/environment/install_genome.sh
fi

##########################################################################################################################################################################################
                                                                         # VEP GENOME INSTALLATION #
##########################################################################################################################################################################################

## Install VEP cache data, if required
if [[ ${bcbio_annotated_species%?} = "yes" ]]; then
   if [[ ! -d ${bcbio_install_path%?}/genomes/${bcbio_species%?}/${bcbio_genome%?}/vep/${bcbio_vep_species%?}/${bcbio_ensembl_ver%?}_${bcbio_vep_assembly%?} ]]; then
      ## Configure VEP --- download cached data for the relevant species
      ${bcbio_install_path%?}/extra3/bin/vep_install -s ${bcbio_vep_species%?} --NO_HTSLIB -a c -c \
         ${bcbio_install_path%?}/genomes/${bcbio_species%?}/${bcbio_genome%?}/vep \
         --NO_UPDATE --VERSION ${bcbio_ensembl_ver%?} --ASSEMBLY ${bcbio_vep_assembly%?}
   fi
fi

echo " --- [$(date +"%F %R")] Creating analysis environment for ${bcbio_workflow%?} workflow analysis"

## Cleanup install scripts
rm ${path_to_scripts}/bcbio_nextgen_install.py*
rm ${path_to_scripts}/Miniconda3-*.sh*
rm ${path_to_scripts}/Miniconda2-*.sh*
