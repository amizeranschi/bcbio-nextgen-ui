#!/bin/bash

##########################################################################################################################################################################################
                                                                         # BCBIO INSTALLATION #
##########################################################################################################################################################################################

## CREATE BCBIO ENVIRONMENT: install bcbio and genomes, if required
## Restriction: use a shorter path --> the home directory
## Avoid a padding error because of the packages built with 80 caracters padded prefix 

## install bcbio in bcbio_nexgen directory with no data
echo " --- [$(date +"%F %R")] Installing bcbio_nextgen into the directory: ${bcbio_install_path:?} with no data"
python3 ${path_to_scripts}/bcbio_nextgen_install.py ${bcbio_install_path:?} --tooldir=${bcbio_install_path:?}/tools --nodata --mamba

## rename conda to differentiate between different versions of the executable
ln -s ${bcbio_install_path:?}/anaconda/bin/conda ${bcbio_install_path:?}/anaconda/bin/bcbio_conda
ln -s ${bcbio_install_path:?}/anaconda/bin/python ${bcbio_install_path:?}/anaconda/bin/bcbio_python

## set biobambam version to be stable 
#${bcbio_install_path:?}/anaconda/bin/bcbio_conda install -y biobambam=2.0.87 -c bioconda
