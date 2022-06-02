#!/bin/bash

##########################################################################################################################################################################################
                                                                        # PYTHON 2 ENVIRONMENT #
##########################################################################################################################################################################################

if [[ ! -d ${bcbio_install_path%?}/extra2 ]]; then

   ## create a Python2 env for bcbio-monitor pytz and python-dateutil and faststructure

   echo " --- [$(date +"%F %R")] Setting up a Python2 environment for legacy utility packages"
   bash ${path_to_scripts}/Miniconda2*.sh -b -p ${bcbio_install_path%?}/extra2
   ln -s ${bcbio_install_path%?}/extra2/bin/conda ${bcbio_install_path%?}/extra2/bin/extra_conda2s

   ${bcbio_install_path%?}/extra2/bin/pip install bcbio-monitor pytz python-dateutil
   ${bcbio_install_path%?}/extra2/bin/extra_conda2 install --yes -c conda-forge -c bioconda faststructure

fi
