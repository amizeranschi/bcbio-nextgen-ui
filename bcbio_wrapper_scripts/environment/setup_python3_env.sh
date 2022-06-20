#!/bin/bash

##########################################################################################################################################################################################
                                                                        # PYTHON 3 ENVIRONMENT #
##########################################################################################################################################################################################

## create a separate Python3 environment for aditional packages

if [[ ! -d ${bcbio_install_path}/extra3 ]]; then
   echo "--- [$(date +"%F %R")] Setting up a Python3 environment for utility packages"

   bash ${path_to_scripts}/Miniconda3-*.sh -b -p ${bcbio_install_path}/extra3
   ln -s ${bcbio_install_path}/extra3/bin/conda ${bcbio_install_path}/extra3/bin/conda_extra3

   ${bcbio_install_path}/extra3/bin/conda_extra3 install --yes -c conda-forge -c bioconda mamba
   
   # r-rcolorbrewer 
   # r-pheatmap 
   # bioconductor-deseq2
   ln -s ${bcbio_install_path}/extra3/bin/mamba ${bcbio_install_path}/extra3/bin/mamba_extra3

   ## installing packages
   ${bcbio_install_path}/extra3/bin/mamba_extra3 install --yes -c conda-forge -c bioconda wget git bedops vcftools sra-tools pyyaml perl entrez-direct tassel beagle plink openssl=1.1.1l ensembl-vep=105 "bcftools>=1.13" biobambam=2.0.87
   ${bcbio_install_path}/extra3/bin/mamba_extra3 install --yes -c conda-forge r-base
   ${bcbio_install_path}/extra3/bin/mamba_extra3 install --yes -c conda-forge r-essentials
   ${bcbio_install_path}/extra3/bin/conda_extra3 install gxx_linux-64

   ${bcbio_install_path}/extra3/bin/conda_extra3 install --yes -c conda-forge -c bioconda bioconductor-annotationhub bioconductor-annotationdbi bioconductor-meshdbi r-ggplot2 bioconductor-clusterprofiler bioconductor-dose bioconductor-meshes bioconductor-chipseeker bioconductor-gosemsim bioconductor-reactomepa bioconductor-variantannotation
fi
