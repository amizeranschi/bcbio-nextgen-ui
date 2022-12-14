#!/bin/bash

##########################################################################################################################################################################################
                                                                        # PYTHON 3 ENVIRONMENT #
##########################################################################################################################################################################################

## create a separate Python3 environment for aditional packages

if [[ ! -d ${bcbio_install_path}/extra ]]; then
   echo "--- [$(date +"%F %R")] Setting up Python3 and Python2 environments for utility packages"

   bash ${path_to_scripts}/Miniconda3-*.sh -b -p ${bcbio_install_path}/extra
   ln -s ${bcbio_install_path}/extra/bin/conda ${bcbio_install_path}/extra/bin/conda_extra

   ${bcbio_install_path}/extra/bin/conda_extra install --yes -c conda-forge -c bioconda "mamba>=0.27"
   
   ln -s ${bcbio_install_path}/extra/bin/mamba ${bcbio_install_path}/extra/bin/mamba_extra
   
   ## create python2 and python3 environments
   
   mamba_extra create --name py3 python=3.9 -y
   mamba_extra create --name py2 python=2.7.18 -y
   

   ${bcbio_install_path}/extra/envs/py2/bin/pip install bcbio-monitor pytz python-dateutil
   ${bcbio_install_path}/extra/bin/extra_mamba install --name py2 --yes -c conda-forge -c bioconda faststructure pathoscope

   ## install the required R packages


   ## removed r-aghmatrix and r-ggroups because of dependency issues (as of 2022-11-16) -- these will be installed later from Bioconductor
   
   ${bcbio_install_path}/extra/bin/mamba_extra install --name py3 --yes -c conda-forge -c bioconda gxx_linux-64 wget git bedops vcftools pyyaml perl-yaml \
   perl-encode-locale "sra-tools>=2.11" perl-net-ssleay entrez-direct tassel beagle gcta "gdal>=2.0.1" plink "bcftools>=1.15" haploview kraken2 blast \
   "ensembl-vep>=105" "perl-bioperl>=1.7.2" "r-base>=4.1" r-essentials r-xml r-knitr r-markdown r-rmarkdown r-shiny r-rcpp r-devtools r-reticulate \
   perl-bio-tools-run-alignment-tcoffee r-terra r-rgdal r-ragg r-ade4 bioconductor-annotationdbi bioconductor-annotationhub r-bglr r-candisc r-car \
   bioconductor-chipseeker bioconductor-clusterprofiler r-cmplot r-corrgram r-corrplot r-data.table bioconductor-deseq2 bioconductor-do.db r-doparallel \
   bioconductor-dose r-dplyr r-dt r-emmeans bioconductor-enrichplot bioconductor-ensembldb bioconductor-ensemblvep bioconductor-gdsfmt \
   bioconductor-genesis r-ggally r-ggplot2 r-ggrepel bioconductor-gosemsim r-gparotation r-gtools r-hierfstat r-lsmeans bioconductor-meshdbi \
   bioconductor-meshes r-misctools r-multcomp r-nnet r-orthopolynom r-outliers bioconductor-pcamethods r-pheatmap r-plotly r-psych r-qqman \
   bioconductor-qvalue r-raster r-rcolorbrewer bioconductor-reactomepa r-readxl bioconductor-regparallel r-reshape2 r-rgdal r-rlist r-rrblup r-rstatix \
   bioconductor-snprelate bioconductor-snpstats r-sp r-stringr r-tidyr r-tidyverse bioconductor-variantannotation r-vcfr r-vegan r-writexl 


   
   ## install the Encode::Locale perl module
   yes | ${bcbio_install_path}/extra/envs/py3/bin/cpan install Encode::Locale
   
   
   ## install R packages 
   echo " --- [$(date +"%F %R")] Installing required R packages"
   ${bcbio_install_path}/extra/envs/py3/bin/Rscript --vanilla ${path_to_scripts}/downstreamAnalysis/install_R_packages.R
   
fi
