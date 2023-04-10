#!/bin/bash

## create a separate Python3 environment for aditional packages

if [[ ! -d ${bcbio_install_path:?}/extra ]]; then
   echo "--- [$(date +"%F %R")] Setting up custom Python3 and Python2 environments for utility packages"
   
   bash ${path_to_scripts}/Miniconda3-*.sh -b -p ${bcbio_install_path:?}/extra
   ln -s ${bcbio_install_path:?}/extra/bin/conda ${bcbio_install_path:?}/extra/bin/extra_conda
   
   ${bcbio_install_path:?}/extra/bin/extra_conda install --yes -c conda-forge -c bioconda "mamba>=0.27"
   
   ln -s ${bcbio_install_path:?}/extra/bin/mamba ${bcbio_install_path:?}/extra/bin/extra_mamba
   
   ## create python2 and python3 environments
   
   ${bcbio_install_path:?}/extra/bin/extra_mamba create --name py3 python=3.10 -y
   ${bcbio_install_path:?}/extra/bin/extra_mamba create --name py2 python=2.7.18 -y
   
   echo "--- [$(date +"%F %R")] Installing packages in the Python2 environment"
   
   ${bcbio_install_path:?}/extra/bin/extra_mamba install --name py2 --yes -c conda-forge -c bioconda faststructure pathoscope
   ${bcbio_install_path:?}/extra/envs/py2/bin/pip install bcbio-monitor pytz python-dateutil
   
   echo "--- [$(date +"%F %R")] Installing packages in the Python3 environment"
   
   ${bcbio_install_path:?}/extra/bin/extra_mamba install --name py3 --yes -c conda-forge -c bioconda kraken2 krakenuniq blast spades masurca quast \
   flye shasta repeatmodeler repeatmasker gxx_linux-64 wget git bedops vcftools pyyaml perl-yaml tectonic gnuplot \
   perl-encode-locale "sra-tools>=2.11" perl-net-ssleay entrez-direct tassel beagle gcta "gdal>=2.0.1" plink "bcftools>=1.15" haploview \
   "ensembl-vep>=109" "perl-bioperl>=1.7.2" "r-base>=4.2" r-essentials r-xml r-knitr r-markdown r-rmarkdown r-shiny r-rcpp r-devtools r-reticulate \
   perl-bio-tools-run-alignment-tcoffee r-terra r-rgdal r-ragg r-ade4 r-aghmatrix bioconductor-ahmeshdbs bioconductor-annotationdbi bioconductor-annotationhub r-bedr \
   r-bglr r-candisc r-car bioconductor-chipseeker bioconductor-clusterprofiler r-cmplot r-corrgram r-corrplot r-data.table bioconductor-deseq2 bioconductor-do.db \
   r-doparallel bioconductor-dose r-dplyr r-dt r-emmeans bioconductor-enrichplot bioconductor-ensembldb bioconductor-gdsfmt bioconductor-geneoverlap \
   bioconductor-genesis r-ggally r-ggplot2 r-ggrepel r-ggroups r-ggupset bioconductor-gosemsim r-gparotation r-gtools r-hierfstat r-lsmeans bioconductor-meshdbi \
   bioconductor-meshes r-misctools r-multcomp r-nnet r-orthopolynom r-outliers bioconductor-pathview bioconductor-pcamethods r-pheatmap r-plotly r-psych r-qqman \
   bioconductor-qvalue r-raster r-rcolorbrewer bioconductor-reactomepa r-readxl r-reshape2 r-rgdal r-rlist r-rrblup r-rstatix bioconductor-sbgnview \
   bioconductor-snprelate bioconductor-snpstats r-sp r-stringr r-tidyr r-tidyverse bioconductor-variantannotation r-vcfr r-vegan r-writexl
   
   ## install the Encode::Locale perl module
   yes | ${bcbio_install_path:?}/extra/envs/py3/bin/cpan install Encode::Locale
   
   ## install R packages 
   echo " --- [$(date +"%F %R")] Installing required R packages"
   ${bcbio_install_path:?}/extra/envs/py3/bin/Rscript --vanilla ${path_to_scripts}/downstreamAnalysis/install_R_packages.R
   
fi
