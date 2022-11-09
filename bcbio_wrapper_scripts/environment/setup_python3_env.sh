#!/bin/bash

##########################################################################################################################################################################################
                                                                        # PYTHON 3 ENVIRONMENT #
##########################################################################################################################################################################################

## create a separate Python3 environment for aditional packages

if [[ ! -d ${bcbio_install_path}/extra3 ]]; then
   echo "--- [$(date +"%F %R")] Setting up a Python3 environment for utility packages"

   bash ${path_to_scripts}/Miniconda3-*.sh -b -p ${bcbio_install_path}/extra3
   ln -s ${bcbio_install_path}/extra3/bin/conda ${bcbio_install_path}/extra3/bin/conda_extra3

   ${bcbio_install_path}/extra3/bin/conda_extra3 install --yes -c conda-forge -c bioconda "mamba>=0.27"
   
   ln -s ${bcbio_install_path}/extra3/bin/mamba ${bcbio_install_path}/extra3/bin/mamba_extra3


   ## install the required R packages


   ## try to install all available CRAN / Bioconductor packages from conda
#   ${bcbio_install_path}/extra3/bin/mamba_extra3 install --yes -c conda-forge -c bioconda gxx_linux-64 wget git bedops vcftools pyyaml perl-yaml perl-encode-locale "sra-tools>=2.11" perl-net-ssleay entrez-direct tassel beagle gcta \
#   "gdal>=2.0.1" plink "bcftools>=1.15" haploview kraken2 blast "ensembl-vep>=105" "perl-bioperl>=1.7.2" "r-base>=4.1" r-essentials r-xml r-knitr r-markdown r-rmarkdown r-shiny r-rcpp r-devtools r-reticulate \
#   perl-bio-tools-run-alignment-tcoffee r-terra r-rgdal r-ragg r-ade4 r-aghmatrix bioconductor-annotationdbi bioconductor-annotationhub r-bglr r-candisc r-car bioconductor-chipseeker bioconductor-clusterprofiler \
#   r-cmplot r-corrgram r-corrplot r-data.table bioconductor-deseq2 bioconductor-do.db r-doparallel bioconductor-dose r-dplyr r-dt r-emmeans bioconductor-enrichplot bioconductor-ensembldb bioconductor-ensemblvep \
#   bioconductor-gdsfmt bioconductor-genesis r-ggally r-ggplot2 r-ggrepel r-ggroups bioconductor-gosemsim r-gparotation r-gtools r-hierfstat r-lsmeans bioconductor-meshdbi bioconductor-meshes r-misctools \
#   r-multcomp r-nnet r-orthopolynom r-outliers bioconductor-pcamethods r-pheatmap r-plotly r-psych r-qqman bioconductor-qvalue r-raster r-rcolorbrewer bioconductor-reactomepa r-readxl bioconductor-regparallel \
#   r-reshape2 r-rgdal r-rlist r-rrblup r-rstatix bioconductor-snprelate bioconductor-snpstats r-sp r-stringr r-tidyr r-tidyverse bioconductor-variantannotation r-vcfr r-vegan r-writexl 


   ## removed r-aghmatrix and r-ggroups because of dependency issues (as of 2022-11-09)
#   ${bcbio_install_path}/extra3/bin/mamba_extra3 install --yes -c conda-forge -c bioconda gxx_linux-64 wget git bedops vcftools pyyaml perl-yaml perl-encode-locale "sra-tools>=2.11" perl-net-ssleay entrez-direct tassel beagle gcta \
#   "gdal>=2.0.1" plink "bcftools>=1.15" haploview kraken2 blast "ensembl-vep>=105" "perl-bioperl>=1.7.2" "r-base>=4.1" r-essentials r-xml r-knitr r-markdown r-rmarkdown r-shiny r-rcpp r-devtools r-reticulate \
#   perl-bio-tools-run-alignment-tcoffee r-terra r-rgdal r-ragg r-ade4 bioconductor-annotationdbi bioconductor-annotationhub r-bglr r-candisc r-car bioconductor-chipseeker bioconductor-clusterprofiler \
#   r-cmplot r-corrgram r-corrplot r-data.table bioconductor-deseq2 bioconductor-do.db r-doparallel bioconductor-dose r-dplyr r-dt r-emmeans bioconductor-enrichplot bioconductor-ensembldb bioconductor-ensemblvep \
#   bioconductor-gdsfmt bioconductor-genesis r-ggally r-ggplot2 r-ggrepel bioconductor-gosemsim r-gparotation r-gtools r-hierfstat r-lsmeans bioconductor-meshdbi bioconductor-meshes r-misctools \
#   r-multcomp r-nnet r-orthopolynom r-outliers bioconductor-pcamethods r-pheatmap r-plotly r-psych r-qqman bioconductor-qvalue r-raster r-rcolorbrewer bioconductor-reactomepa r-readxl bioconductor-regparallel \
#   r-reshape2 r-rgdal r-rlist r-rrblup r-rstatix bioconductor-snprelate bioconductor-snpstats r-sp r-stringr r-tidyr r-tidyverse bioconductor-variantannotation r-vcfr r-vegan r-writexl 



   ## only install some important CRAN packages from conda; the rest of the needed R packages will be installed later from within R itself
   ${bcbio_install_path}/extra3/bin/mamba_extra3 install --yes -c conda-forge -c bioconda gxx_linux-64 wget git bedops vcftools pyyaml perl-yaml perl-encode-locale "sra-tools>=2.11" perl-net-ssleay entrez-direct \
   tassel beagle gcta "gdal>=2.0.1" plink "bcftools>=1.15" haploview gatk4 kraken2 blast "ensembl-vep>=107" "perl-bioperl>=1.7.2" "r-base>=4.2" r-essentials r-raster r-xml r-knitr r-rmarkdown r-shiny r-rcpp \
   r-devtools r-reticulate perl-bio-tools-run-alignment-tcoffee r-terra r-rgdal r-ragg
   
   
   
   ## install the Encode::Locale perl module
   yes | ${bcbio_install_path}/extra3/bin/cpan install Encode::Locale
   
   
   ## install R packages 
   echo " --- [$(date +"%F %R")] Installing required R packages"
   ${bcbio_install_path}/extra3/bin/Rscript --vanilla ${path_to_scripts}/downstreamAnalysis/install_R_packages.R
   
fi
