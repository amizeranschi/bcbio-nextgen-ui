#!/bin/bash

## Give excecution rights to every script
chmod +x *.sh

# CREATE MAIN WORKSPACE, VARIABLES AND SET PATHS
## Parse YAML configuration file to get all user's analsis data
source parse_yaml.sh $1
## Set environment variables, links for installers and paths
source set_environment_variables.sh

## Call the module for the setup of  bcbio_nextgen and the environment
## either install tool or create symlinks to already installed version
## Create Python environments to install the following tools
##                      * bioconda 
##                      * r-base 
##                      * r-xml 
##                      * bedops 
##                      * vcftools
##                      * sra-tools 
##                      * perl-net-ssleay 
##                      * gitentrez-direct
##                      * tassel
##                      * beagle
##                      * plink 
##                      * openssl=1.1.1l 
##                      * ensembl-vep=105
##                      * bcftools>=1.13
##                      * bcbio-monitor
##                      * pytz
##                      * python-dateutil
##                      * faststructure

bash setup_environment_module.sh
# bash samples_module.sh
