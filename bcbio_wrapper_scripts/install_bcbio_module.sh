#!/bin/bash

# CREATE MAIN WORKSPACE, VARIABLES AND SET PATHS
## POSSIBLY ADDING THEM TO THE .YAML FILE OR utils.sh

miniconda3="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
miniconda2="https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh"
install_script_bcbio="https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py"

chmod +x *.sh

# PARSE THE INPUT YAML CONFIG FILE

## Function for reading from a simple YAML file, adapted from: https://stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script/21189044#21189044
## Modified the printf line in the awk script code below to add the word "export ", such that it creates environment variables accessible from the rest of the analysis scripts
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("export %s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

echo "using configuration file: " $1


## read the config file and create the variables using the "bcbio_" prefix
eval $(parse_yaml $1 "bcbio_")

echo "" 
echo "## printing the whole contents of the file:"

parse_yaml $1 "bcbio_"

# echo "## set the \$PATH variable like this:"
# echo "export PATH=${bcbio_install_path:?}/tools/bin:${bcbio_install_path:?}/anaconda/bin:${bcbio_install_path:?}/anaconda/envs/bcbiovm/bin:${bcbio_install_path:?}/extra3/bin:${bcbio_install_path:?}/extra2/bin:${PATH}"
# echo "" 

# CHECK IF BCBIO IS ALREADY INSTALLED
echo "--- [$(date +"%F %R")] check if bcbio is already installed"
## IF DECIDE TO KEEP THE OLDER VERSION, CHECK SYMLINKS OR CREATE NEW ONES
if [[ ${bcbio_existing_version%?} = "no" ]]; then
   echo ${bcbio_install_path%?} 

   cd ${HOME}

   #  GET UTILITARIES AND SCRIPTS
   ## check if wget is available. if not, use curl instead

   if [ -x "$(command -v wget)" ]; then
      wget ${miniconda3}
      wget ${miniconda2}
      wget ${install_script_bcbio}
   else
      curl -L -C - -O ${miniconda3}
      curl -L -C - -O ${miniconda2}
      curl -L -C - -O ${install_script_bcbio}
   fi
   
   # CREATE BCBIO ENVIRONMENT: install bcbio and genomes, if required
   ## install bcbio in bcbio_nexgen directory with no data
   ## Restriction: use a shorter path yada yada

   python3 bcbio_nextgen_install.py ${bcbio_install_path%?} --tooldir=${bcbio_install_path%?}/tools --nodata --mamba

   ## create symlink for bcbio
   echo "export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:\$PATH" >> ~/.bashrc
   source ~/.bashrc

   ## upgrade bcbio if requested
   bcbio_nextgen.py upgrade -u ${bcbio_development_branch%?} --tools

   ## make sure to have a recent samtools in the main conda environment
   # mamba install -c bioconda "samtools=1.13" --yes

   ## install genomic data as described in the config file
   if [[ ${bcbio_annotated_species%?} = "yes" ]]; then
      # TODO check the arguments
      echo " --- [$(date +"%F %R")] Installing bcbio_nextgen into the directory: ${bcbio_install_path%?}"

      bcbio_nextgen.py upgrade -u skip --genomes ${bcbio_genome%?} --datatarget variation --datatarget rnaseq --datatarget smallrna \
			  --aligners bwa --aligners bowtie2 --aligners star --isolate --cores ${bcbio_total_cores%?}
   else
      # install from custom genome, bcbio needs to have an already installed genome in order to install a custom genome
      echo " --- [$(date +"%F %R")] Installing Bcbio-nextgen with genome and transcriptome annotations for the sacCer3 reference"
		echo " --- [$(date +"%F %R")] This uses a small amount of storage and is required for creating a baseline Bcbio-nextgen configuration"
      # install custom genome
      bcbio_nexgen.py upgrade -u skip --genomes sacCer3 --datatarget variation --datatarget rnaseq --datatarget smallrna \
			  --aligners bwa --aligners bowtie2 --aligners hisat2 --aligners star --isolate --cores ${bcbio_total_cores} --mamba
      echo " --- [$(date +"%F %R")] Installing custom genome and transcriptome annotations from the user-provided FASTA and GTF files for the ${bcbio_genome} reference"
		bcbio_setup_genome.py -f ${bcbio_genome_fasta%?} -g ${bcbio_transcriptome_gtf%?} -i bwa bowtie2 star seq \
			  -n ${bcbio_species%?} -b ${bcbio_genome%?} -c ${bcbio_total_cores%?} --buildversion ${bcbio_genome%?}
   fi

   ## rename conda to differentiate between different versions of the executable
	ln -s ${bcbio_install_path%?}/anaconda/bin/conda ${bcbio_install_path%?}/anaconda/bin/bcbio_conda
	ln -s ${bcbio_install_path%?}/anaconda/bin/python ${bcbio_install_path%?}/anaconda/bin/bcbio_python
fi

# CREATE ENVIRONMENT FOR OTHER UTILS

## create a separate Python3 environment for aditional packages for downstream
## analysis and other

# echo " --- [$(date +"%F %R")] Setting up a Python3 environment for utility packages"

# cd ${HOME}

# bash Miniconda3-latest-Linux-x86_64.sh -b -p ${bcbio_install_path%?}/extra3

# ln -s ${bcbio_install_path%?}/extra3/bin/conda ${bcbio_install_path%?}/extra3/bin/conda_extra3

# ${bcbio_install_path%?}/extra3/bin/conda_extra3 install --yes -c conda-forge -c bioconda mamba

# ln -s ${bcbio_install_path%?}/extra3/bin/mamba ${bcbio_install_path%?}/extra3/bin/mamba_extra3

# ## installing packages
# ${bcbio_install_path%?}/extra3/bin/mamba_extra3 install --yes -c conda-forge -c bioconda wget git bedops vcftools sra-tools perl-net-ssleay entrez-direct tassel beagle plink
# echo "export PATH=${bcbio_install_path%?}/extra3/bin:\$PATH" >> ~/.bashrc
# source ~/.bashrc
# # TODO create symlinks for some packages

# ## create a Python2 env for bcbio-monitor pytz and python-dateutil and faststructure
# echo " --- [$(date +"%F %R")] Setting up a Python2 environment for legacy utility packages"
# bash Miniconda2-latest-Linux-x86_64.sh -b -p ${bcbio_install_path%?}/extra2
# ln -s ${bcbio_install_path%?}/extra2/bin/conda ${bcbio_install_path%?}/extra2/bin/extra_conda2

# ${bcbio_install_path%?}/extra2/bin/extra_conda2 install --yes -c conda-forge -c bioconda mamba
# ln -s ${bcbio_install_path%?}/extra2/bin/mamba ${bcbio_install_path%?}/extra2/bin/mamba_extra2

# ${bcbio_install_path%?}/extra2/bin/pip install bcbio-monitor pytz python-dateutil
# ${bcbio_install_path%?}/extra2/bin/mamba_extra2 install --yes -c conda-forge -c bioconda faststructure

# TODO:
# ## edit a script for bcbio monitor to make it compatible with latest dependencies
# sed -i "s/from gevent.wsgi import WSGIServer/from gevent.pywsgi import WSGIServer/g" ${bcbio_path:?}/extra2/lib/python2.7/site-packages/bcbio_monitor/cli.py



## set up a temporary PATH variable to include the previously installed python; keep a backup of the old $PATH
# old_PATH=$PATH
# PATH=$PATH:${bcbio_install_path:?}/extra3/bin

# ## restore PATH
# PATH=$old_PATH

# GENERATE TREE: bcbio_runs/workflow_name_project/-> input, config, work, final
# echo " --- [$(date +"%F %R")] Creating analysis environment for ${bcbio_workflow%?} workflow analysis"

# bcbio_runs="${HOME}/bcbio_runs/"
# workflow_name="workflow_${bcbio_workflow}"
# bcbio_workflow="${bcbio_runs}${workflow_name%?}"
# bcbio_runs_input="${bcbio_runs}${workflow_name%?}/input"
# bcbio_runs_config="${bcbio_runs}${workflow_name%?}/config"
# bcbio_runs_work="${bcbio_runs}${workflow_name%?}/work"
# bcbio_runs_final="${bcbio_runs}${workflow_name%?}/final"

# mkdir ${bcbio_runs}
# mkdir ${bcbio_workflow}
# mkdir ${bcbio_runs_input}
# mkdir ${bcbio_runs_config}
# mkdir ${bcbio_runs_work}
# mkdir ${bcbio_runs_final}

# cd ${bcbio_workflow}

# SAMPLES MODULE
## if user wants to download data 
# if [[ ${bcbio_download_samples%?} == "yes" ]]; then
#    echo "${bcbio_samples%?}"
# fi
# TODO 5    
# TODO 6    GENERATE TREE FOR WORKFLOW WORKSPACE
# TODO 7
