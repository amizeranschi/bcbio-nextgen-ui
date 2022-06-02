#!/bin/bash

##########################################################################################################################################################################################
                                                                         # BCBIO INSTALLATION #
##########################################################################################################################################################################################

## CREATE BCBIO ENVIRONMENT: install bcbio and genomes, if required
## Restriction: use a shorter path --> the home directory
## Avoid a padding error because of the packages built with 80 caracters padded prefix 

## install bcbio in bcbio_nexgen directory with no data
echo " --- [$(date +"%F %R")] Installing bcbio_nextgen into the directory: ${bcbio_install_path%?} with no data"
python3 ${path_to_scripts}/bcbio_nextgen_install.py ${bcbio_install_path%?} --tooldir=${bcbio_install_path%?}/tools --nodata --mamba

# export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:${PATH}

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

## rename conda to differentiate between different versions of the executable
ln -s ${bcbio_install_path%?}/anaconda/bin/conda ${bcbio_install_path%?}/anaconda/bin/bcbio_conda
ln -s ${bcbio_install_path%?}/anaconda/bin/python ${bcbio_install_path%?}/anaconda/bin/bcbio_python

## set biobambam version to be stable 
bcbio_conda install -y biobambam=2.0.87 -c bioconda
