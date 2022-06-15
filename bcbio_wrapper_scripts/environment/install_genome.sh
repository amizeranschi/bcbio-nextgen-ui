#!/bin/bash

##########################################################################################################################################################################################
                                                                         # INSTALL GENOME #
##########################################################################################################################################################################################

## create symlink for python2 in bcbio in order to avoid any errors
ln -s ${bcbio_install_path}/anaconda/envs/python2/bin/python ${bcbio_install_path}/tools/bin/python2

## install genomic data as described in the config file
if [[ ${bcbio_annotated_species} = "yes" ]]; then
    # TODO check the arguments
    echo " --- [$(date +"%F %R")] Installing genome ${bcbio_genome} in the ${bcbio_install_path}/genomes path"
    bcbio_nextgen.py upgrade -u skip --genomes ${bcbio_genome} --datatarget variation --datatarget rnaseq --datatarget smallrna \
            --aligners bwa --aligners bowtie2 --aligners star --isolate --cores ${bcbio_total_cores}
else
    ## install from custom genome, bcbio needs to have an already installed genome in order to install a custom genome
    echo " --- [$(date +"%F %R")] Installing Bcbio-nextgen with genome and transcriptome annotations for the sacCer3 reference"
    echo " --- [$(date +"%F %R")] This uses a small amount of storage and is required for creating a baseline Bcbio-nextgen configuration"
    
    ## install custom genome
    bcbio_nextgen.py upgrade -u skip --genomes sacCer3 --datatarget variation --datatarget rnaseq --datatarget smallrna \
            --aligners bwa --aligners bowtie2 --aligners hisat2 --aligners star --isolate --cores ${bcbio_total_cores} --mamba
    
    echo " --- [$(date +"%F %R")] Installing custom genome and transcriptome annotations from the user-provided FASTA and GTF files for the ${bcbio_genome} reference"
    
    bcbio_setup_genome.py -f ${bcbio_genome_fasta} -g ${bcbio_transcriptome_gtf} -i bwa bowtie2 star seq \
            -n ${bcbio_species} -b ${bcbio_genome} -c ${bcbio_total_cores} --buildversion ${bcbio_genome}
fi
