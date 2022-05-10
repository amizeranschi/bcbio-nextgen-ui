# bcbio_nextgen_usability_improvements
Usability improvements for the Bcbio-nextgen data analysis pipeline: Thesis project on improving usability of the bcbio nextgen pipeline on various analysis workflows

The structure of the wrapper:

├── bcbio_wrapper_scripts
│  deploy.sh       
│  parse_yaml.sh     
│  set_environment_variables.sh
│  setup_environment_module.sh    
│  install_bcbio_nextgen.sh
│  samples_module.sh
│  add_to_yaml.py    
│  config_module.sh
│  run_atac_seq.sh  
│  run_variant_calling.sh
│  run_variant_calling.sh
│  small_variants_annotation.sh
│  structural_variants_annotation.sh
│  infos_atac.yaml  
│  infos.yaml  


To run an anallysis it is required to configure the yaml configuration file infos.yaml.

***** CONFIGURATION STEPS *****

existing_version:       ----> choose to install bcbio from scratch or from an existing install
path_to_existing:       ----> the usage of an older installation of bcbio requires the path to the install
development_branch:     ----> choose the development branch to upgrade and install bcbio to
total_cores:            ----> number of cores to run bcbio with
main_cores:             ----> number of cores to run bcbio with
install_path:           ----> the isnatllation path must be specified
                                WARNING: the path shall be located in the home directory of the system
                                for the packages installed, paths that exceed 80 characters can not be processed
upgrade:                ----> choose to upgrade bcbio_nextgen or not
annotated_species:      ----> choose if the analysis will run on an existing genome in bcbio or a custom genome
genome_fasta:           ----> the path toward the .fa file of the custom genome
                              sort gtf if no annotated species
transcriptome_gtf:      ----> the path to the transcriptome of the custom genome
species:                ----> the annotated species
genome:                 ----> the annotated genome
vep_species:            ----> species for usage of vep tool
vep_assembly:           ----> genome for usage of vep tool
ensembl_ver:            ----> vep tool version
workflow:               ----> name of the workflow
                              convention available:
                                        * variant_calling for Variant calling and variant annotation
                                        * atac_seq for ATAC seq workflow

variant_annotation:     ----> when running variant calling workflow, there is the choice of running variant annotation also

download_samples:       ----> choose if download samples or get them from local system
path_to_samples_on_sys: ----> the choice to use locally stored samples requires the path to the samples
samples:                ----> id of the samples to download
samples_fastq:          ----> name convention for each sample id in order without extension
                                if a sample has more than 1 file the names will be placed in order for _1 _2 or _3
                                if the samples are already on the system, write their names without the extension
csv_file_path:          ----> the path toward the csv file for the analysis


RUN IN A SHELL LIKE THIS

    $ bash deploy.sh <your_yaml_configuration_file>.yaml