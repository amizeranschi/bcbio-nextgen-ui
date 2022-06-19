# bcbio_nextgen_usability_improvements
Usability improvements for the Bcbio-nextgen data analysis pipeline: Thesis project on improving usability of the bcbio nextgen pipeline on various analysis workflows

The structure of the wrapper:

├── bcbio_wrapper_scripts <br />
├── deploy.sh <br />
├── downstreamAnalysis <br />
│   ├── bulk_rna_seq-downstream_analysis.R <br />
│   ├── metadata.csv <br />
│   ├── tryDownstreamAnalysis <br />
│   └── tximport-counts.csv <br />
├── downstreamAnalysisVariantCalling <br />
│   ├── gene_annotation_variant_calling.R <br />
│   ├── set_packages.R <br />
│   ├── small_variants_annotation.sh <br />
│   ├── structural_variants_annotation.sh <br />
│   └── variant_annotation.sh <br />
├── envrionment <br />
│   ├── install_bcbio_nextgen.sh<br />
│   ├── install_genome.sh <br />
│   ├── parse_yaml.sh <br />
│   ├── set_environment_variables.sh <br />
│   ├── setup_environment_module.sh <br />
│   ├── setup_python2_env.sh <br />
│   └── setup_python3_env.sh <br />
├── install_dependencies_interface.sh <br />
├── main.py <br />
├── result.yaml <br />
├── utils <br />
│   └── add_to_yaml.py <br />
├── web <br />
│   ├── about.html <br />
│   ├── data <br />
│   │   └── report_context.json <br />
│   ├── downstream_report.html <br />
│   ├── help.html <br />
│   ├── hero <br />
│   │   ├── banner.png <br />
│   │   └── favicon.ico <br /> 
│   ├── images <br /> 
│   │   ├── gridspec_ex.webp <br /> 
│   │   ├── plot_2.png <br />
│   │   └── plot_3.png <br />
│   ├── index.html <br /> 
│   ├── multiqc_report.html <br />
│   ├── run_config.html <br />
│   ├── script.js <br />
│   └── style.css <br />
├── workflows <br /> 
│   ├── config_module.sh <br />
│   ├── run_atac_seq.sh <br />
│   ├── run_bulk_rna_seq.sh <br />
│   ├── run_variant_calling.sh <br />
│   └── samples_module.sh <br />
└── yaml_to_table.py <br />



To run an anallysis it is required to configure the yaml configuration file infos.yaml.

 ### ***** CONFIGURATION STEPS *****

* existing_version:       ----> choose to install bcbio from scratch or from an existing install <br />
* path_to_existing:       ----> the usage of an older installation of bcbio requires the path to the install <br />
* development_branch:     ----> choose the development branch to upgrade and install bcbio to <br />
* total_cores:            ----> number of cores to run bcbio with <br />
* main_cores:             ----> number of cores to run bcbio with <br />
* install_path:           ----> the isnatllation path must be specified
                                WARNING: the path shall be located in the home directory of the system <br />
                                for the packages installed, paths that exceed 80 characters can not be processed <br />
* upgrade:                ----> choose to upgrade bcbio_nextgen or not <br />
* annotated_species:      ----> choose if the analysis will run on an existing genome in bcbio or a custom genome <br />
* genome_fasta:           ----> the path toward the .fa file of the custom genome <br />
                              sort gtf if no annotated species <br />
* transcriptome_gtf:      ----> the path to the transcriptome of the custom genome <br />
* species:                ----> the annotated species <br />
* genome:                 ----> the annotated genome <br />
* vep_species:            ----> species for usage of vep tool <br />
* vep_assembly:           ----> genome for usage of vep tool <br />
* ensembl_ver:            ----> vep tool version <br />
* workflow:               ----> name of the workflow <br />
                              convention available: <br />
                                        ** variant_calling for Variant calling and variant annotation <br />
                                        ** atac_seq for ATAC seq workflow <br />
* variant_annotation:     ----> when running variant calling workflow, there is the choice of running variant annotation also <br />
* exclude_lcr:            ----> when running variant calling workflow, there is the choice of performing exclusion of low complexity regions <br />
* download_samples:       ----> choose if download samples or get them from local system <br />
* path_to_samples_on_sys: ----> the choice to use locally stored samples requires the path to the samples <br />
* samples:                ----> id of the samples to download <br />
* samples_fastq:          ----> name convention for each sample id in order without extension <br />
                                if a sample has more than 1 file the names will be placed in order for _1 _2 or _3 <br />
                                if the samples are already on the system, write their names without the extension <br />
* csv_file_path:          ----> the path toward the csv file for the analysis <br />


RUN IN A SHELL LIKE THIS

    $ bash deploy.sh <your_yaml_configuration_file>.yaml
