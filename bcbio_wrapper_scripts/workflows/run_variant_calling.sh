#!/bin/bash

echo ""
echo " --- [$(date +"%F %R")] Starting the Variant Calling Workflow"

##########################################################################################################################################################################################
                                                                            # VARIANT CALLING WORKFLOW #
##########################################################################################################################################################################################

# Perform variant calling
## Perform variant calling if the output (final) VCF file does not already exist:

if [ ! -f ${bcbio_runs_final}/*_${bcbio_exp_name}/*-gatk-haplotype*.vcf.gz ]; then
    echo " --- [$(date +"%F %R")] Running variant calling"
    cd ${bcbio_runs_input}
    ## If there's no variant_regions.bed file provided, create one with the full chromosomes of the bcbio genome
    if [[ ! -f full_chr_${bcbio_genome:?}.bed ]]; then
        if [[ ! -f ${genome_dir}/seq/${bcbio_genome:?}.fa.fai ]]; then
            echo "--- [$(date +"%F %R")] Creating an index for the genome ${bcbio_genome:?}"
            faidx ${genome_dir}/seq/${bcbio_genome:?}.fa
        fi
        
        echo " --- [$(date +"%F %R")] Creating a BED file for the full chromosomes of genome ${bcbio_genome:?}"
        awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${genome_dir}/seq/${bcbio_genome:?}.fa.fai > full_chr_${bcbio_genome:?}.bed
    fi

    if [[ ! -f variant_regions.bed ]]; then
        cp full_chr_${bcbio_genome:?}.bed variant_regions.bed
    fi
    cd ${bcbio_runs_input}
    bcbio_nextgen.py -w template gatk-variant.yaml ${bcbio_exp_name}.csv *.gz

    # exclusion of low complexity regions
    if [[ ${bcbio_exclude_lcr} = "yes" ]]; then
        ## Create a BED file with repetitive sequences (microsatellites and short tandem repeats) taken from UCSC
        if [[ ! -f ${bcbio_install_path:?}/genomes/${bcbio_species:?}/${bcbio_genome:?}/seq/lcr-${bcbio_genome:?}.bed ]]; then
            echo " --- [$(date +"%F %R")] Preparing a file with low complexity regions for genome ${bcbio_genome:?}"
            ## Download annotation files with repetitive regions from UCSC (http://hgdownload.soe.ucsc.edu/goldenPath/${bcbio_genome:?}/database/):
            ## microsatellite file:
            wget http://hgdownload.soe.ucsc.edu/goldenPath/${bcbio_genome:?}/database/microsat.txt.gz
            gunzip -f microsat.txt.gz
            ## simple repeats file:
            wget http://hgdownload.soe.ucsc.edu/goldenPath/${bcbio_genome:?}/database/simpleRepeat.txt.gz
            gunzip -f simpleRepeat.txt.gz
            ## Extract relevant columns from these files into BED files and concatenate them

            cut -f2-4 microsat.txt > microsat.bed
            cut -f2-4 simpleRepeat.txt > simpleRepeat.bed

            ## Concatenate the files, sort and merge the files
            cat *.bed | sort-bed - | bedops --merge - > lcr-${bcbio_genome:?}.bed

            ## Remove unneeded files
            rm -f microsat.txt microsat.bed simpleRepeat.txt simpleRepeat.bed
            ## move the file with lcr into the bcbio genome directory
            mv lcr-${bcbio_genome:?}.bed ${bcbio_install_path:?}/genomes/${bcbio_species:?}/${bcbio_genome:?}/seq/
            mv lcr-${bcbio_genome:?}.bed ${bcbio_install_path:?}/genomes/${bcbio_species:?}/${bcbio_genome:?}/seq/lcr-sorted-uniq.bed
        fi
        
        ## Check to see if a variant_regions.bed file was provided; if not, use the full chromosome lengths
        if [[ -f variant_regions.bed ]]; then
            echo " --- [$(date +"%F %R")] Subtracting low-complexity regions from the file variant_regions.bed"
            ## Create a back-up of the file variant_regions.bed and restore it when bcbio is finished running
            mv variant_regions.bed variant_regions.bed.bak
            
            ## Subtract the repetitive sequences from the variant_regions.bed file
            #bedtools subtract -a variant_regions.bed -b ${bcbio_path:?}/genomes/${bcbio_species:?}/${bcbio_genome:?}/seq/lcr-${bcbio_genome:?}.bed > variant_regions.bed
            bedops --difference variant_regions.bed.bak ${bcbio_install_path:?}/genomes/${bcbio_species:?}/${bcbio_genome:?}/seq/lcr-${bcbio_genome:?}.bed > variant_regions.bed
            #ls -l lcr-${bcbio_genome:?}.bed
        else
            echo " --- [$(date +"%F %R")] Subtracting low-complexity regions from the full chromosomes of genome ${bcbio_genome:?}"
            ## Create a BED file from the genome FASTA file, listing the full chromosome lengths
            if [[ ! -f full_chr_${bcbio_genome:?}.bed ]]; then
                echo " --- [$(date +"%F %R")] Creating a BED file for the full chromosomes of genome ${bcbio_genome:?}"
                awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${bcbio_install_path:?}/genomes/${bcbio_species:?}/${bcbio_genome:?}/seq/${bcbio_genome:?}.fa.fai > full_chr_${bcbio_genome:?}.bed
            fi
            
            ## Subtract the repetitive sequences from the variant_regions.bed file
            bedops --difference full_chr_${bcbio_genome:?}.bed ${bcbio_install_path:?}/genomes/${bcbio_species:?}/${bcbio_genome:?}/seq/lcr-${bcbio_genome:?}.bed > variant_regions.bed
        fi	
    fi

    # copy the variant_regiond.bed file to the config directory
    cp ${bcbio_runs_input}/variant_regions.bed ${bcbio_workflow_config}

    # go to work dir
    cd ${bcbio_workflow_work}
    # TODO cluster usage or not

    echo " --- [$(date +"%F %R")] Running bcbio-nextgen locally, using ${bcbio_total_cores} CPU cores"
    bcbio_nextgen.py ${bcbio_workflow_config}/${bcbio_exp_name}.yaml -n ${bcbio_total_cores}

else
    echo " --- [$(date +"%F %R")] Skipping variant calling because the VCF file exists: "${bcbio_runs_final}/*_${bcbio_exp_name}/*-gatk-haplotype*.vcf.gz

fi

## clean work directory
rm -rf ${bcbio_workflow_work}

## copy multiqc report
cp ${bcbio_runs_final}/*${bcbio_exp_name}/multiqc/multiqc_report.html ${path_to_web}

## print message for workflow completed
echo " --- [$(date +"%F %R")] Variant calling workflow is finished."
