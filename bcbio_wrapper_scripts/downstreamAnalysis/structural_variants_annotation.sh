#!/bin/bash

##########################################################################################################################################################################################
                                                                            # STRUCTURAL VARIANTS #
##########################################################################################################################################################################################

# process structural variants, similarly to small variants (with slight differences between scenarios involving bcbio_nextgen)

# only the Manta structural variant caller is currently supported

if [ `ls -1 ${bcbio_runs_final}/*/*-manta.vcf.gz 2>/dev/null | wc -l ` -gt 0 ]; then
    echo " --- [$(date +"%F %R")] Performing variant annotation on structural variants."

    ## copy all the manta vcf files and rename them according to their parent directory's name:
    cd ${bcbio_runs_final}
    ## iterate through all the directories (which are named according to the samples)
    ## check if Manta VCF files exist

    for dir_file in */*-manta.vcf.gz; do
        ## extract the directory name
        dir=$( dirname ${dir_file} )
        cp ${dir_file} ${path_downstream_analysis}/${dir}-manta.vcf.gz
        ## copy the indexes as well
        cp ${dir_file}.tbi ${path_downstream_analysis}/${dir}-manta.vcf.gz.tbi
    done
        
    cd ${path_downstream_analysis}
        
    ## create a list of VCF.gz files to be merged
    rm -f vcf_files.list
    for file_name in *-manta.vcf.gz; do
        echo ${file_name} >> vcf_files.list
    done
        
    ## merge the VCF files with bcftools (the VCF files must be gzipped and indexed):
    file_size=$(wc -c < vcf_files.list)
    if (( "$file_size" > "0" )); then
        bcftools merge -Ov --file-list vcf_files.list > ${action_name}-struct-var.vcf
        ## reorder the samples according to the sample-order file
        bcftools view -Ov --force-samples -S ${bcbio_runs_input}/${action_name}/${action_name}-samples.txt ${action_name}-struct-var.vcf > ${action_name}-struct-var-sorted.vcf
        mv ${action_name}-struct-var-sorted.vcf ${action_name}-struct-var.vcf
    fi
        
    ## remove the single-sample VCF files and indexes
    for file_name in *-manta.vcf.gz*; do
        rm -f ${file_name}
    done
    rm -f vcf_files.list
fi

## extract the files
for vcf_file in *.vcf.gz; do
    if [ -f "${vcf_file}" ]; then
        gunzip -f ${vcf_file}
    fi
done


## process structural variants if these have been (successfully) called:
## check the file size (in bytes) of the file testingVC-struct-var.vcf
if [ -f "${action_name}-struct-var.vcf" ]; then
    file_size=$(wc -c < ${action_name}-struct-var.vcf)
else
    file_size=0
fi

if (( "$file_size" > "0" )); then
    vcf_file=${action_name}-struct-var.vcf
    ## extract the file name without the extension
    vcf_file_name=$(echo "${vcf_file}" | cut -f 1 -d '.')
    ## filter out variants that do not PASS the filters
    bcftools view -Ov -f .,PASS ${vcf_file_name}.vcf > ${vcf_file_name}-fil.vcf
    mv ${vcf_file_name}-fil.vcf ${vcf_file_name}.vcf
    bcftools view -Ov --force-samples -S ${bcbio_runs_input}/${action_name}/${action_name}-samples.txt ${vcf_file_name}.vcf > ${vcf_file_name}-ordered.vcf
    rm -f ${vcf_file_name}.vcf
    mv ${vcf_file_name}-ordered.vcf ${vcf_file_name}.vcf
    ## create stats for each VCF file
    bcftools stats ${vcf_file_name}.vcf > ${vcf_file_name}-bcf.stats.txt
    
    if [[ ${bcbio_annotated_species} = "yes" ]]; then
        ## annotate variants with Ensembl VEP
        vep --fork 4 --vcf --biotype --check_existing --distance 5000 --species ${bcbio_vep_species} --symbol --cache --dir_cache ${genome_dir}/vep --input_file ${vcf_file_name}.vcf --output_file ${vcf_file_name}-vep.vcf --force_overwrite --stats_file ${vcf_file_name}-vep.stats --stats_text
        mv ${vcf_file_name}-vep.vcf ${vcf_file_name}.vcf
        ## run VEP again and output a tab-separated file with only the most severe consequence per variant, skipping stats generation this time
        vep --fork 4 --tab --pick --no_stats --biotype --check_existing --distance 5000 --species ${bcbio_vep_species} --cache --dir_cache ${genome_dir}/vep --input_file ${vcf_file_name}.vcf --output_file ${vcf_file_name}-vep.table --force_overwrite
    else
        ## check if there exists a bgzipped and tabix-indexed GTF file for the genome
        if [ -f "${genome_dir}/rnaseq/ref-transcripts.gtf" ]; then
            if [ -f "${genome_dir}/rnaseq/ref-transcripts.gtf" ]; then
                ## check if it is compressed with bgzip; otherwise, remove and repackage
                if [[ $(grabix check ${genome_dir}/rnaseq/ref-transcripts.gtf) = "no" ]]; then
                    rm ${genome_dir}/rnaseq/ref-transcripts.gtf
                    ## see: https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html
                    grep -v "#" ${genome_dir}/rnaseq/ref-transcripts.gtf | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -@ ${bcbio_main_cores} -c > ${genome_dir}/rnaseq/ref-transcripts.gtf
                    tabix -p gff ${genome_dir}/rnaseq/ref-transcripts.gtf
                fi
            else
                ## see: https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html
                grep -v "#" ${genome_dir}/rnaseq/ref-transcripts.gtf | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -@ ${bcbio_main_cores} -c > ${genome_dir}/rnaseq/ref-transcripts.gtf
                tabix -p gff ${genome_dir}/rnaseq/ref-transcripts.gtf
            fi
        fi
        ## remove BND variants because VEP can't handle them (as of June 2020): https://github.com/Ensembl/ensembl-vep/issues/782#issuecomment-644031020
        sed "/SVTYPE=BND/d" ${vcf_file_name}.vcf > ${vcf_file_name}-no-BND.vcf
        ## run VEP with custom annotations from the GTF file
        vep --fork 4 --vcf --biotype --check_existing --distance 5000 --symbol --fasta ${genome_dir}/seq/${bcbio_genome}.fa --custom ${genome_dir}/rnaseq/ref-transcripts.gtf,ref-transcripts,gtf,overlap,0 --input_file ${vcf_file_name}-no-BND.vcf --output_file ${vcf_file_name}-no-BND-vep.vcf --force_overwrite --stats_file ${vcf_file_name}-no-BND-vep.stats --stats_text --max_sv_size 1000000000
        ## run VEP again and output a tab-separated file with only the most severe consequence per variant, skipping stats generation this time
        vep --fork 4 --tab --pick --no_stats --biotype --check_existing --distance 5000 --fasta ${genome_dir}/seq/${bcbio_genome}.fa --custom ${genome_dir}/rnaseq/ref-transcripts.gtf,ref-transcripts,gtf,overlap,0 --input_file ${vcf_file_name}-no-BND.vcf --output_file ${vcf_file_name}-no-BND-vep.table --force_overwrite --max_sv_size 1000000000
    fi
    ## extract relevant columns (GT and PR,SR) from a VCF file into tab-separated files
    gatk VariantsToTable -R ${genome_dir}/seq/${bcbio_genome}.fa -V ${vcf_file_name}.vcf -F CHROM -F POS -F ID -F REF -F ALT -GF GT -O ${vcf_file_name}-GT.table
else
    echo " --- [$(date +"%F %R")] Structural variants have not been called, or there was an error while calling them."
    rm -f ${action_name}-struct-var.vcf
fi
