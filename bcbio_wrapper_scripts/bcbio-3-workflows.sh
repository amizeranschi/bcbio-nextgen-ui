#!/bin/bash

## Install SRA tools in a Conda environment

## SRA-tools main page: https://hpc.nih.gov/apps/sratoolkit.html
## set the path where sra-tools should be installed
SRA_PATH=${HOME}/sra-tools

mkdir ${SRA_PATH}
cd ${HOME}
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $SRA_PATH/anaconda
cd ${SRA_PATH}

## create symlinks for easy access
ln -s ${SRA_PATH}/anaconda/bin/conda ${SRA_PATH}/anaconda/bin/sra_conda
ln -s ${SRA_PATH}/anaconda/bin/python ${SRA_PATH}/anaconda/bin/sra_python

## add the path to sra-tools binaries to the system's $PATH
echo "export PATH=$SRA_PATH/anaconda/bin:\${PATH}" >> $HOME/.bashrc
source ${HOME}/.bashrc

## install Mamba (a faster replacement of Conda)
sra_conda install -c conda-forge -c bioconda mamba --yes
ln -s ${SRA_PATH}/anaconda/bin/mamba ${SRA_PATH}/anaconda/bin/sra_mamba

## install sra-tools
mamba install -c conda-forge -c bioconda sra-tools --yes





## install bcbio-nextgen
## online documentation: https://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html
cd ${HOME}
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
python3 bcbio_nextgen_install.py ${HOME}/bcbio-nextgen --tooldir=${HOME}/bcbio-nextgen/tools --nodata --mamba
## add the locations to the bcbio-nextgen dependency executables to the system's $PATH variable
echo "export PATH=${HOME}/bcbio-nextgen/anaconda/bin:${HOME}/bcbio-nextgen/tools/bin:\$PATH" >> ~/.bashrc
source ~/.bashrc
## upgrade bcbio to the latest development version
#bcbio_nextgen.py upgrade -u development
## download genome data for sacCer3
bcbio_nextgen.py upgrade -u skip --genomes sacCer3 --aligners bwa --aligners bowtie2 --aligners hisat2 --aligners star





## prepare FAIRE-seq and RNA-seq data



## FAIRE-seq
mkdir ${HOME}/bcbio-runs && cd ${HOME}/bcbio-runs
mkdir faire-seq
mkdir rna-seq
cd faire-seq
fasterq-dump --split-files -O . -t . SRR6059150 SRR6059151 SRR6783014 SRR6783015 SRR6783016 SRR6784354
mv SRR6059150.fastq 500-F-Rep3.fastq
gzip -c 500-F-Rep3.fastq > 500-F-Rep3.fastq.gz
mv SRR6059151.fastq 500-I-Rep3.fastq
gzip -c 500-I-Rep3.fastq > 500-I-Rep3.fastq.gz
mv SRR6783014.fastq 500-I-Rep1.fastq
gzip -c 500-I-Rep1.fastq > 500-I-Rep1.fastq.gz
mv SRR6783015.fastq 500-F-Rep2.fastq
gzip -c 500-F-Rep2.fastq > 500-F-Rep2.fastq.gz
mv SRR6783016.fastq 500-I-Rep2.fastq
gzip -c 500-I-Rep2.fastq > 500-I-Rep2.fastq.gz
mv SRR6784354.fastq 500-F-Rep1.fastq
gzip -c 500-F-Rep1.fastq > 500-F-Rep1.fastq.gz
rm *.fastq




## RNA-seq -- downloald FASTQ, map reads to genome, extract reads that map to chrII, convert BAM to FASTQ and remove original FASTQ
cd ${HOME}/bcbio-runs/rna-seq


## main for loop, iterating over the 6 SRA run IDs for the RNA-seq data (Aerial 1, 2, 3 and Root 1, 2, 3)

for SRA_run_ID in SRR5482575 SRR5482576 SRR5482577 SRR5482578 SRR5482579 SRR5482580
do
  ## download and compress reads
  fasterq-dump --split-files -O . -t . ${SRA_run_ID}
  mv ${SRA_run_ID}_1.fastq ${SRA_run_ID}-R1.fastq
  bgzip ${SRA_run_ID}-R1.fastq
  mv ${SRA_run_ID}_2.fastq ${SRA_run_ID}-R2.fastq
  bgzip ${SRA_run_ID}-R2.fastq
  ## map reads to the sacCer3 genome using hisat2 (a splice-aware read mapper is required for RNA-seq reads)
  ## see: http://daehwankimlab.github.io/hisat2/manual/
  ## convert the output of hisat2 to a sorted and indexed BAM file using samtools (hisat2 outputs SAM)
  ## see: https://www.biostars.org/p/377961/#377980
  hisat2 -x /home/user/bcbio/genomes/Scerevisiae/sacCer3/hisat2/sacCer3 -1 ${SRA_run_ID}-R1.fastq.gz -2 ${SRA_run_ID}-R2.fastq.gz | \
    tee >(samtools flagstat - > ${SRA_run_ID}.flagstat) | \
    samtools sort -O BAM | \
    tee ${SRA_run_ID}.bam | \
    samtools index - ${SRA_run_ID}.bam.bai
  ## extract reads from chrII
  samtools view -b -o ${SRA_run_ID}-chrII.bam ${SRA_run_ID}.bam chrII
  ## sort reads by read name and convert the BAM file to FASTQ using bedtools bamtofastq
  ## see: https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html
  samtools sort -n -o ${SRA_run_ID}-chrII.qsort.bam ${SRA_run_ID}-chrII.bam
  ## remove the unsorted BAM and the original BAM (with all the chromosomes), as well as its flagstat
  rm ${SRA_run_ID}.bam* ${SRA_run_ID}.flagstat ${SRA_run_ID}-chrII.bam
  ## run bedtools bamtofastq and redirect the stderr stream to /dev/null to suppress warnings about improperly paired reads
  bedtools bamtofastq -i ${SRA_run_ID}-chrII.qsort.bam \
                        -fq ${SRA_run_ID}-chrII-R1.fastq \
                        -fq2 ${SRA_run_ID}-chrII-R2.fastq 2>/dev/null
  rm ${SRA_run_ID}-chrII.qsort.bam
  ## compress the new FASTQ files with bgzip
  bgzip ${SRA_run_ID}-chrII-R1.fastq
  bgzip ${SRA_run_ID}-chrII-R2.fastq
  ## replace the original FASTQ files with the reduced ones
  mv ${SRA_run_ID}-chrII-R1.fastq.gz ${SRA_run_ID}-R1.fastq.gz
  mv ${SRA_run_ID}-chrII-R2.fastq.gz ${SRA_run_ID}-R2.fastq.gz
done



## rename the RNA-seq fastq.gz files
mv SRR5482575-R1.fastq.gz AE1-R1.fastq.gz
mv SRR5482575-R2.fastq.gz AE1-R2.fastq.gz
mv SRR5482576-R1.fastq.gz AE2-R1.fastq.gz
mv SRR5482576-R2.fastq.gz AE2-R2.fastq.gz
mv SRR5482577-R1.fastq.gz AE3-R1.fastq.gz
mv SRR5482577-R2.fastq.gz AE3-R2.fastq.gz
mv SRR5482578-R1.fastq.gz RT1-R1.fastq.gz
mv SRR5482578-R2.fastq.gz RT1-R2.fastq.gz
mv SRR5482579-R1.fastq.gz RT2-R1.fastq.gz
mv SRR5482579-R2.fastq.gz RT2-R2.fastq.gz
mv SRR5482580-R1.fastq.gz RT3-R1.fastq.gz
mv SRR5482580-R2.fastq.gz RT3-R2.fastq.gz





## configure and run a FAIRE-seq analysis

cd ${HOME}/bcbio-runs/faire-seq
## download a template YAML file, describing the analysis
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/templates/illumina-chipseq.yaml
## edit the settings in the illumina-chipseq.yaml file, to make the analysis work for our files:
sed -i 's/genome_build: mm10/genome_build: sacCer3/' illumina-chipseq.yaml
sed -i 's/options: \["--broad", "-q 0.05"\]/options: \["--broad", "-q 0.05", "--nomodel", "--extsize 200"\]/' illumina-chipseq.yaml
## create a CSV file describing the samples
echo "samplename,description,batch,phenotype" > faire-seq-analysis.csv
echo "500-F-Rep1.fastq.gz,Rep1-500-F,batch1,chip" >> faire-seq-analysis.csv 
echo "500-I-Rep1.fastq.gz,Rep1-500-I,batch1,input" >> faire-seq-analysis.csv 
echo "500-F-Rep2.fastq.gz,Rep2-500-F,batch2,chip" >> faire-seq-analysis.csv 
echo "500-I-Rep2.fastq.gz,Rep2-500-I,batch2,input" >> faire-seq-analysis.csv 
echo "500-F-Rep3.fastq.gz,Rep3-500-F,batch3,chip" >> faire-seq-analysis.csv 
echo "500-I-Rep3.fastq.gz,Rep3-500-I,batch3,input" >> faire-seq-analysis.csv 

## configure the analysis
bcbio_nextgen.py -w template illumina-chipseq.yaml faire-seq-analysis.csv *.gz

## run bcbio
cd ${HOME}/bcbio-runs/faire-seq/faire-seq-analysis/work
bcbio_nextgen.py ../config/faire-seq-analysis.yaml -n 8

## copy the debug log to the final directory
cp log/bcbio-nextgen-debug.log ../final/*_faire-seq-analysis/

## remove the work directory
cd ../final/
rm -rf ../work/




## configure and run an RNA-seq analysis

## upgrade bcbio (including its dependencies) to the latest development version
#bcbio_nextgen.py upgrade -u development --tools

cd ${HOME}/bcbio-runs/rna-seq
## download a template YAML file, describing the analysis
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/templates/illumina-rnaseq.yaml
## edit the settings in the illumina-rnaseq.yaml file, to make the analysis work for our files:
sed -i 's/genome_build: GRCh37/genome_build: sacCer3/' illumina-rnaseq.yaml
sed -i 's/aligner: star/aligner: hisat2/' illumina-rnaseq.yaml
sed -i '/aligner: hisat2/ a \      tools_on: bcbiornaseq' illumina-rnaseq.yaml
sed -i '/tools_on: bcbiornaseq/ a \      bcbiornaseq:' illumina-rnaseq.yaml
sed -i '/bcbiornaseq:/ a \        organism: saccharomyces cerevisiae' illumina-rnaseq.yaml
sed -i '/organism: saccharomyces cerevisiae/ a \        interesting_groups: panel' illumina-rnaseq.yaml
## create a CSV file describing the samples
echo "samplename,description,panel" > rna-seq-analysis.csv 
echo "AE1,AE1,AE" >> rna-seq-analysis.csv 
echo "AE2,AE2,AE" >> rna-seq-analysis.csv 
echo "AE3,AE3,AE" >> rna-seq-analysis.csv 
echo "RT1,RT1,RT" >> rna-seq-analysis.csv 
echo "RT2,RT2,RT" >> rna-seq-analysis.csv 
echo "RT3,RT3,RT" >> rna-seq-analysis.csv 

## configure the analysis
bcbio_nextgen.py -w template illumina-rnaseq.yaml rna-seq-analysis.csv *.gz

## run bcbio
cd ${HOME}/bcbio-runs/rna-seq/rna-seq-analysis/work
bcbio_nextgen.py ../config/rna-seq-analysis.yaml -n 8

## copy the debug log to the final directory
cp log/bcbio-nextgen-debug.log ../final/*_rna-seq-analysis/

## remove the work directory
cd ../final/
rm -rf ../work/





## configure and run a variant calling analysis

mkdir ${HOME}/bcbio-runs/variant-calling
cd ${HOME}/bcbio-runs/variant-calling

## create soft links (similar to shortcuts on Windows) to the FASTQ files
ln -s /home/user/bcbio-runs/faire-seq/500-I-Rep1.fastq.gz 500-I-Rep1.fastq.gz 
ln -s /home/user/bcbio-runs/faire-seq/500-I-Rep2.fastq.gz 500-I-Rep2.fastq.gz 
ln -s /home/user/bcbio-runs/faire-seq/500-I-Rep3.fastq.gz 500-I-Rep3.fastq.gz 

## download a template YAML file, describing the analysis
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/templates/gatk-variant.yaml
## edit the settings in the illumina-rnaseq.yaml file, to make the analysis work for our files:
sed -i 's/genome_build: hg38/genome_build: sacCer3/' gatk-variant.yaml
sed -i 's/aligner: bwa/aligner: bowtie2/' gatk-variant.yaml
sed -i 's/recalibrate: gatk/# recalibrate: gatk/' gatk-variant.yaml

## create a CSV file describing the samples
echo "samplename,description,batch" > varcall-analysis.csv 
echo "500-I-Rep1.fastq.gz,Rep1-500-I,500-I-allReps" >> varcall-analysis.csv 
echo "500-I-Rep2.fastq.gz,Rep2-500-I,500-I-allReps" >> varcall-analysis.csv 
echo "500-I-Rep3.fastq.gz,Rep3-500-I,500-I-allReps" >> varcall-analysis.csv 

## configure the analysis
bcbio_nextgen.py -w template gatk-variant.yaml varcall-analysis.csv *.gz

## run bcbio
cd ${HOME}/bcbio-runs/variant-calling/varcall-analysis/work
bcbio_nextgen.py ../config/varcall-analysis.yaml -n 8

## copy the debug log to the final directory
cp log/bcbio-nextgen-debug.log ../final/*_varcall-analysis/

## remove the work directory
cd ../final/
rm -rf ../work/


## split the multi-sample VCF file
cd ${HOME}/bcbio-runs/variant-calling/varcall-analysis/final/*_varcall-analysis/

for file in *.vcf.gz; do
  for sample in `bcftools query -l $file`; do
    bcftools view -c1 -Ov -s ${sample} -o ${sample}.vcf ${file}
  done
done


