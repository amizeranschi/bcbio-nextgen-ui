rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
my_args <- head(args)

## Clear all objects from the memory:
# Set working directory
work_dir=my_args[1]
# work_dir="/home/maria/bcbio_nextgen_usability_improvements/bcbio_wrapper_scripts/downstreamAnalysisBulk-RNA-seq"
setwd(work_dir)

# set files with genomic data
peaks_file=my_args[2]

# species name

# vep_species stores the original string provided by the user in the config file
vep_species = my_args[3]
# vep_species="Saccharomyces_cerevisiae"
# replace the underscore with a space to get the name of the organism
my_species = sub("_", " ", vep_species)

# get gtf file location
gtf_location = my_args[4]

path_to_scripts = my_args[5]

print(" --- Preparing data...")
peaks = ChIPseeker::readPeakFile(peaks_file, header = F)

## check the object type
class(peaks)

## print the first 6 entries
head(peaks)

## check the structure of the peaks object (its attributes and their types)
str(peaks)

## create a coverage plot using the R package 
ChIPseeker::covplot(peaks, weightCol="V5")

## annotate peaks
peakAnno = ChIPseeker::annotatePeak(peaks, tssRegion = c(-3000, 3000), TxDb = txdb)
peakAnno

## visualize genomic annotations of peaks
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)


## create a data frame with peak annotations and write it to a file
peakAnnoDF = data.frame(peakAnno)
write.table(peakAnnoDF, file = "peakAnnotations.txt", sep = "\t", row.names = F, quote = F)

## check the first 6 entries
head(peakAnnoDF)

## get the column names from peakAnnoDF
colnames(peakAnnoDF)

## check the structure of the peakAnnoDF data frame (its columns and their types)
str(peakAnnoDF)

## what was the total number of peaks?
val1 = nrow(peakAnnoDF)
val1

## how many peaks were annotated as "Promoter (<=1kb)"?
val2 = nrow( peakAnnoDF[ peakAnnoDF$annotation == "Promoter (<=1kb)", ] )
val2

## what percentage of peaks were annotated as "Promoter (<=1kb)"?
percentage = val2 * 100 / val1
percentage

## compare this with the original summary from peakAnno
peakAnno


## extract a vector of gene IDs (e.g. YAL069W, YAL068C etc.) and convert them to Entrez IDs using org.Sc.sgd.db
## see also: https://bioconductor.org/packages/release/workflows/vignettes/annotation/inst/doc/Annotation_Resources.html#orgdb-objects

## current gene IDs are of type "ORF"; see an example of "ORF" keys in org.Sc.sgd.db
head(keys(org.Sc.sgd.db, keytype=species_keytype, pattern = "YAL"))

## check some gene IDs from our annotationDF data frame
head(peakAnnoDF$geneId)

## get a list of Entrez IDs corresponding to the gene IDs from peakAnnoDF
entrezIDs = AnnotationDbi::select(org.Sc.sgd.db, keys = peakAnnoDF$geneId, keytype=species_keytype, columns = "ENTREZID")

workflow_name = "chip_seq"
system(paste(paste0("Rscript --vanilla ", path_to_scripts, "/downstreamAnalysis/computeMetrics.R"), entrezIDs, null, my_species, workflow_name))

# ## functional enrichment analysis using ReactomePA
# pathwayAnno = ReactomePA::enrichPathway(entrezIDs$ENTREZID, organism = "yeast")

# ## convert to pathwayAnno to data frame
# pathwayAnnoDF = data.frame(pathwayAnno)
# str(pathwayAnnoDF)

# ##remove the last 2 columns
# pathwayAnnoDF = pathwayAnnoDF[, -c(8, 9)]
# str(pathwayAnnoDF)

# nrow(pathwayAnnoDF)

# pathwayAnnoDF

# pathwayAnnoDF$Description

