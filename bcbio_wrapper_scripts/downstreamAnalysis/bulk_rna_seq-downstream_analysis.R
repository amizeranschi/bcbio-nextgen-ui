## Clear all objects from the memory:
rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
my_args <- head(args)

# Set working directory
work_dir=my_args[1]
# work_dir="/home/maria/bcbio_nextgen_usability_improvements/bcbio_wrapper_scripts/downstreamAnalysisBulk-RNA-seq"
setwd(work_dir)

# set files with genomic data
counts_file=my_args[2]
# counts_file="/home/maria/bcbio_nextgen_usability_improvements/test_data_bulk/tximport-counts.csv"
metadata_file=my_args[3]
# metadata_file="/home/maria/bcbio_nextgen_usability_improvements/test_data_bulk/metadata.csv"
# species name

# vep_species stores the original string provided by the user in the config file
vep_species = my_args[4]
# vep_species="Saccharomyces_cerevisiae"
# replace the underscore with a space to get the name of the organism
my_species = sub("_", " ", vep_species)

# get gtf file location
gtf_location = my_args[5]

path_to_scripts = my_args[6]

print(" --- Preparing data...")

## Preparing count data
# counts_data <- read.csv('tximport-counts.csv', header = TRUE, sep = ",")
counts_data <- read.csv(counts_file, header = TRUE, sep = ",")

head(counts_data)

## preparing samples info or metadata
# metaData <- read.table('metadata.csv', header = TRUE, sep = ",")
metaData <- read.table(metadata_file, header = TRUE, sep = ",")

metaData

## Get condition for DESeq object
# condition

## Construct a DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = metaData,
                              design = ~ panel, tidy = TRUE)

dds

## Check size of dds rows x columns
dim(dds)

## Remove genes with less than 10 reads in total for the 6 samples
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

## Check size of dds rows x columns
dim(dds)

## Run DESeq
dds <- DESeq(dds)

## Extract the results
res <- results(dds, lfcThreshold = 0.585, altHypothesis = "greaterAbs", alpha = 0.05)

## remove rows with NA values on the padj column
res <- res[! is.na(res$padj), ]

## Explore the results
summary(res)

## Retain the signifficantly differentially expressed genes
res_sgnf <- res[res$padj <= 0.05, ]
sgnfGenes <- rownames(res_sgnf)
res_sgnf

print(" --- Retain the significantly DE genes that are upregulated in RT vs. AE")

## Retain the significantly DE genes that are upregulated in RT vs. AE
res_sgnf_RT = res_sgnf[res_sgnf$log2FoldChange > 0, ]
sgnfGenes_RT = rownames(res_sgnf_RT)

res_sgnf_AE = res_sgnf[res_sgnf$log2FoldChange < 0, ]
sgnfGenes_AE = rownames(res_sgnf_AE)

## Apply a variance stabilizing transformation on the read count data, for plotting PCA and sample distances
vsd = vst(dds, blind = FALSE)

print(" --- Creating MA plot")

## MA Plot
png("bulk_rna_seq_MA.png", width = 950, height = 650)
DESeq2::plotMA(res, alpha = 0.05, colSig = "blue", cex = 1, main="MA Plot")
dev.off()

print(" --- Computing sample distances")

## Heatmap of sample distances
sampleDists = dist(t(assay(vsd)))
sampleDists
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(vsd$panel, vsd$X, sep = "-")
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

print(" --- Creating heatmap plot for sample distances")

png("bulk_rna_seq_heatmap_distances.png", width = 850, height = 900)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, main="Sample distances")
dev.off()

## sample PCA plot
# first method
# plotPCA(vsd, intgroup = c("panel", "X"))
# second method

pcaData <- plotPCA(vsd, intgroup = c("panel", "X"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

png("bulk_rna_seq_PCA-plot.png", width = 700, height = 500)
print(" --- Creating PA plot")

ggplot(pcaData, aes(PC1, PC2, color=panel, shape=X), title="PCA") +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

## Heatmap for the top 30 significantly differentially expressed genes
## Order the data by row variance in decreasing order and retain the top 30 rows
topDEGenes = head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)

## Select the variance stabilized data for the topDEGenes
mat1 = assay(vsd)[ topDEGenes, ]

## Normalize the data by subtracting row means
mat1 = mat1 - rowMeans(mat1)

## Create a data frame with metadata for the heatmap
metadata_heatmap = as.data.frame(colData(vsd)[, "panel"])
rownames(metadata_heatmap) = colnames(mat1)

print(" --- Creating top significantly differentially expressed pheatmap plot")

png("bulk_rna_seq_top_genes.png", width = 800, height = 900)
pheatmap(mat1, annotation_col = metadata_heatmap, main="Top significantlly expressed genes")
dev.off()

workflow_name = "bulk_rna_seq"
system(paste("Rscript --vanilla ", path_to_scripts, "/downstreamAnalysis/computeMetrics.R", sgnfGenes_RT, sgnfGenes_AE, my_species, workflow_name))
