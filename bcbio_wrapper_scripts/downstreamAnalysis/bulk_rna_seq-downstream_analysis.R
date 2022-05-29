## Clear all objects from the memory:
rm(list = ls())

## The following function tests whether a list of packages is installed and installs them if they are not available:
check_install = function(packages) {
  not_installed = setdiff(packages, rownames(installed.packages()))
  if(length(not_installed) > 0) {
    write(paste("The libraries", not_installed, "are not available, so they are being installed now.", sep=" "), stdout())
    
    ## use BiocManager::install() instead of install.packages()
    if (!requireNamespace("BiocManager", quietly = T))
      install.packages("BiocManager")
    
    BiocManager::install(not_installed, ask = F)
  }
}

## Install and load the required packages:
packages = c("pheatmap", "RColorBrewer", "ChIPseeker", "ReactomePA", "org.Sc.sgd.db", "DESeq2", "bedr", "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", "BSgenome.Scerevisiae.UCSC.sacCer3", "VariantAnnotation", "SBGNview", "pathview", "gage", "gageData", "clusterProfiler", "topGO", "biomaRt", "AnnotationHub", "AnnotationDbi")
check_install(packages)
installed = lapply(packages, library, character.only = T)

## TODO
## For wrapper integration:
##  * give the R script arguments for:
##      * work directory,
##      * DESQ2 condition,
##      * path to final or directly the csv files paths
##  * copy csv/input files to work directory


## Set working directory
setwd("C:/Users/40724/Documents/downstreamAnalysis")

## Preparing count data
counts_data <- read.csv('tximport-counts.csv', header = TRUE, sep = ",")
head(counts_data)

## preparing samples info or metadata
metaData <- read.table('metadata.csv', header = TRUE, sep = ",")
metaData

## Get condition for DESeq object
# condition = colnames(metaData)[2]
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


## Retain the significantly DE genes that are upregulated in RT vs. AE
res_sgnf_RT = res_signif[res_signif$log2FoldChange > 0, ]
sgnfGenes_RT = rownames(res_sgnf_RT)

res_sgnf_AE = res_sgnif[res_sgnf$log2FoldChange < 0, ]
sgnfGenes_AE = rownames(res_sgnf_AE)

## Apply a variance stabilizing transformation on the read count data, for plotting PCA and sample distances
vsd = vst(dds, blind = FALSE)

## MA Plot
DESeq2::plotMA(res, alpha = 0.05, colSig = "blue", cex = 1)

## Heatmap of sample distances
sampleDists = dist(t(assay(vsd)))
sampleDists
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(vsd$panel, vsd$X, sep = "-")
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

## sample PCA plot
# first method
# plotPCA(vsd, intgroup = c("panel", "X"))
# second method
pcaData <- plotPCA(vsd, intgroup = c("panel", "X"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=panel, shape=X)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


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

pheatmap(mat1, annotation_col = metadata_heatmap)

