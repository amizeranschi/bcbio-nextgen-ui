## Clear all objects from the memory:
rm(list = ls())
check_install = function(packages) {
  not_installed = setdiff(packages, rownames(installed.packages()))
  if(length(not_installed) > 0) {
    write(paste("The libraries", not_installed, "are not available, so they are being installed now.", sep=" "), stdout())
    
    ## use BiocManager::install() instead of install.packages()
    if (!requireNamespace("BiocManager", quietly = T))
      install.packages("BiocManager", repos='https://cloud.r-project.org')
    # install.packages("Rcpp")
    # , repos='https://cloud.r-project.org')
    
    BiocManager::install(not_installed, ask = F)
  }
}

packages = c("AnnotationHub", "AnnotationDbi", "MeSHDbi", "ggplot2", "clusterProfiler", "DOSE", "meshes", "ChIPseeker", "GOSemSim","ReactomePA",  "VariantAnnotation", "pheatmap", "RColorBrewer",  "DESeq2", "enrichplot")
install.packages('GOplot', repos='https://cloud.r-project.org')
check_install(packages)

installed = lapply(packages, library, character.only = T)
library(GOplot)
## TODO
## For wrapper integration:
##  * give the R script arguments for:
##      * work directory,
##      * DESQ2 condition,
##      * path to final or directly the csv files paths
##  * copy csv/input files to work directory
args = commandArgs(trailingOnly=TRUE)
my_args <- head(args)

## Clear all objects from the memory:
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
# gtf_location="/home/maria/bcbio_nextgen/genomes/Scerevisiae/sacCer3/rnaseq/ref-transcripts.gtf"

# get database for the species using my_species to search AnnotationHub
ah_species <- AnnotationHub()
# query(ah_species, "OrgDb")
qu = query(ah_species, c("OrgDb", my_species))
if (length(qu) > 0) {
  database_species <- qu[[1]]
} else {
    database_species = NULL
    print(" --- No OrgDb found for the specified organism in AnnotationHub! Exiting gene annotation.")
    quit()
}

species_keytype = "ORF"
# create a custom Transcript database based on the GTF file from Bcbio's genome
# create a file name for the TxDB file, which is an SQLite DB, by appending the
# extension ".sqlite" to gtf_location
txdb_file = paste0(gtf_location, ".sqlite")

# use gtf_location to create an SQLite file called txdb_file, if the latter not
# exist already
if(!file.exists(txdb_file))
{
  txdb_gtf = GenomicFeatures::makeTxDbFromGFF(gtf_location, organism = my_species)
  AnnotationDbi::saveDb(txdb_gtf, txdb_file)
}

# create a TxDB object from txdb_file that will later be used as a instance
my_txdb = loadDb(txdb_file)

## Set working directory
# setwd("C:/Users/40724/Documents/downstreamAnalysis")

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
png("MA.png", width = 950, height = 650)
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

png("heatmap_distances.png", width = 850, height = 900)

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

png("PCA-plot.png", width = 700, height = 500)
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

png("top_genes.png", width = 800, height = 900)
pheatmap(mat1, annotation_col = metadata_heatmap, main="Top significantlly expressed genes")
dev.off()

## gene annotation
# use the database_species to get a list of Entrez IDs corresponding to the
# significant gene IDs
entrezIDs_RT = AnnotationDbi::select(database_species, keys = sgnfGenes_RT, keytype=species_keytype, columns = "ENTREZID")
entrezIDs_AE = AnnotationDbi::select(database_species, keys = sgnfGenes_AE, keytype=species_keytype, columns = "ENTREZID")

# use cluster profile for representation of gene ontology on each list of
# interest genes
print(" --- Cluster profile for representation of gene ontology")

ego_RT <- enrichGO(gene       = entrezIDs_RT$ENTREZID,
                OrgDb         = database_species,
                keyType       = 'ENTREZID',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

head(summary(ego_RT), n=3)
head(ego_RT)
ego_AE = enrichGO(gene = entrezIDs_AE$ENTREZID,
               keyType       = "ENTREZID",
               OrgDb         = database_species,
               ont           = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)
head(ego_AE)

print(" --- Cluster profile for representation of gene ontology plots")

# plot the results
if (!is.null(ego_RT)) {
  png("ego_RT.png", width = 750, height = 750)
  barplot(ego_RT, showCategory=20, title="Enrichment Analysis - representation of gene ontology on ROOT samples") 
  dev.off()
} else {
    print("Failed to return results of Gene Ontology for ROOT!")
}
if (!is.null(ego_AE)) {
  png("ego_AE.png", width = 750, height = 950)
  barplot(ego_AE, showCategory=20, title="Enrichment Analysis - representation of gene ontology on AERIAL samples") 
  dev.off()
} else {
    print("Failed to return results of Gene Ontology for AERIAL!")
}

# search the organism for kegg analysis
kegg_organism = search_kegg_organism(my_species)$kegg_code

print(" --- KEGG pathway over-representation analysis on the extracted genes")

# KEGG pathway over-representation analysis on the extracted genes
kegg_RT <- enrichKEGG(gene           = sgnfGenes_RT,
                        organism     = kegg_organism,
                        pvalueCutoff = 0.05)
head(kegg_RT)

kegg_AE <- enrichKEGG(gene               = sgnfGenes_AE,
                            organism     = kegg_organism,
                            pvalueCutoff = 0.05)
head(kegg_AE)

print(" --- KEGG pathway over-representation analysis on the extracted genes plots")

# plot the results
if (!is.null(kegg_RT)) {
  png("kegg_RT.png", width = 600, height = 650)
  barplot(kegg_RT, showCategory=20, title="KEGG pathway over-representation analysis on the ROOT genes")
  dev.off()  
} else {
  print("Failed to return results of KEGG analysis for the given list of ROOT significant genes!")
  
}

if (!is.null(kegg_AE)) {
  png("kegg_AE.png", width = 600, height = 650)
  barplot(kegg_AE, showCategory=20, title="KEGG pathway over-representation analysis on the AERIAL genes")
  dev.off()
} else {
  print("Failed to return results of KEGG analysis for the given list of AERIAL significant genes!")
}

print(" --- Over-representation analysis for disease ontology")

# Over-representation analysis for disease ontology
entrez_RT <- c(entrezIDs_RT$ENTREZID)
head(entrez_RT)
# remove <- NA
# if (any (is.na (remove))) 
#   entrez_RT <- entrez_RT [! is.na (entrez_RT)]
# head(entrez_RT)
# typeof(entrez_RT)

x_RT <- enrichDO(gene          = entrez_RT,
                 ont           = "DO",
                 pvalueCutoff  = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff  = 0.05,
                 readable      = FALSE)
head(x_RT)

x_AE <- enrichDO(gene                = entrezIDs_AE$ENTREZID,
                 ont           = "DO",
                 pvalueCutoff  = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff  = 0.05,
                 readable      = FALSE)
head(x_AE)

print(" --- Over-representation analysis for disease ontology plots")

# plot the results
if (!is.null(x_RT)) {
  png("x_RT.png", width = 600, height = 650)
  barplot(x_RT, showCategory=20, title="Over-representation analysis for disease ontology in ROOT genes")
  dev.off()
} else {
  print("Failed to return results of Disease Ontology for the given list of ROOT significant genes!")
}

if (!is.null(x_AE)) {
  png("x_AE.png", width = 600, height = 650)
  barplot(x_AE, showCategory=20, title="Over-representation analysis for disease ontology in AERIAL genes")
  dev.off()
} else {
  print("Failed to return results of Disease Ontology for the given list of AERIAL significant genes!")
}



# MeSH enrichment analysis
qu <- query(ah_species, c("MeSHDbi", my_species))
if (length(qu) > 0) {
  print(" --- MeSH enrichment analysis...")

  file_qu <- qu[[1]]
  db <- MeSHDbi::MeSHDb(file_qu)

  de_RT <- names(sgnfGenes_RT)[1:100]
  mesh_RT <- enrichMeSH(de_RT, MeSHDb = db, database='gendoo', category = 'C')

  de_AE <- names(sgnfGenes_AE)[1:100]
  mesh_AE <- enrichMeSH(de_AE, MeSHDb = db, database='gendoo', category = 'C')

  # plot the results
  if (!is.null(de_RT)) {
    png("mesh_RT.png", width = 600, height = 650)
    barplot(mesh_RT)
    dev.off()
  } else {
      print("Failed to return MeSH enrichment analysis results for the given list of ROOT significant genes!")
  }

  if (!is.null(de_AE)) {
    png("mesh_AE.png", width = 600, height = 650)
    barplot(mesh_AE)
    dev.off()
  } else {
      print("Failed to return MeSH enrichment analysis results for the given list of AERIAL significant genes!")
  }
}