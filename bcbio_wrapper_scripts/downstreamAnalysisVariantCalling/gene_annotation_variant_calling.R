rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
my_args <- head(args)

## Clear all objects from the memory:
# Set working directory
work_dir=my_args[1]
setwd(work_dir)

# set file with genomic data
file=my_args[2]
# species name
# vep_species stores the original string provided by the user in the config file
vep_species = my_args[3]
# replace the underscore with a space to get the name of the organism
my_species = sub("_", " ", vep_species)

# organism_type = my_args[4]
# get gtf file location
gtf_location = my_args[4]
# check if all packages are installed, if not install them
source("set_packages.R")

# get database for the species using my_species to search AnnotationHub
ah_species <- AnnotationHub()
query(ah_species, "OrgDb")
database_species <- query(ah_species, c("OrgDb", my_species))[[1]]

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

# get data
gatk_data = read.table(file, header = F, stringsAsFactors = F)

# set the column names
# TODO check if mammal genome or  yeast
colnames(gatk_data) = c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "IMPACT", "DISTANCE", "STRAND", "FLAGS", "BIOTYPE", "CLIN_SIG", "SOMATIC", "PHENO")

# select variants with "HIGH" impact, "MODERATE" impact, or either of them
high_impact_vars = gatk_data[(gatk_data$IMPACT == "HIGH"), ]
moderate_impact_vars = gatk_data[(gatk_data$IMPACT == "MODERATE"), ]
high_moderate_impact_vars = gatk_data[(gatk_data$IMPACT == "HIGH" | gatk_data$IMPACT == "MODERATE"), ]

# retain the genes corresponding to high- and moderate-impact variants
# genes high-impact variants
high_impact_genes = unique(high_impact_vars[, c("Gene")])

# genes with moderate-impact variants
moderate_impact_genes = unique(moderate_impact_vars[, c("Gene")])

# genes with high-impact OR moderate-impact variants
high_moderate_impact_genes = unique(high_moderate_impact_vars[, c("Gene")])

# TODO figure out how to choose the key type
# use the database_species to get a list of Entrez IDs corresponding to the
# significant gene IDs
entrezIDs_high = AnnotationDbi::select(database_species, keys = high_impact_genes, keytype=species_keytype, columns = "ENTREZID")
entrezIDs_moderate = AnnotationDbi::select(database_species, keys = moderate_impact_genes, keytype=species_keytype, columns = "ENTREZID")
entrezIDs_moderate_high = AnnotationDbi::select(database_species, keys = high_moderate_impact_genes, keytype=species_keytype, columns = "ENTREZID")


# use cluster profile for representation of gene ontology on each list of
# interest genes
ego_high = enrichGO(gene       = entrezIDs_high$ENTREZID,
                 keyType       = "ENTREZID",
                 OrgDb         = database_species,
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = FALSE)
head(ego_high)

ego_moderate = enrichGO(gene = entrezIDs_moderate$ENTREZID,
               keyType       = "ENTREZID",
               OrgDb         = database_species,
               ont           = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05,
               readable      = FALSE)
head(ego_moderate)

ego_moderate_high = enrichGO(gene     = entrezIDs_moderate_high$ENTREZID,
                        keyType       = "ENTREZID",
                        OrgDb         = database_species,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = FALSE)
head(ego_moderate_high)

# plot the results
if (!is.null(ego_high)) {
  jpeg("ego_high.jpg", width = 350, height = 350)
  ggplot(ego_high)
  dev.off()
} else {
    print("Failed to return results of Gene Ontology for the given list of HIGH impact genes!")
}

if (!is.null(moderate_impact_genes)) {
  jpeg("ego_moderate.jpg", width = 350, height = 350)
  ggplot(ego_moderate)
  dev.off()
} else {
    print("Failed to return results of Gene Ontology for the given list of MODERATE impact genes!")
}

if (!is.null(ego_moderate_high)) {
  jpeg("ego_moderate_high.jpg", width = 350, height = 350)
  ggplot(ego_moderate_high)
  dev.off()
} else {
    print("Failed to return results of Gene Ontology for the given list of HIGH and MODERATE impact genes!")
}

# search the organism for kegg analysis
kegg_organism = search_kegg_organism(my_species)$kegg_code

# KEGG pathway over-representation analysis on the extracted genes with the
kegg_high <- enrichKEGG(gene        = high_impact_genes,
                       organism     = kegg_organism,
                       pvalueCutoff = 0.05)
head(kegg_high)

kegg_moderate <- enrichKEGG(gene         = moderate_impact_genes,
                            organism     = kegg_organism,
                            pvalueCutoff = 0.05)
head(kegg_moderate)

kegg_high_moderate <- enrichKEGG(gene         = high_moderate_impact_genes,
                                 organism     = kegg_organism,
                                 pvalueCutoff = 0.05)
head(kegg_high_moderate)

# plot the results
if (!is.null(kegg_high)) {
  jpeg("kegg_high.jpg", width = 350, height = 350)
  ggplot(kegg_high)
  dev.off()  
} else {
    print("Failed to return results of KEGG analysis for the given list of HIGH impact genes!")
  
}

if (!is.null(kegg_moderate)) {
  jpeg("kegg_moderate.jpg", width = 350, height = 350)
  ggplot(kegg_moderate)
  dev.off()
} else {
    print("Failed to return results of KEGG analysis for the given list of MODERATE impact genes!")
}


if(!is.null(kegg_high_moderate)) {
  jpeg("kegg_high_moderate.jpg", width = 350, height = 350)
  ggplot(kegg_moderate_high)
  dev.off()
} else {
    print("Failed to return results of KEGG analysis for the given list of HIGH and MODERATE impact genes!")
}

# Over-representation analysis for disease ontology
x_high <- enrichDO(gene     = entrezIDs_high$ENTREZID,
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff  = 0.05,
              readable      = FALSE)database_species
head(x_high)

x_moderate <- enrichDO(gene          = entrezIDs_moderate$ENTREZID,
                   pvalueCutoff  = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = FALSE)
head(x_moderate)

x_high_moderate <- enrichDO(gene          = entrezIDs_moderate_high$ENTREZID,
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05,
                       readable      = FALSE)
head(x_high_moderate)

# plot the results
if (!is.null(x_high)) {
  jpeg("x_high.jpg", width = 350, height = 350)
  ggplot(x_high)
  dev.off()
} else {
    print("Failed to return results of Disease Ontology for the given list of HIGH impact genes!")
}

if (!is.null(x_moderate)) {
  jpeg("x_moderate.jpg", width = 350, height = 350)
  ggplot(x_moderate)
  dev.off()
} else {
    print("Failed to return results of Disease Ontology for the given list of MODERATE impact genes!")
}

if (!is.null(x_high_moderate)) {
  jpeg("x_moderate_high.jpg", width = 350, height = 350)
  ggplot(x_high_moderate)
  dev.off()
} else {
    print("Failed to return results of Disease Ontology for the given list of HIGH and MODERATE impact genes!")
}

# MeSH enrichment analysis
# ah <- AnnotationHub(localHub=TRUE)
qu <- query(ah_species, c("MeSHDbi", my_species))
if (length(qu) > 0) {
  file_qu <- qu[[1]]
  db <- MeSHDbi::MeSHDb(file_qu)

  de_high <- names(high_impact_genes)[1:100]
  mesh_high <- enrichMeSH(de_high, MeSHDb = db, database='gendoo', category = 'C')

  de_moderate <- names(moderate_impact_genes)[1:100]
  mesh_moderate <- enrichMeSH(de_moderate, MeSHDb = db, database='gendoo', category = 'C')

  de_high_moderate <- names(high_moderate_impact_genes)[1:100]
  mesh_high_moderate <- enrichMeSH(de_high_moderate, MeSHDb = db, database='gendoo', category = 'C')

  # plot the results
  if (!is.null(de_high)) {
    jpeg("mesh_high.jpg", width = 350, height = 350)
    ggplot(mesh_high)
    dev.off()
  } else {
      print("Failed to return MeSH enrichment analysis results for the given list of HIGH impact genes!")
  }

  if (!is.null(de_moderate)) {
    jpeg("mesh_moderate.jpg", width = 350, height = 350)
    ggplot(mesh_moderate)
    dev.off()
  } else {
      print("Failed to return MeSH enrichment analysis results for the given list of MODERATE impact genes!")
  }

  if (!is.null(de_high_moderate)) {
    jpeg("mesh_high_moderate.jpg", width = 350, height = 350)
    ggplot(mesh_high_moderate)
    dev.off()
  } else {
      print("Failed to return MeSH enrichment analysis results for the given list of HIGH and MODERATE impact genes!")
  }
} else {
    print("No MeSHDb for the specified genome!")
}

#############################################################
# possibly giving up on this
#############################################################
# Reactome pathway over-representation analysis
# head(entrezIDs_high)
# reactome_high <- enrichPathway(gene=entrezIDs_high$ENTREZID, organism="yeast", pvalueCutoff = 0.05, readable=TRUE)
# head(reactome_high)

# reactome_moderate <- enrichPathway(gene=entrezIDs_moderate$ENTREZID, organism="yeast", pvalueCutoff = 0.05, readable=TRUE)
# head(reactome_moderate)

# reactome_high_moderate <- enrichPathway(gene=entrezIDs_high_moderate$ENTREZID, organism="yeast", pvalueCutoff = 0.05, readable=TRUE)
# head(reactome_high_moderate)

