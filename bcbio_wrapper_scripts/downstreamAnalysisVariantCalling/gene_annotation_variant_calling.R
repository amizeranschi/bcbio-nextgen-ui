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
# query(ah_species, "OrgDb")
qu = query(ah_species, c("OrgDb", my_species))
if (length(qu) > 0) {
  database_species <- qu[[1]]
} else {
    database_species = NULL
    print("No OrgDb found for the specified organism in AnnotationHub! Exiting gene annotation.")
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

# get data
gatk_data = read.table(file, header = F, stringsAsFactors = F)

print(" --- Preparing data...")
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

print(" --- Cluster profile for representation of gene ontology")
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

print(" --- Cluster profile for representation of gene ontology plots")

# plot the results
if (!is.null(ego_high)) {
  jpeg("variant_calling_ego_high.jpg", width = 350, height = 350)
  barplot(ego_high, showCategory=20, title="Enrichment Analysis - representation of gene ontology on HIGH impact genes") 
  dev.off()
} else {
    print("Failed to return results of Gene Ontology for the given list of HIGH impact genes!")
}

if (!is.null(moderate_impact_genes)) {
  jpeg("variant_calling_ego_moderate.jpg", width = 350, height = 350)
  barplot(ego_moderate, showCategory=20, title="Enrichment Analysis - representation of gene ontology on MODERATE impact genes") 
  dev.off()
} else {
    print("Failed to return results of Gene Ontology for the given list of MODERATE impact genes!")
}

if (!is.null(ego_moderate_high)) {
  jpeg("variant_calling_ego_moderate_high.jpg", width = 350, height = 350)
  barplot(ego_moderate_high, showCategory=20, title="Enrichment Analysis - representation of gene ontology on HIGH and MODERATE impact genes") 
  dev.off()
} else {
    print("Failed to return results of Gene Ontology for the given list of HIGH and MODERATE impact genes!")
}

print(" --- KEGG pathway over-representation analysis on the extracted genes")

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

print(" --- KEGG pathway over-representation analysis on the extracted genes plots")

# plot the results
if (!is.null(kegg_high)) {
  jpeg("variant_calling_kegg_high.jpg", width = 350, height = 350)
  barplot(kegg_high, showCategory=20, title="KEGG pathway over-representation analysis on the HIGH impact genes")
  dev.off()  
} else {
    print("Failed to return results of KEGG analysis for the given list of HIGH impact genes!")
  
}

if (!is.null(kegg_moderate)) {
  jpeg("variant_calling_kegg_moderate.jpg", width = 350, height = 350)
  barplot(kegg_moderate, showCategory=20, title="KEGG pathway over-representation analysis on the MODERATE impact genes")
  dev.off()
} else {
    print("Failed to return results of KEGG analysis for the given list of MODERATE impact genes!")
}


if(!is.null(kegg_high_moderate)) {
  jpeg("variant_calling_kegg_high_moderate.jpg", width = 350, height = 350)
  barplot(kegg_high_moderate, showCategory=20, title="KEGG pathway over-representation analysis on the HIGH and MODERATE impact genes")
  dev.off()
} else {
    print("Failed to return results of KEGG analysis for the given list of HIGH and MODERATE impact genes!")
}

print(" --- Over-representation analysis for disease ontology")

# Over-representation analysis for disease ontology
x_high <- enrichDO(gene     = entrezIDs_high$ENTREZID,
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff  = 0.05,
              readable      = FALSE)
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

print(" --- Over-representation analysis for disease ontology plots")

# plot the results
if (!is.null(x_high)) {
  jpeg("variant_calling_x_high.jpg", width = 350, height = 350)
  barplot(x_high, showCategory=20, title="Over-representation analysis for disease ontology in HIGH impact genes")
  dev.off()
} else {
    print("Failed to return results of Disease Ontology for the given list of HIGH impact genes!")
}

if (!is.null(x_moderate)) {
  jpeg("variant_calling_x_moderate.jpg", width = 350, height = 350)
  barplot(x_moderate, showCategory=20, title="Over-representation analysis for disease ontology in MODERATE impact genes")
  dev.off()
} else {
    print("Failed to return results of Disease Ontology for the given list of MODERATE impact genes!")
}

if (!is.null(x_high_moderate)) {
  jpeg("variant_calling_x_moderate_high.jpg", width = 350, height = 350)
  barplot(x_high_moderate, showCategory=20, title="Over-representation analysis for disease ontology in HIGH and MODERATE impact genes")
  dev.off()
} else {
    print("Failed to return results of Disease Ontology for the given list of HIGH and MODERATE impact genes!")
}

# MeSH enrichment analysis
# ah <- AnnotationHub(localHub=TRUE)
qu <- query(ah_species, c("MeSHDbi", my_species))
if (length(qu) > 0) {
  print(" --- MeSH enrichment analysis...")

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
    jpeg("variant_calling_mesh_high.jpg", width = 350, height = 350)
    barplot(mesh_high)
    dev.off()
  } else {
      print("Failed to return MeSH enrichment analysis results for the given list of HIGH impact genes!")
  }

  if (!is.null(de_moderate)) {
    jpeg("variant_calling_mesh_moderate.jpg", width = 350, height = 350)
    barplot(mesh_moderate)
    dev.off()
  } else {
      print("Failed to return MeSH enrichment analysis results for the given list of MODERATE impact genes!")
  }

  if (!is.null(de_high_moderate)) {
    jpeg("variant_calling_mesh_high_moderate.jpg", width = 350, height = 350)
    barplot(mesh_high_moderate)
    dev.off()
  } else {
      print("Failed to return MeSH enrichment analysis results for the given list of HIGH and MODERATE impact genes!")
  }
} else {
    print("No MeSHDb found for the specified organism in AnnotationHub!")
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

