
# load the required packages
packages = c("AnnotationHub", "AnnotationDbi", "MeSHDbi", "ggplot2", "clusterProfiler", "DOSE", "meshes", "ChIPseeker", "GOSemSim", "GOplot", "ReactomePA",  "VariantAnnotation", "pheatmap", "RColorBrewer",  "DESeq2", "enrichplot")
installed = lapply(packages, library, character.only = T)

args = commandArgs(trailingOnly=TRUE)
my_args <- head(args)

sgnfGenes_set1 = args[1]
sgnfGenes_set2 = args[2]
my_species = args[3]
my_species = paste(my_species, args[4])
workflow_name = args[5]
# gtf_location = args[6]

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
# txdb_file = paste0(gtf_location, ".sqlite")

# use gtf_location to create an SQLite file called txdb_file, if the latter not
# exist already
# if(!file.exists(txdb_file))
# {
#   txdb_gtf = GenomicFeatures::makeTxDbFromGFF(gtf_location, format = "gtf", organism = my_species)
#   AnnotationDbi::saveDb(txdb_gtf, txdb_file)
# }

# # create a TxDB object from txdb_file that will later be used as a instance
# my_txdb = loadDb(txdb_file)

## gene annotation
# use the database_species to get a list of Entrez IDs corresponding to the
# significant gene IDs
entrezIDs_set1 = AnnotationDbi::select(database_species, keys = sgnfGenes_set1, keytype=species_keytype, columns = "ENTREZID")
entrezIDs_set2 = AnnotationDbi::select(database_species, keys = sgnfGenes_set2, keytype=species_keytype, columns = "ENTREZID")

# use cluster profile for representation of gene ontology on each list of
# interest genes
print(" --- Cluster profile for representation of gene ontology")

ego_set1 <- enrichGO(gene       = entrezIDs_set1$ENTREZID,
                OrgDb         = database_species,
                keyType       = 'ENTREZID',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
name_set1 = paste(workflow_name,"ego_set1",sep="_")
file_name_set1 = paste(name_set1, ".txt")
figure_name_set1 = paste(name_set1, ".png")
write.table(ego_set1, file = file_name_set1, sep = "\t", quote = F, row.names = F, na = "")

head(summary(ego_set1), n=3)
head(ego_set1)
ego_set2 = enrichGO(gene = entrezIDs_set2$ENTREZID,
               keyType       = "ENTREZID",
               OrgDb         = database_species,
               ont           = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)
name_set2 = paste(workflow_name,"ego_set2",sep="_")
file_name_set2 = paste(name_set2, ".txt")
figure_name_set2 = paste(name_set2, ".png")
write.table(ego_set2, file = file_name_set2, sep = "\t", quote = F, row.names = F, na = "")

head(ego_set2)

print(" --- Cluster profile for representation of gene ontology plots")

# plot the results
if (!is.null(ego_set1)) {
  png(figure_name_set1, width = 750, height = 750)
  barplot(ego_set1, showCategory=20, title="Enrichment Analysis - representation of gene ontology on set1 significant samples") 
  dev.off()
} else {
    print("Failed to return results of Gene Ontology for set1 significant!")
}
if (!is.null(ego_set2)) {
  png(figure_name_set2, width = 750, height = 950)
  barplot(ego_set2, showCategory=20, title="Enrichment Analysis - representation of gene ontology on set2 significant samples") 
  dev.off()
} else {
    print("Failed to return results of Gene Ontology for set2 significant!")
}

# search the organism for kegg analysis
kegg_organism = search_kegg_organism(my_species)$kegg_code

print(" --- KEGG pathway over-representation analysis on the extracted genes")

# KEGG pathway over-representation analysis on the extracted genes
kegg_set1 <- enrichKEGG(gene           = sgnfGenes_set1,
                        organism     = kegg_organism,
                        pvalueCutoff = 0.05)

kegg_name_set1 = paste(workflow_name,"kegg_set1",sep="_")
file_name_kegg_set1 = paste(kegg_name_set1, ".txt")
figure_name_kegg_set1 = paste(kegg_name_set1, ".png")
write.table(kegg_set1, file = file_name_kegg_set1, sep = "\t", quote = F, row.names = F, na = "")

head(kegg_set1)

kegg_set2 <- enrichKEGG(gene               = sgnfGenes_set2,
                            organism     = kegg_organism,
                            pvalueCutoff = 0.05)

kegg_name_set2 = paste(workflow_name,"kegg_set2",sep="_")
file_name_kegg_set2 = paste(kegg_name_set2, ".txt")
figure_name_kegg_set2 = paste(kegg_name_set2, ".png")
write.table(kegg_set2, file = file_name_kegg_set2, sep = "\t", quote = F, row.names = F, na = "")

head(kegg_set2)

print(" --- KEGG pathway over-representation analysis on the extracted genes plots")

# plot the results
if (!is.null(kegg_set1)) {
  png(figure_name_kegg_set1, width = 600, height = 650)
  barplot(kegg_set1, showCategory=20, title="KEGG pathway over-representation analysis on the set1 significant genes")
  dev.off()  
} else {
  print("Failed to return results of KEGG analysis for the given list of set1 significant significant genes!")
  
}

if (!is.null(kegg_set2)) {
  png(figure_name_kegg_set2, width = 600, height = 650)
  barplot(kegg_set2, showCategory=20, title="KEGG pathway over-representation analysis on the set2 significant genes")
  dev.off()
} else {
  print("Failed to return results of KEGG analysis for the given list of set2 significant significant genes!")
}

print(" --- Over-representation analysis for disease ontology")

# Over-representation analysis for disease ontology
entrez_set1 <- c(entrezIDs_set1$ENTREZID)
head(entrez_set1)
# remove <- NA
# if (any (is.na (remove))) 
#   entrez_set1 <- entrez_set1 [! is.na (entrez_set1)]
# head(entrez_set1)
# typeof(entrez_set1)

x_set1 <- enrichDO(gene          = entrez_set1,
                 ont           = "DO",
                 pvalueCutoff  = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff  = 0.05,
                 readable      = FALSE)

x_name_set1 = paste(workflow_name,"disease_ont_set1",sep="_")
file_name_x_set1 = paste(x_name_set1, ".txt")
figure_name_x_set1 = paste(x_name_set1, ".png")
write.table(x_set1, file = file_name_x_set1, sep = "\t", quote = F, row.names = F, na = "")

head(x_set1)

x_set2 <- enrichDO(gene                = entrezIDs_set2$ENTREZID,
                 ont           = "DO",
                 pvalueCutoff  = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff  = 0.05,
                 readable      = FALSE)
x_name_set2 = paste(workflow_name,"disease_ont_set2",sep="_")
file_name_x_set2 = paste(x_name_set2, ".txt")
figure_name_x_set2 = paste(x_name_set2, ".png")
write.table(x_set2, file = file_name_x_set2, sep = "\t", quote = F, row.names = F, na = "")

head(x_set2)

print(" --- Over-representation analysis for disease ontology plots")

# plot the results
if (!is.null(x_set1)) {
  png(figure_name_x_set1, width = 600, height = 650)
  barplot(x_set1, showCategory=20, title="Over-representation analysis for disease ontology in set1 significant genes")
  dev.off()
} else {
  print("Failed to return results of Disease Ontology for the given list of set1 significant significant genes!")
}

if (!is.null(x_set2)) {
  png(figure_name_x_set2, width = 600, height = 650)
  barplot(x_set2, showCategory=20, title="Over-representation analysis for disease ontology in set2 significant genes")
  dev.off()
} else {
  print("Failed to return results of Disease Ontology for the given list of set2 significant significant genes!")
}



# MeSH enrichment analysis
qu <- query(ah_species, c("MeSHDbi", my_species))
if (length(qu) > 0) {
  print(" --- MeSH enrichment analysis...")

  file_qu <- qu[[1]]
  db <- MeSHDbi::MeSHDb(file_qu)

  de_set1 <- names(sgnfGenes_set1)[1:100]
  mesh_set1 <- enrichMeSH(de_set1, MeSHDb = db, database='gendoo', category = 'C')
 
  mesh_name_set1 = paste(workflow_name,"mesh_set1",sep="_")
  file_name_mesh_set1 = paste(mesh_name_set1, ".txt")
  figure_name_mesh_set1 = paste(mesh_name_set1, ".png")
  write.table(mesh_set1, file = file_name_mesh_set1, sep = "\t", quote = F, row.names = F, na = "")

  de_set2 <- names(sgnfGenes_set2)[1:100]
  mesh_set2 <- enrichMeSH(de_set2, MeSHDb = db, database='gendoo', category = 'C')

  mesh_name_set2 = paste(workflow_name,"mesh_set2",sep="_")
  file_name_mesh_set2 = paste(mesh_name_set2, ".txt")
  figure_name_mesh_set2 = paste(mesh_name_set2, ".png")
  write.table(mesh_set2, file = file_name_mesh_set2, sep = "\t", quote = F, row.names = F, na = "")

  # plot the results
  if (!is.null(de_set1)) {
    png(figure_name_mesh_set1, width = 600, height = 650)
    barplot(mesh_set1)
    dev.off()
  } else {
      print("Failed to return MeSH enrichment analysis results for the given list of set1 significant significant genes!")
  }

  if (!is.null(de_set2)) {
    png(figure_name_mesh_set2, width = 600, height = 650)
    barplot(mesh_set2)
    dev.off()
  } else {
      print("Failed to return MeSH enrichment analysis results for the given list of set2 significant significant genes!")
  }
}
