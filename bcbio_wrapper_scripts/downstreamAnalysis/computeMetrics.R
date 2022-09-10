
# load the required packages
packages = c("AnnotationHub", "AnnotationDbi", "MeSHDbi", "ggplot2", "clusterProfiler", "DOSE", "meshes", "ChIPseeker", "GOSemSim", "GOplot", "ReactomePA",  "VariantAnnotation", "pheatmap", "RColorBrewer",  "DESeq2", "enrichplot")
installed = lapply(packages, library, character.only = T)

args = commandArgs(trailingOnly=TRUE)
my_args <- head(args)

sgnfGenes_HI = args[1]
sgnfGenes_ME = args[2]
my_species = args[3]
workflow_name = args[4]

# get database for the species using my_species to search AnnotationHub
ah_species = AnnotationHub()

## search for a suitable OrgDb
qu = query(ah_species, c("OrgDb", sub("_", " ", my_species)))
if (length(qu) > 0)
{
  orgdb_species = qu[[1]]
} else 
{
  orgdb_species = NULL
  print("No OrgDb found for the specified organism in AnnotationHub! Exiting gene annotation.")
  quit()
}

## search for a suitable MeSHDb
qu = query(ah_species, c("MeSHDb", sub("_", " ", my_species)))
if (length(qu) > 0)
{
  meshdb_species = MeSHDbi::MeSHDb(qu[[1]])
} else 
{
  meshdb_species = NULL
  print("No OrgDb found for the specified organism in AnnotationHub! Exiting gene annotation.")
  quit()
}



## view the available key types in orgdb_species
# keytypes(orgdb_species)

## set the appropriate key type according to species
if(sub("_", " ", my_species) == "Saccharomyces cerevisiae")
{
  species_keytype = "ORF"
} else 
{
  species_keytype = "SYMBOL"
}

print(paste0(" --- The current species is ", my_species, " and the annotation key type used is ", species_keytype))


## gene annotation
# use orgdb_species to get a list of Entrez IDs corresponding to the significant gene IDs
entrezIDs_HI = AnnotationDbi::select(orgdb_species, keys = sgnfGenes_HI, keytype=species_keytype, columns = "ENTREZID")
entrezIDs_ME = AnnotationDbi::select(orgdb_species, keys = sgnfGenes_ME, keytype=species_keytype, columns = "ENTREZID")



## function for replacing ORFs (for yeast) or SYMBOLs (for any other organism) with GENENAMEs in KEGG and GO enrichment results
## if an NA is returned for a gene, use the ORF instead

makeReadable = function(geneIDs, orig_keytype)
{
  ## check if length(geneIDs) is greater than 0
  if(length(geneIDs) == 0)
  {
    return(character(0))
  }
  ## convert ORFs to GENENAMEs, for yeast
  # genenames = lapply(geneIDs, function(x) { AnnotationDbi::mapIds(org.Sc.sgd.db, keys = unlist(strsplit(x, split = "/")), keytype="ORF", column = "GENENAME", multiVals = "first") })
  
  ## convert ORFs or ENTREZIDs to COMMON gene names
  genenames = lapply(geneIDs, function(x) { AnnotationDbi::mapIds(orgdb_species, keys = unlist(strsplit(x, split = "/")), keytype=orig_keytype, column = "COMMON", multiVals = "first") })
  
  ## replace NA values with original species_keytype (SYMBOL or ORF)
  for(someiter in 1:length(genenames))
  {
    genenames[[someiter]][is.na(genenames[[someiter]])] = names(genenames[[someiter]][is.na(genenames[[someiter]])])
  }
  
  ## paste together the genenames and return the concatenated strings
  genenames = sapply(genenames, paste, collapse = "/")
}




# use clusterProfiler for representation of gene ontology on each list of interest genes
print(" --- Gene Ontology overrepresentation analysis")

GO_HI = enrichGO(gene     = entrezIDs_HI$ENTREZID,
                 OrgDb         = orgdb_species,
                 keyType       = 'ENTREZID',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
name_HI = paste(workflow_name,"GO_HI",sep="_")
file_name_HI = paste(name_HI, ".txt")
figure_name_HI = paste(name_HI, ".png")
write.table(GO_HI, file = file_name_HI, sep = "\t", quote = F, row.names = F, na = "")

head(summary(GO_HI), n=3)
head(GO_HI)
GO_ME = enrichGO(gene = entrezIDs_ME$ENTREZID,
                 keyType       = "ENTREZID",
                 OrgDb         = orgdb_species,
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
name_ME = paste(workflow_name,"GO_ME",sep="_")
file_name_ME = paste(name_ME, ".txt")
figure_name_ME = paste(name_ME, ".png")
write.table(GO_ME, file = file_name_ME, sep = "\t", quote = F, row.names = F, na = "")

head(GO_ME)

print(" --- Cluster profile for representation of gene ontology plots")

# plot the results
if (!is.null(GO_HI)) {
  png(figure_name_HI, width = 750, height = 750)
  barplot(GO_HI, showCategory=20, title="Enrichment Analysis - representation of gene ontology on HI significant samples") 
  dev.off()
} else {
  print("Failed to return results of Gene Ontology for HI significant!")
}
if (!is.null(GO_ME)) {
  png(figure_name_ME, width = 750, height = 950)
  barplot(GO_ME, showCategory=20, title="Enrichment Analysis - representation of gene ontology on ME significant samples") 
  dev.off()
} else {
  print("Failed to return results of Gene Ontology for ME significant!")
}



# search the organism for kegg analysis
kegg_organism = search_kegg_organism(sub("_", " ", my_species))$kegg_code

print(" --- KEGG pathway over-representation analysis on the extracted genes")

# KEGG pathway over-representation analysis on the extracted genes
kegg_HI = enrichKEGG(gene           = sgnfGenes_HI,
                     organism     = kegg_organism,
                     pvalueCutoff = 0.05)

kegg_name_HI = paste(workflow_name,"kegg_HI",sep="_")
file_name_kegg_HI = paste(kegg_name_HI, ".txt")
figure_name_kegg_HI = paste(kegg_name_HI, ".png")
write.table(kegg_HI, file = file_name_kegg_HI, sep = "\t", quote = F, row.names = F, na = "")

head(kegg_HI)

kegg_ME = enrichKEGG(gene                = sgnfGenes_ME,
                     organism     = kegg_organism,
                     pvalueCutoff = 0.05)
kegg_ME_readable = kegg_ME
kegg_ME_readable@result$geneID = makeReadable(kegg_ME@result$geneID, orig_keytype = "ORF")
kegg_name_ME = paste(workflow_name,"kegg_ME",sep="_")
file_name_kegg_ME = paste0(kegg_name_ME, ".txt")
figure_name_kegg_ME = paste0(kegg_name_ME, ".png")
write.table(kegg_ME, file = file_name_kegg_ME, sep = "\t", quote = F, row.names = F, na = "")

head(kegg_ME)

print(" --- KEGG pathway over-representation analysis on the extracted genes plots")

# plot the results
if (!is.null(kegg_HI)) {
  png(figure_name_kegg_HI, width = 600, height = 650)
  barplot(kegg_HI, showCategory=20, title="KEGG pathway over-representation analysis on the HI significant genes")
  dev.off()  
} else {
  print("Failed to return results of KEGG analysis for the given list of HI significant significant genes!")
  
}

if (!is.null(kegg_ME)) {
  png(figure_name_kegg_ME, width = 600, height = 650)
  barplot(kegg_ME, showCategory=20, title="KEGG pathway over-representation analysis on the ME significant genes")
  dev.off()
} else {
  print("Failed to return results of KEGG analysis for the given list of ME significant significant genes!")
}








print(" --- Over-representation analysis for disease ontology")

# Over-representation analysis for disease ontology
entrez_HI = c(entrezIDs_HI$ENTREZID)
head(entrez_HI)
# remove = NA
# if (any (is.na (remove))) 
#   entrez_HI = entrez_HI [! is.na (entrez_HI)]
# head(entrez_HI)
# typeof(entrez_HI)

x_HI = enrichDO(gene          = entrez_HI,
                ont           = "DO",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = FALSE)

x_name_HI = paste(workflow_name,"disease_ont_HI",sep="_")
file_name_x_HI = paste(x_name_HI, ".txt")
figure_name_x_HI = paste(x_name_HI, ".png")
write.table(x_HI, file = file_name_x_HI, sep = "\t", quote = F, row.names = F, na = "")

head(x_HI)

x_ME = enrichDO(gene                = entrezIDs_ME$ENTREZID,
                ont           = "DO",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = FALSE)
x_name_ME = paste(workflow_name,"disease_ont_ME",sep="_")
file_name_x_ME = paste0(x_name_ME, ".txt")
figure_name_x_ME = paste0(x_name_ME, ".png")
write.table(x_ME, file = file_name_x_ME, sep = "\t", quote = F, row.names = F, na = "")

head(x_ME)

print(" --- Over-representation analysis for disease ontology plots")

# plot the results
if (!is.null(x_HI)) {
  png(figure_name_x_HI, width = 600, height = 650)
  barplot(x_HI, showCategory=20, title="Over-representation analysis for disease ontology in HI significant genes")
  dev.off()
} else {
  print("Failed to return results of Disease Ontology for the given list of HI significant significant genes!")
}

if (!is.null(x_ME)) {
  png(figure_name_x_ME, width = 600, height = 650)
  barplot(x_ME, showCategory=20, title="Over-representation analysis for disease ontology in ME significant genes")
  dev.off()
} else {
  print("Failed to return results of Disease Ontology for the given list of ME significant significant genes!")
}



# MeSH enrichment analysis
qu = query(ah_species, c("MeSHDb", sub("_", " ", sub("_", " ", my_species))))
if (length(qu) > 0) {
  print(" --- MeSH enrichment analysis...")
  
  file_qu = qu[[1]]
  db = MeSHDbi::MeSHDb(file_qu)
  
  # de_HI = names(sgnfGenes_HI)[1:100]
  # mesh_HI = enrichMeSH(de_HI, MeSHDb = db, database='gendoo', category = 'C')
  mesh_HI = tryCatch(enrichMeSH(na.omit(entrezIDs_HI$ENTREZID), MeSHDb = meshdb_species, database='gendoo', category = 'C'),
                     error=function(cond) {
                       return(NA)
                     },
                     warning=function(cond) {
                       # return(NA)
                     })
  
  mesh_name_HI = paste(workflow_name,"mesh_HI",sep="_")
  file_name_mesh_HI = paste0(mesh_name_HI, ".txt")
  figure_name_mesh_HI = paste0(mesh_name_HI, ".png")
  write.table(mesh_HI, file = file_name_mesh_HI, sep = "\t", quote = F, row.names = F, na = "")
  
  # de_ME = names(sgnfGenes_ME)[1:100]
  # mesh_ME = enrichMeSH(de_ME, MeSHDb = db, database='gendoo', category = 'C')
  mesh_ME = tryCatch(enrichMeSH(na.omit(entrezIDs_ME$ENTREZID), MeSHDb = meshdb_species, database='gendoo', category = 'C'),
                     error=function(cond) {
                       return(NA)
                     },
                     warning=function(cond) {
                       # return(NA)
                     })
  
  mesh_name_ME = paste(workflow_name,"mesh_ME",sep="_")
  file_name_mesh_ME = paste0(mesh_name_ME, ".txt")
  figure_name_mesh_ME = paste0(mesh_name_ME, ".png")
  write.table(mesh_ME, file = file_name_mesh_ME, sep = "\t", quote = F, row.names = F, na = "")
  
  # plot the results
  if (!is.null(mesh_HI)) {
    png(figure_name_mesh_HI, width = 600, height = 650)
    barplot(mesh_HI)
    dev.off()
  } else {
    print("Failed to return MeSH enrichment analysis results for the given list of HI significant significant genes!")
  }
  
  if (!is.null(mesh_ME)) {
    png(figure_name_mesh_ME, width = 600, height = 650)
    barplot(mesh_ME)
    dev.off()
  } else {
    print("Failed to return MeSH enrichment analysis results for the given list of ME significant significant genes!")
  }
} else
{
  print(paste0("No MeSHDb available for ", my_species))
}
