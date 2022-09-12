
# load the required packages
packages = c("AnnotationHub", "AnnotationDbi", "MeSHDbi", "ggplot2", "clusterProfiler", "DOSE", "meshes", "ChIPseeker", "GOSemSim", "GOplot", "ReactomePA",  "VariantAnnotation", "pheatmap", "RColorBrewer",  "DESeq2", "enrichplot")
installed = lapply(packages, library, character.only = T)

args = commandArgs(trailingOnly=TRUE)
my_args <- head(args)

sgnfGenes_HI = args[1]
sgnfGenes_MO = args[2]
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
  species_keytype = "ENSEMBL"
}



# print(paste0(" --- The current species is ", my_species, " and the annotation key type used is ", species_keytype))


## gene annotation
# use orgdb_species to get a list of Entrez IDs corresponding to the significant gene IDs
entrezIDs_HI = AnnotationDbi::select(orgdb_species, keys = sgnfGenes_HI, keytype = species_keytype, columns = "ENTREZID")
entrezIDs_MO = AnnotationDbi::select(orgdb_species, keys = sgnfGenes_MO, keytype = species_keytype, columns = "ENTREZID")



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
  genenames = lapply(geneIDs, function(x) { AnnotationDbi::mapIds(orgdb_species, keys = unlist(strsplit(x, split = "/")), keytype = orig_keytype, column = "COMMON", multiVals = "first") })
  
  ## replace NA values with original SYMBOL
  for(someiter in 1:length(genenames))
  {
    genenames[[someiter]][is.na(genenames[[someiter]])] = names(genenames[[someiter]][is.na(genenames[[someiter]])])
  }
  
  ## paste together the genenames and return the concatenated strings
  genenames = sapply(genenames, paste, collapse = "/")
}




# use clusterProfiler for representation of gene ontology on each list of interest genes
print(" --- Gene Ontology overrepresentation analysis")

GO_HI = enrichGO(gene          = entrezIDs_HI$ENTREZID,
                 OrgDb         = orgdb_species,
                 keyType       = "ENTREZID",
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
name_HI = paste(workflow_name, "GO_HI", sep = "_")
file_name_HI = paste(name_HI, ".txt")
figure_name_HI = paste(name_HI, ".png")
write.table(GO_HI, file = file_name_HI, sep = "\t", quote = F, row.names = F, na = "")

head(summary(GO_HI), n=3)
head(GO_HI)
GO_MO = enrichGO(gene          = entrezIDs_MO$ENTREZID,
                 OrgDb         = orgdb_species,
                 keyType       = "ENTREZID",
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
name_MO = paste(workflow_name, "GO_MO", sep = "_")
file_name_MO = paste(name_MO, ".txt")
figure_name_MO = paste(name_MO, ".png")
write.table(GO_MO, file = file_name_MO, sep = "\t", quote = F, row.names = F, na = "")

head(GO_MO)





# plot the results
if (!is.null(GO_HI)) {
  png(figure_name_HI, width = 750, height = 750)
  barplot(GO_HI, showCategory=20, title="Enrichment Analysis - representation of gene ontology for genes associated with HIGH impact variants") 
  dev.off()
} else {
  print("Failed to return results of Gene Ontology for HIGH impact variants!")
}
if (!is.null(GO_MO)) {
  png(figure_name_MO, width = 750, height = 950)
  barplot(GO_MO, showCategory=20, title="Enrichment Analysis - representation of gene ontology for genes associated with MODERATE impact variants") 
  dev.off()
} else {
  print("Failed to return results of Gene Ontology for MODERATE impact variants!")
}



# search the organism for kegg analysis
kegg_organism = search_kegg_organism(sub("_", " ", my_species))$kegg_code

print(" --- KEGG pathway over-representation analysis on the extracted genes")

# KEGG pathway over-representation analysis on the extracted genes
kegg_HI = enrichKEGG(gene         = sgnfGenes_HI,
                     organism     = kegg_organism,
                     pvalueCutoff = 0.05)

print(" --- kegg_HI done")

kegg_HI_readable = kegg_HI
kegg_HI_readable@result$geneID = makeReadable(kegg_HI@result$geneID, orig_keytype = species_keytype)
kegg_name_HI = paste(workflow_name,"kegg_HI",sep="_")
file_name_kegg_HI = paste(kegg_name_HI, ".txt")
figure_name_kegg_HI = paste(kegg_name_HI, ".png")
write.table(kegg_HI_readable, file = file_name_kegg_HI, sep = "\t", quote = F, row.names = F, na = "")

# head(kegg_HI)

kegg_MO = enrichKEGG(gene         = sgnfGenes_MO,
                     organism     = kegg_organism,
                     pvalueCutoff = 0.05)

print(" --- kegg_ME done")

kegg_MO_readable = kegg_MO
kegg_MO_readable@result$geneID = makeReadable(kegg_MO@result$geneID, orig_keytype = species_keytype)
kegg_name_MO = paste(workflow_name,"kegg_MO",sep="_")
file_name_kegg_MO = paste0(kegg_name_MO, ".txt")
figure_name_kegg_MO = paste0(kegg_name_MO, ".png")
write.table(kegg_MO_readable, file = file_name_kegg_MO, sep = "\t", quote = F, row.names = F, na = "")

# head(kegg_MO)

print(" --- KEGG pathway over-representation analysis on the extracted genes plots")

# plot the results
if (!is.null(kegg_HI)) {
  png(figure_name_kegg_HI, width = 600, height = 650)
  barplot(kegg_HI, showCategory=20, title="KEGG pathway over-representation analysis on the HI significant genes")
  dev.off()  
} else {
  print("Failed to return results of KEGG analysis for the given list of HI significant significant genes!")
  
}

if (!is.null(kegg_MO)) {
  png(figure_name_kegg_MO, width = 600, height = 650)
  barplot(kegg_MO, showCategory=20, title="KEGG pathway over-representation analysis on the ME significant genes")
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

x_MO = enrichDO(gene                = entrezIDs_MO$ENTREZID,
                ont           = "DO",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = FALSE)
x_name_MO = paste(workflow_name,"disease_ont_MO",sep="_")
file_name_x_MO = paste0(x_name_MO, ".txt")
figure_name_x_MO = paste0(x_name_MO, ".png")
write.table(x_MO, file = file_name_x_MO, sep = "\t", quote = F, row.names = F, na = "")

head(x_MO)

print(" --- Over-representation analysis for disease ontology plots")

# plot the results
if (!is.null(x_HI)) {
  png(figure_name_x_HI, width = 600, height = 650)
  barplot(x_HI, showCategory=20, title="Over-representation analysis for disease ontology in HI significant genes")
  dev.off()
} else {
  print("Failed to return results of Disease Ontology for the given list of HI significant significant genes!")
}

if (!is.null(x_MO)) {
  png(figure_name_x_MO, width = 600, height = 650)
  barplot(x_MO, showCategory=20, title="Over-representation analysis for disease ontology in ME significant genes")
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
  
  # de_MO = names(sgnfGenes_MO)[1:100]
  # mesh_MO = enrichMeSH(de_MO, MeSHDb = db, database='gendoo', category = 'C')
  mesh_MO = tryCatch(enrichMeSH(na.omit(entrezIDs_MO$ENTREZID), MeSHDb = meshdb_species, database='gendoo', category = 'C'),
                     error=function(cond) {
                       return(NA)
                     },
                     warning=function(cond) {
                       # return(NA)
                     })
  
  mesh_name_MO = paste(workflow_name,"mesh_MO",sep="_")
  file_name_mesh_MO = paste0(mesh_name_MO, ".txt")
  figure_name_mesh_MO = paste0(mesh_name_MO, ".png")
  write.table(mesh_MO, file = file_name_mesh_MO, sep = "\t", quote = F, row.names = F, na = "")
  
  # plot the results
  if (!is.na(de_HI)) {
    png(figure_name_mesh_HI, width = 600, height = 650)
    barplot(mesh_HI)
    dev.off()
  } else {
    print("Failed to return MeSH enrichment analysis results for the given list of HI significant significant genes!")
  }
  
  if (!is.na(de_MO)) {
    png(figure_name_mesh_MO, width = 600, height = 650)
    barplot(mesh_MO)
    dev.off()
  } else {
    print("Failed to return MeSH enrichment analysis results for the given list of ME significant significant genes!")
  }
} else
{
  print(paste0("No MeSHDb available for ", my_species))
}
