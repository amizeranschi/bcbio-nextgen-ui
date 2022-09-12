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
# my_species = sub("_", " ", vep_species)

# organism_type = my_args[4]
# get gtf file location
gtf_location = my_args[4]

path_to_scripts = my_args[5]

# get data
gatk_data = read.table(file, header = F, stringsAsFactors = F)

print(" --- Preparing data...")
# set the column names
colnames(gatk_data) = c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "IMPACT", "DISTANCE", "STRAND", "FLAGS", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE", "CLIN_SIG", "SOMATIC", "PHENO")

# select variants with "HIGH" impact, "MODERATE" impact, or either of them
high_impact_vars = gatk_data[(gatk_data$IMPACT == "HIGH"), ]
moderate_impact_vars = gatk_data[(gatk_data$IMPACT == "MODERATE"), ]
# high_moderate_impact_vars = gatk_data[(gatk_data$IMPACT == "HIGH" | gatk_data$IMPACT == "MODERATE"), ]

# retain the genes corresponding to high- and moderate-impact variants
# genes high-impact variants
high_impact_genes = unique(high_impact_vars[, c("Gene")])
print(" --- genes with HIGH impact variants: ")
print(high_impact_genes)

# genes with moderate-impact variants
moderate_impact_genes = unique(moderate_impact_vars[, c("Gene")])
print(" --- genes with MODERATE impact variants: ")
print(moderate_impact_genes)

# genes with high-impact OR moderate-impact variants
# high_moderate_impact_genes = unique(high_moderate_impact_vars[, c("Gene")])

workflow_name = "variant_calling"

system(paste(paste0("Rscript --vanilla ", path_to_scripts, "/downstreamAnalysis/computeMetrics.R"), high_impact_genes, moderate_impact_genes, vep_species, workflow_name))
