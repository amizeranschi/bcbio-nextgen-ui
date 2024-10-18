## set up a default CRAN mirror
chooseCRANmirror(ind = 1)

## The following function tests whether a list of packages is installed and installs them if they are not available:
check_install = function(packages) 
{
  not_installed = setdiff(packages, rownames(installed.packages()))
  if(length(not_installed) > 0) 
  {
    write(paste("The libraries", not_installed, "are not available, so they are being installed now",
                sep=" "), stdout())
    
    ## use BiocManager::install() instead of install.packages()
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install(not_installed, quiet = T, ask = F, update = F)
  }
}

## Install and load the required packages:
packages = c("ade4", "AGHmatrix", "AHMeSHDbs", "AnnotationDbi", "AnnotationHub", "bedr", "BGLR", "candisc", "car", "ChIPseeker", "clusterProfiler", "CMplot", "corrgram", "corrplot", "data.table", "DESeq2", "DO.db", "doParallel", "DOSE", "dplyr", "DT", "emmeans", "enrichplot", "ensembldb", "gdsfmt", "GeneOverlap", "GENESIS", "GGally", "ggplot2", "ggrepel", "ggroups", "ggupset", "GOSemSim", "GPArotation", "gtools", "hierfstat", "karyoploteR", "lsmeans", "MeSHDbi", "meshes", "miscTools", "multcomp", "nnet", "orthopolynom", "outliers", "pathview", "pcaMethods", "pheatmap", "plotly", "poolr", "psych", "qqman", "qvalue", "raster", "RColorBrewer", "ReactomePA", "readxl", "rehh", "reshape2", "rgdal", "rlist", "rrBLUP", "rstatix", "SBGNview", "SNPRelate", "snpStats", "sp", "stringr", "tidyr", "tidyverse", "VariantAnnotation", "vcfR", "vegan", "writexl")


check_install(packages)
installed = lapply(packages, library, character.only = TRUE)


## if needed, install the argyle package from Github sources
if(!require(argyle))
{
  if(!require(remotes)){install.packages("remotes")}
  library(remotes)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  remotes::install_github("andrewparkermorgan/argyle", force = T, upgrade = F)
  library(argyle)
}


## if needed, install the IntAssoPlot package from Github sources
if(!require(IntAssoPlot))
{
  if(!require(remotes)){install.packages("remotes")}
  library(remotes)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  remotes::install_github("whweve/IntAssoPlot", force = T, upgrade = F)
  library(IntAssoPlot)
}


## if needed, install the DRP package from Github sources
if(!require(DRP))
{
  if(!require(remotes)){install.packages("remotes")}
  library(remotes)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  remotes::install_github("camult/DRP", force = T, upgrade = F)
  library(DRP)
}


## if needed, install the lassosum package from Github sources
if(!require(lassosum))
{
  if(!require(remotes)){install.packages("remotes")}
  library(remotes)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  remotes::install_github("tshmak/lassosum", force = T, upgrade = F)
  library(lassosum)
}


## if needed, install the bigsnpr package from Github sources
if(!require(bigsnpr))
{
  if(!require(remotes)){install.packages("remotes")}
  library(remotes)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  remotes::install_github("privefl/bigsnpr", force = T, upgrade = F)
  library(bigsnpr)
}


## if needed, install the synbreed package from Github sources
if(!require(synbreed))
{
  if(!require(remotes)){install.packages("remotes")}
  library(remotes)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  remotes::install_url("https://cran.r-project.org/src/contrib/Archive/regress/regress_1.3-15.tar.gz")
  remotes::install_url("https://cran.r-project.org/src/contrib/Archive/synbreed/synbreed_0.12-9.tar.gz")
  library(synbreed)
}


## if needed, install the synbreedData package from Github sources
if(!require(synbreedData))
{
  if(!require(remotes)){install.packages("remotes")}
  library(remotes)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  remotes::install_url("https://cran.r-project.org/src/contrib/Archive/synbreedData/synbreedData_1.5.tar.gz")
  library(synbreedData)
}


## if needed, install the DRP package from Github sources
if(!require(GAPIT)){
  if(!require(remotes)){install.packages("remotes")}
  library(remotes)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  if(!require(LDheatmap)){
    install.packages("https://cran.r-project.org/src/contrib/Archive/LDheatmap/LDheatmap_1.0-6.tar.gz")
  }
  remotes::install_github("jiabowang/GAPIT@GAPIT3.4", upgrade = "never")
  library(GAPIT)
}

