
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
    
    BiocManager::install(not_installed, ask = F, update = F)
  }
}

## Install and load the required packages:
packages = c("ade4", "AGHmatrix", "AnnotationDbi", "AnnotationHub", "BGLR", "candisc", "car", "ChIPseeker", "clusterProfiler", "CMplot", "corrgram", "corrplot", "data.table", "detectRUNS", "doParallel", "DOSE", "dplyr", "DT", "emmeans", "enrichplot", "ensembldb", "ensemblVEP", "gdsfmt", "GGally", "ggplot2", "ggrepel", "ggroups", "GOSemSim", "GPArotation", "gtools", "hierfstat", "lsmeans", "MeSHDbi", "meshes", "miscTools", "multcomp", "nnet", "orthopolynom", "outliers", "pcaMethods", "pheatmap", "plotly", "poolr", "psych", "qqman", "qvalue", "raster", "RColorBrewer", "ReactomePA", "readxl", "rehh", "reshape2", "rgdal", "rlist", "rrBLUP", "rstatix", "SNPRelate", "snpStats", "sp", "stringr", "tidyr", "tidyverse", "VariantAnnotation", "vcfR", "vegan", "writexl", "WVPlots")



check_install(packages)
installed = lapply(packages, library, character.only = TRUE)


## if needed, install the argyle package from Github sources
if(!require(argyle))
{
  if(!require(devtools)){install.packages("devtools")}
  library(devtools)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  devtools::install_github("andrewparkermorgan/argyle")
  library(argyle)
}


## if needed, install the IntAssoPlot package from Github sources
if(!require(IntAssoPlot))
{
  if(!require(devtools)){install.packages("devtools")}
  library(devtools)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  devtools::install_github("whweve/IntAssoPlot")
  library(IntAssoPlot)
}


## if needed, install the DRP package from Github sources
if(!require(DRP))
{
  if(!require(devtools)){install.packages("devtools")}
  library(devtools)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  devtools::install_github("camult/DRP")
  library(DRP)
}


## if needed, install the lassosum package from Github sources
if(!require(lassosum))
{
  if(!require(devtools)){install.packages("devtools")}
  library(devtools)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  devtools::install_github("tshmak/lassosum")
  library(lassosum)
}


## if needed, install the bigsnpr package from Github sources
if(!require(bigsnpr))
{
  if(!require(devtools)){install.packages("devtools")}
  library(devtools)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  devtools::install_github("privefl/bigsnpr")
  library(bigsnpr)
}


## if needed, install the synbreed package from Github sources
if(!require(synbreed))
{
  if(!require(devtools)){install.packages("devtools")}
  library(devtools)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  devtools::install_url("https://cran.r-project.org/src/contrib/Archive/regress/regress_1.3-15.tar.gz")
  devtools::install_url("https://cran.r-project.org/src/contrib/Archive/synbreed/synbreed_0.12-9.tar.gz")
  library(synbreed)
}


## if needed, install the synbreedData package from Github sources
if(!require(synbreedData))
{
  if(!require(devtools)){install.packages("devtools")}
  library(devtools)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  devtools::install_url("https://cran.r-project.org/src/contrib/Archive/synbreedData/synbreedData_1.5.tar.gz")
  library(synbreedData)
}


## if needed, install the GAPIT3 package from Github sources
if(!require(GAPIT3))
{
  if(!require(devtools)){install.packages("devtools")}
  library(devtools)
  
  ## allow R to look for pacakges in both CRAN and Bioconductor
  setRepositories(ind = 1:2)
  devtools::install_github("jiabowang/GAPIT3", force=TRUE)
  library(GAPIT3)
}

