## The following function tests whether a list of packages is installed and installs them if they are not available:
# install.packages('plyr', repos = "http://cran.us.r-project.org")
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
packages = c("MeSHDbi", "ggplot2","ChIPseeker","clusterProfiler", "DOSE", "meshes", "GOSemSim","ReactomePA", "AnnotationHub", "AnnotationDbi", "VariantAnnotation", "pheatmap", "RColorBrewer",  "DESeq2", "bedr",  "SBGNview", "pathview", "gage", "gageData",  "topGO", "biomaRt")
check_install(packages)
installed = lapply(packages, library, character.only = T)

