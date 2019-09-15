source("./R/utils.R")

#installation of required R packages by CytoMeth
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
install.packages("data.table", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages("ggplot2", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages("RColorBrewer", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages("rjson", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages(c("yaml"),lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages(c("stringr"),lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')

# install biocunductor packages
bioPckg <- c("methylKit","GenomicRanges","genomation","genomationData")
if(getRVer() < 3.5){
  source("https://bioconductor.org/biocLite.R")
  biocLite(bioPckg)
}else{
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install(bioPckg)
}
