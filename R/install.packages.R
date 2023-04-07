source("./R/utils.R")

#installation of required R packages by CytoMeth
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
install.packages("data.table", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages("ggplot2", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages("RColorBrewer", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages("rjson", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages(c("yaml"),lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages(c("stringr"),lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages(c("benchmarkme"),lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')

# install biocunductor packages
bioPckg <- c("Rhtslib","Rsamtools","rtracklayer","methylKit","GenomicRanges","genomation","genomationData")
if(getRVer() < 3.5){
  source("https://bioconductor.org/biocLite.R")
  for(i in 1:length(bioPckg)){
    biocLite(bioPckg[i])
  }
}else{
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  for(i in 1:length(bioPckg)){
    BiocManager::install(bioPckg[i])
  }
}
