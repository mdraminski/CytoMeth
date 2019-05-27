#installation of required R packages by CytoMeth
dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
install.packages(c("yaml"),lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
#install.packages("tools", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages("data.table", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages("rjson", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
install.packages("ggplot2", lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')
