source("./R/main.R")
conf <- readConfig()

# This scrip takes set of sample_QC_summary.yml files as an input
# It creates 1 csv table and 2 figures
CytoMethQCReport(conf)
