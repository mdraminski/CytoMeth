source("./R/main.R")
conf <- readConfig()

# batch processing of all files located in /input/ directory
CytoMeth(conf)

# for single sample processing you need to define R1 and R2 files
#CytoMeth(conf, file.path(conf$input_path,"small_FAKE01_R1.fastq"), file.path(conf$input_path,"small_FAKE01_R2.fastq"))
#CytoMeth(conf, file.path(conf$input_path,"small_FAKE02_R1.fastq"), file.path(conf$input_path,"small_FAKE02_R2.fastq"))
#CytoMeth(conf, file.path(conf$input_path,"small_FAKE03_R1.fastq"), file.path(conf$input_path,"small_FAKE03_R2.fastq"))
