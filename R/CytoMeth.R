source("./R/main.R")
#read default config from the config.yml file
conf <- readConfig()
#check if the machine can run hardware parameters defined in config.yml [mem_max in GB]
conf <- fixMachineConfig(conf, thread_max = 12, mem_max = 24)
#set up required parameters e.g. input path
#conf$input_path <- "./input/"

# run batch processing of all files located in conf$input_path folder
CytoMeth(conf)

# For single sample processing it is required to define R1 and R2 files
# CytoMeth(conf, file.path(conf$input_path,"small_FAKE03_R1.fastq"), file.path(conf$input_path,"small_FAKE03_R2.fastq"))
