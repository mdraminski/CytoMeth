
#########################################
# readConfig
#########################################
# file = "conf.yml"; tools.file = "./tools/tools.conf.yml";
readConfig <- function(file = "config.yml", tools.file = "./tools/tools.conf.yml"){
  conf <- yaml::yaml.load_file(file)
  dir.create(file.path(conf$input_path), showWarnings = F)
  dir.create(file.path(conf$results_path), showWarnings = F)
  conf$input_path <- normalizePath(file.path(conf$input_path))
  conf$results_path <- normalizePath(file.path(conf$results_path))
  
  conf_tools <- yaml::yaml.load_file(tools.file)
  conf_tools$tools_path <- normalizePath(file.path(conf_tools$tools_path))
  conf_tools$ref_data_path <- normalizePath(file.path(conf_tools$ref_data_path))
  
  condaPath <- readCondaPath("conda.info")
  if(!is.na(condaPath) & is.null(conf$anaconda_bin_path))
    conf$anaconda_bin_path <- condaPath
  
  conf <- c(conf, conf_tools)
  return(conf)
}

#########################################
# createResultDirs
#########################################
readCondaPath <- function(file = "conda.info"){
  json_data <- NA
  if(file.exists(file)){
    json_data <- rjson::fromJSON(paste(readLines("conda.info"), collapse=""))
    json_data <- file.path(json_data$conda_prefix, "bin")
  }
  else{
    stop("File 'conda.info' does not exist! Please make sure anaconda is installed on your system and run: 'conda info --json > conda.info' in the CytoMeth directory.")
  }
  return(json_data)
}

#########################################
# createResultDirs
#########################################
# config = conf;
createResultDirs <- function(config){
  config_tools <- read.csv(file.path(config$tools_path, config$tools_config), stringsAsFactors = FALSE)
  results_path <- config$results_path
  
  #create main results directory
  if(!dir.exists(file.path(results_path))){
    dir.create(file.path(results_path), showWarnings = T)
  }else{
    #print(paste0("Directory: ",file.path(results_path), " already exists!"))
  }
  #create target results directory
  for(i in 1:length(config_tools$temp_results_dirs)){
    if(!dir.exists(file.path(results_path, config_tools$temp_results_dirs[i]))){
      dir.create(file.path(results_path, config_tools$temp_results_dirs[i]), showWarnings = T, recursive = T)
    }
  }
  #create main logs directory
  if(!dir.exists(file.path(results_path, "logs"))){
    dir.create(file.path(results_path, "logs"), showWarnings = T)
  }else{
    #print(paste0("Warning! Directory: ",file.path(results_path,"logs"), " already exists!"))
  }
  #create QC report directory
  if(!dir.exists(file.path(results_path, "QC_report"))){
    dir.create(file.path(results_path, "QC_report"), showWarnings = T)
  }else{
    #print(paste0("Warning! Directory: ",file.path(results_path, "QC_report"), " already exists!"))
  }
  return(T)
}

#########################################
# checkRequiredTools
#########################################
checkRequiredTools <- function(config){
  cytoMethTools <- c(
    file.path(config$anaconda_bin_path, config$bedtools),
    file.path(config$anaconda_bin_path, config$samtools),
    file.path(config$anaconda_bin_path, config$bamtools),
    file.path(config$anaconda_bin_path, config$bamUtil),
    file.path(config$anaconda_bin_path, config$bsmap),
    file.path(config$anaconda_bin_path, config$methratio),
    file.path(config$tools_path, config$trimmomatic),
    file.path(config$tools_path, config$picard),
    file.path(config$tools_path, config$gatk))
  
  for(i in 1:length(cytoMethTools)){
    if(!file.exists(cytoMethTools[i])){
      stop(paste0("Required tool: ",  cytoMethTools[i], 
                  " is not available! Please install Required Tools! (See 'Required Tools' section in CytoMeth manual."))
    }
  }
  print("All required tools are available. OK.")
  return(T)
}

#########################################
# checkRequiredFiles
#########################################
checkRequiredFiles <- function(config){
  cytoMethRefFiles <- c(
    file.path(config$ref_data_path, config$reference_sequence_file),
    file.path(config$ref_data_path, paste0(config$reference_sequence_file,".fai")),
    file.path(config$ref_data_path, paste0(tools::file_path_sans_ext(config$reference_sequence_file),".dict")),
    file.path(config$ref_data_path, config$intervals_file),
    file.path(config$tools_path, config$trimmomatic_adapter))
  
  for(i in 1:length(cytoMethRefFiles)){
    if(!file.exists(cytoMethRefFiles[i])){
      stop(paste0("Required Reference File: ",  cytoMethRefFiles[i], 
                  " is not available! Please download Reference Files! (See 'Reference Files' section in CytoMeth manual."))
    }
  }
  print("All required reference files are available. OK.")
  return(T)
}

#########################################
# checkIfFileExists
#########################################
checkIfFileExists <- function (file){
  if(!file.exists(file))
    stop(paste0("Error! File: '",file,"' does not exists! Cannot continue the processing!"))
}

#########################################
# runSystemCommand
#########################################
runSystemCommand <- function(app_name, process_name, subprocess_name, command){
  print(paste0(app_name, ": Running ", process_name, " - '", subprocess_name, "'..."))
  start_time <- Sys.time()
  system(command)
  stop_time <- Sys.time()
  print(paste0(app_name, ": Process ",process_name," - '",subprocess_name,"' is finished. [", format(stop_time - start_time, digits=3) ,"]"))
}

#########################################
# getLinesNumber
#########################################
getLinesNumber <- function(filepath) {
  con <- file(filepath, "r")
  lines_cnt <- 0
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    lines_cnt <- lines_cnt + 1
  }
  close(con)
  return(lines_cnt)
}
