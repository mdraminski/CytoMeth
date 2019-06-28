
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
# readCondaPath
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
# readPython2Path
#########################################
readPython2Path <- function(file = "python2.info"){
  python2_path <- NA
  if(file.exists(file)){
    python2_path <- data.frame(data.table::fread(file, header = F))
    if(is.null(python2_path))
      stop("File 'python2.info' is empty! Please make sure python 2.* is installed on your system and run: 'whereis python2 > python2.info' in the CytoMeth directory.")
    
    python2_path <- as.character(python2_path[1,])
    python2_path <- python2_path[endsWith(python2_path,'bin/python2')]
    python2_path <- python2_path[1]
  }
  else{
    stop("File 'python2.info' does not exist! Please make sure python 2.* is installed on your system and run: 'whereis python2 > python2.info' in the CytoMeth directory.")
  }
  return(python2_path)
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
    file.path(config$tools_path, config$methratio),
    file.path(config$tools_path, config$trimmomatic),
    file.path(config$tools_path, config$picard),
    file.path(config$tools_path, config$gatk))
  
  for(i in 1:length(cytoMethTools)){
    if(!file.exists(cytoMethTools[i])){
      stop(paste0("Required tool: ",  cytoMethTools[i], 
                  " is not available! Please install Required Tools! (See 'Required Tools' section in CytoMeth manual."))
    }
  }
  cat("All required tools are available. OK.\n")
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
  cat("All required reference files are available. OK.\n")
  return(T)
}

#########################################
# checkIfFileExists
#########################################
checkIfFileExists <- function (file){
  if(!file.exists(file)){
    print(paste0("Error! File: '",file,"' does not exists! Cannot continue the processing!"))
    return(F)
  }else{
    return(T)
  }
}

#########################################
# skipProcess
#########################################
skipProcess <- function(app_name, process_name, subprocess_name, dir_name){
  dir_name <- gsub("//","/",dir_name)
  cat(paste0(app_name,": Process ",process_name, " - '",subprocess_name,"' is skipped. Result exists in the directory: ", dir_name, "\n"))
}

#########################################
# runSystemCommand
#########################################
runSystemCommand <- function(app_name, process_name, subprocess_name, command, verbose = FALSE){
  cat(paste0(app_name, ": Running ", process_name, " - '", subprocess_name, "'...\n"))
  start_time <- Sys.time()
  if(verbose)
    cat(paste0(command,"\n"))
  system(command)
  stop_time <- Sys.time()
  cat(paste0(app_name, ": Process ",process_name," - '",subprocess_name,"' is finished. [", format(stop_time - start_time, digits=3) ,"]\n"))
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

#########################################
# qc2dataframe
#########################################
qc2dataframe <- function(qc_report){
  df <- t(as.data.frame(qc_report))
  df <- data.frame(metric = row.names(df), value = as.character(df[,1]))
  return(df)
}

#########################################
# getVersionR
#########################################
getRVer <- function(){
  rver <- as.numeric(paste0(R.version$major,".",unlist(strsplit(R.version$minor,"\\."))[1]))
  return(rver)
}

###############################
#file.ext
###############################
fileExt <- function(x){
  ext <- unlist(strsplit(basename(x), '[.]'))
  if(length(ext) > 1)
    ext <- tolower(tail(ext, 1))
  else
    ext <- ''
  return (ext)
}
###############################
#open.plot.file
###############################
openPlotFile <- function(filename, width = 10, height = 6, res = 72)
{
  dev.flush()
  ext <- fileExt(filename)
  if (ext == "png") {
    png(filename, width=width, height = height, units = 'in', res = 72)
  } else if (ext == "jpg") {
    # jpg size is set by default
    jpeg(filename, width = width, height = height, units = 'in', res = 72)
  } else if (ext == "pdf") {
    # pdf size is set by default
    pdf(filename, width = width, height = height)
  } else if (ext == "svg") {
    # svg size is set by default
    svg(filename, width = width, height = height)
  } else if (ext == "ps") {
    # ps size is set by default
    postscript(filename, width = width, height = height)
  } else if (ext == "wmf") {
    # wmf size is set by default
    win.metafile(filename, width = width, height = height)
  } else{ # pdf by default
    pdf(filename, width = width, height = height)
  }
  return(T)
}
