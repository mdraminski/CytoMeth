#########################################
# readConfig
#########################################
# file = "conf.yml"; tools.file = "./tools/tools.conf.yml";
readConfig <- function(file = "config.yml", tools.file = "./tools/tools.conf.yml"){
  ### READ MAIN CONFIG
  conf <- yaml::yaml.load_file(file)
  dir.create(file.path(conf$input_path), showWarnings = F)
  dir.create(file.path(conf$results_path), showWarnings = F)
  conf$input_path <- normalizePath(file.path(conf$input_path))
  conf$results_path <- normalizePath(file.path(conf$results_path))
  conf$ref_data_path <- normalizePath(file.path(conf$ref_data_path))

  ### READ TOOLS CONFIG
  conf_tools <- yaml::yaml.load_file(tools.file)
  conf_tools$tools_path <- normalizePath(file.path(conf_tools$tools_path))
  condaPath <- readCondaPath("conda.info")
  if(!is.na(condaPath) & is.null(conf$anaconda_bin_path))
    conf$anaconda_bin_path <- condaPath
  #setup picard parser
  conf_tools$picard_parser <- ""
  if(conf_tools$picard_ver == 2){
    conf_tools$picard_parser <- " -Dpicard.useLegacyParser=true "
  }
  
  conf <- c(conf, conf_tools)
  return(conf)
}

#########################################
# fixMachineConfig
#########################################
#check if the machine can run hardware parameters defined in config.yml [mem_max in GB]
#config <- conf;
fixMachineConfig <- function(config, thread_max = 16, mem_max = 32){
  mem_mult <- 1e+9

  defined_memory <- si2f(config$memory)
  total_ram <- as.numeric(benchmarkme::get_ram())
  opt_memory <- min(round((total_ram * 0.75)), mem_max*mem_mult)
  if(is.na(defined_memory)){
    warning(paste0("Memory detected: ",f2si(total_ram), ". Memory definition in 'config.yml' is broken: ", config$memory, ". Switching to: ", f2si(opt_memory)))
    config$memory <- opt_memory
  }else if(total_ram < si2f(config$memory)){
    warning(paste0("Memory detected: ",f2si(total_ram), " is lower than defined: ", f2si(defined_memory), " in 'config.yml' file. Switching to: ", f2si(opt_memory)))
    config$memory <- opt_memory
  }else if(si2f(config$memory) < (6*mem_mult)){
    warning(paste0("Memory detected: ",f2si(total_ram), ". Memory defined: ", f2si(defined_memory), " in 'config.yml' file is too low. It is highly recommended to use at least 6G for this configuration."))
    config$memory <- defined_memory
  }else{
    config$memory <- defined_memory
  }

  cpu_threads <- benchmarkme::get_cpu()$no_of_cores
  if(cpu_threads < config$threads){
    new_threads <- min(cpu_threads - 1, thread_max)
    warning(paste0("Number of detected CPU threads: ",cpu_threads," is lower than set up threads: ",config$threads," in 'config.yml' file. Switching to: ",new_threads, " threads"))
    config$threads <- new_threads
  }
  return(config)
}

#########################################
# si2f
#########################################
si2f <- function(x){
  sufix <- c('k','kb','m','mb','g','gb','t','tb')
  mult <- c(1e+3,1e+3,1e+6,1e+6,1e+9,1e+9,1e+12,1e+12)
  curr_sufix <- endsWith(tolower(x), sufix)
  if(any(curr_sufix)){
    sufix_length <- nchar(sufix)[curr_sufix]
    x <- as.numeric(stringr::str_trim(substr(x,start = 1, nchar(x)-sufix_length))) * mult[curr_sufix]
  }else{
    x <- as.numeric(x)
  }
  return(x)
}
#########################################
# f2si
#########################################
f2si <- function(x){
  if (x < 1e+06)
    return (paste0(floor(x/1000),'k'))
  else if (x < 1e+09)
    return (paste0(floor(x/1e+06),'M'))
  else if (x < 1e+12)
    return (paste0(floor(x/1e+09),'G'))
  else if (x < 1e+15)
    return (paste0(floor(x/1e+12),'T'))
  else
    return (x)
}
#########################################
# readCondaPath()
#########################################
readCondaPath <- function(file = "conda.info"){
  json_data <- NA
  if(file.exists(file)){
    json_data <- tryCatch({
      json_data <- rjson::fromJSON(paste(readLines("conda.info"), collapse=""))
      json_data <- file.path(json_data$conda_prefix, "bin")
    }, warning = function(w) {
      warning(w)
    }, error = function(e) {
      stop("File 'conda.info' is broken! Please make sure anaconda is installed on your system and run: 'conda info --json > conda.info' in the CytoMeth directory.")
    })
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
  #create tmp report directory
  if(!dir.exists(file.path(results_path, "tmp"))){
    dir.create(file.path(results_path, "tmp"), showWarnings = T)
  }else{
    #print(paste0("Warning! Directory: ",file.path(results_path, "tmp"), " already exists!"))
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
  
  if(!dir.exists(config$ref_data_path)){
    stop(paste0("Required Reference directory: ",  config$ref_data_path, 
                " does not exist. (See 'Reference Files' section in CytoMeth manual."))
  }
  
  cytoMethRefFiles <- c(
    file.path(config$ref_data_path, config$ref_data_sequence_file),
    file.path(config$ref_data_path, paste0(config$ref_data_sequence_file,".fai")),
    file.path(config$ref_data_path, paste0(tools::file_path_sans_ext(config$ref_data_sequence_file),".dict")),
    file.path(config$tools_path, config$trimmomatic_adapter))
  
  config$ref_data_intervals_file <- str_trim(config$ref_data_intervals_file)
  if(str_trim(config$ref_data_intervals_file) != ''){
    cytoMethRefFiles <- c(cytoMethRefFiles, file.path(config$ref_data_path, config$ref_data_intervals_file))
  }
  
  for(i in 1:length(cytoMethRefFiles)){
    if(!file.exists(cytoMethRefFiles[i])){
      stop(paste0("Required Reference File: ",  cytoMethRefFiles[i], 
                  " is not available! Please download Reference Files! (See 'Reference Files' section in CytoMeth manual)."))
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
# checkFormat
#########################################
checkFormat <- function(format, supported = c('bed','rds')){
  format <- tolower(format[1])
  if(!tolower(format) %in% supported){
    warning(paste0("Unsupported format: '", format, "'. Trying to read from 'bed' files."))
    format <- 'bed'
  }
  return(format)
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
# checkLog
#########################################
checkLog <- function(logfile, success_text, process_name){  
  log_text <- readLines(logfile)
    if(!any(grepl(success_text, log_text))){
      print(paste0("Error! Process ", process_name, " is not completed successfully! See the log file: ", logfile))
      print(log_text)
      return(F)
    }else{
      return(T)
    }
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
# makeDataFrame
#########################################
makeDataFrame <- function(data, header){
  if(length(data) != length(header))
    stop("length(data) != length(header)")
  df <- data.frame(matrix(nrow=1, ncol=length(data)))
  colnames(df) <- header
  df[1,] <- data
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
#fileExt
###############################
fileExt <- function(x){
  ext <- unlist(strsplit(basename(x), '[.]'))
  if(length(ext) > 1)
    ext <- tolower(tail(ext, 1))
  else
    ext <- ''
  return (ext)
}

######################################
######## file.remove.if.exists  ####
fileRemoveIfExists <- function(x){
  if(file.exists(x)){
    file.remove(x)
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#########################################
# basename_noext
#########################################
basename_noext <- function(filename){
  filename <- tools::file_path_sans_ext(basename(filename))
  return(filename)
}

#########################################
# basename_sample
#########################################
basename_sample <- function(filename){
  filename <- tools::file_path_sans_ext(basename(filename))
  mymask <- endsWith(tolower(filename),'_r1') | endsWith(tolower(filename),'_r2')
  #remove last three characters if there is _r1 or _r2 at the end
  filename[mymask] <- stringr::str_sub(filename[mymask], end=-4)
  return(filename)
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

#########################################
# qc2dataframe
#########################################
qc2dataframe <- function(qc_report){
  # df <- t(as.data.frame(qc_report))
  # df <- data.frame(metric = row.names(df), value = as.character(df[,1]))
  df <- data.frame(value = unlist(qc_report))
  return(df)
}

###############################
#getMethylDataHeader
###############################
getMethylDataHeader <- function(version = c(1,2), size = 9){
  methHeader <- switch(version[1],
                         c("chr","start","end","strand","context","betaVal","coverage", "numCs", "numTs"),
                         c("chr","start","end","context","betaVal","strand","coverage", "numCs", "numTs", "posCs"))
  return(methHeader[1:size])
}

###############################
#readMethResult
###############################
readMethResult <- function(methyl_result_file, sample_name = NULL, version = c(1,2)){
  methyl_data <- NULL
  if(checkIfFileExists(methyl_result_file)){
    if(endsWith(tolower(methyl_result_file),".rds")){
      methyl_data <- readRDS(methyl_result_file)
    }else{
      methyl_data <- data.table::fread(methyl_result_file, header = F, sep = '\t')
      colnames(methyl_data) <- getMethylDataHeader(version, size = ncol(methyl_data))
    }
    methyl_data <- as.data.frame(methyl_data)
    methyl_data <- methyl_data[,getMethylDataHeader(2, size = ncol(methyl_data))]
  }else{
    methyl_data <- NULL    
  }
  
  if(is.null(sample_name)){
    sample_name <- unlist(str_split(basename(methyl_result_file), "\\."))[1]
  }
  attr(methyl_data, 'sample_name') <- sample_name
  
  return(methyl_data)
}

###############################
# filterMethResult
###############################
# methyl_data <- methyl_result_data
# methyl_data <- methData[[i]]
# ref_sequence_name <- config$ref_control_sequence_name
filterMethResult <- function(methyl_data, ref_sequence_name = NULL, context = c("CG","CHG","CHH"), context_label = NULL, min_coverage = 10, max_coverage = NA){
  methyl_data <- as.data.frame(methyl_data)
  if(!is.null(ref_sequence_name)){
    methyl_data <- methyl_data[methyl_data$chr != ref_sequence_name,]
  }
  if(is.na(max_coverage))
    methyl_data <- methyl_data[methyl_data$context %in% context & methyl_data$coverage >= min_coverage,]
  else
    methyl_data <- methyl_data[methyl_data$context %in% context & methyl_data$coverage < max_coverage,]
  
  if(is.null(context_label)){
    context_label <- paste0(context, collapse = ",")
  }
  attr(methyl_data, "context") <- context_label
  attr(methyl_data, "min_coverage") <- min_coverage
  
  return(methyl_data)
}

###############################
# convertMethResult
###############################
#methyl_data <- methData[[i]]
convertMethResult <- function(methyl_data, assembly='hg38', resolution='base'){
  methyl_data <- new("methylRaw", methyl_data %>% dplyr::select(-context, -betaVal) %>% data.table(), 
      sample.id=attr(methyl_data, 'sample_name'), assembly=assembly, context=attr(methyl_data, "context"), resolution=resolution)
  return(methyl_data)
}

##############################################
######## Read data to MethylKit ############
#conf <- config; context = c("CHG","CHH"); context_label = 'non-CpG'; min_coverage = 10;result_format = c('bed','rds');
readMethData <- function(config, context = c("CG","CHG","CHH"), context_label = NULL, min_coverage = 10, result_format = c('bed','rds')){
  result_format <- checkFormat(result_format, supported = c('bed','rds'))
  result_dir <- file.path(config$results_path, "methyl_results")
  result_files <- list.files(result_dir, full.names = T)
  result_files <- result_files[endsWith(tolower(result_files), paste0('.methylation_results.', result_format))]
  sample_id <- lapply(strsplit(basename(result_files), '\\.'), function(x){x[1]})
  
  if(is.null(context_label)){
    context_label <- paste0(context, collapse = ",")
  }
  
  # read function from methylKit - it does not filter by context 
  # methData <- methylKit::methRead(as.list(result_files),
  #                                 sample.id = sample_id,
  #                                 treatment = c(rep(0, length(sample_id))), # 0/1 control/test samples
  #                                 assembly = "hg38",
  #                                 header = TRUE,
  #                                 context = methcontext,
  #                                 resolution = "base",
  #                                 pipeline = list(fraction=TRUE, chr.col=1, start.col=2, end.col=3, coverage.col=7, strand.col=4, freqC.col=6, context.col=5),
  #                                 mincov = min_coverage)
  # #filtering ref_control_sequence_name
  # methData <- lapply(methData, function(x) {x[x$chr != config$ref_control_sequence_name,]})

  methData <- list()
  for(i in 1:length(result_files)){
    cat(paste0('Reading file: ', result_files[i], '\n'))
    methData[[i]] <- readMethResult(result_files[i], version = 2)
    methData[[i]] <- filterMethResult(methData[[i]], ref_sequence_name = config$ref_control_sequence_name, context, context_label, min_coverage)
    methData[[i]] <- convertMethResult(methData[[i]], assembly='hg38', resolution='base')
  }

  methData <- new("methylRawList", methData, treatment = c(rep(0, length(sample_id))))
  attr(methData, 'context') <- context_label
  attr(methData, 'mincov') <- min_coverage
  
  return(methData)
}

##############################################
######## Read data using MethylKit ############
#config <- conf
getCovSummary <- function(config, min_coverage = c(7,8,9,10,11,12,13), result_format = c('bed','rds')){
  result_format <- checkFormat(result_format, supported = c('bed','rds'))
  result_dir <- file.path(config$results_path, "methyl_results")
  result_files <- list.files(result_dir, full.names = T)
  result_files <- result_files[endsWith(tolower(result_files), paste0('.methylation_results.', result_format))]
  sample_id <- unlist(lapply(strsplit(basename(result_files), '\\.'), function(x){x[1]}))
  covSummaryDF <- list()
  for(i in 1:length(result_files)){
    cat(paste0('Reading file: ', result_files[i], '\n'))
    methyl_result_data <- readMethResult(result_files[i], version = 2)
    cat(paste0('Calculation Coverages: ', result_files[i], '\n'))
    for(c in min_coverage){
      number_of_Cs_in_panel_CpG_cov_minC <- nrow(filterMethResult(methyl_result_data, config$ref_control_sequence_name, context = c('CG'), min_coverage = c))
      number_of_Cs_in_panel_CpG_cov_maxC <- nrow(filterMethResult(methyl_result_data, config$ref_control_sequence_name, context = c('CG'), max_coverage = c))
      number_of_Cs_in_panel_non_CpG_cov_minC <- nrow(filterMethResult(methyl_result_data, config$ref_control_sequence_name, context = c('CHG','CHH'), min_coverage = c))
      
      covSummaryDF[[length(covSummaryDF)+1]] <- data.frame(SampleID = sample_id[i], value = number_of_Cs_in_panel_CpG_cov_minC, context = 'CpG', min_coverage = c, max_coverage = NA)  
      covSummaryDF[[length(covSummaryDF)+1]] <- data.frame(SampleID = sample_id[i], value = number_of_Cs_in_panel_CpG_cov_maxC, context = 'CpG', min_coverage = NA, max_coverage = c)  
      covSummaryDF[[length(covSummaryDF)+1]] <- data.frame(SampleID = sample_id[i], value = number_of_Cs_in_panel_non_CpG_cov_minC, context = 'non-CpG', min_coverage = c, max_coverage = NA) 
    }
  }
  covSummaryDF <- rbindlist(covSummaryDF)
  return (covSummaryDF)
}

#############################################################
######## Read unique chromosomes from sample bed file  #####
getSampleChr <- function(filepath, buff = 100000, chr_filter = F) {
  con <- file(filepath, "r")
  chrlines <- c()
  b <- 0
  while (TRUE) {
    lines <- readLines(con, n = buff)
    if (length(lines) == 0) {
      break
    }
    chrline <- unique(unlist(lapply(str_split(lines, pattern="\t", n=2), function(x) return(x[1]))))
    if(chr_filter){
      chrline <- chrline[startsWith(chrline,"chr")]
    }
    if(length(chrline)>0 | is.null(chrlines)){
      chrlines <- unique(c(chrlines, chrline))
    }
    cat(".")
    b <- b + 1
    if(b%%100 == 0)
      cat("\n")  
  }
  if(isOpen(con, rw = "")){
    close(con)
  }
  cat("\n")
  return(chrlines)
}

#############################################################
######## Read unique chromosomes from reference fa file  ####
getRefChr <- function(filepath, buff = 100000) {
  con <- file(filepath, "r")
  chrlines <- c()
  b <- 0
  while (TRUE) {
    lines <- readLines(con, n = buff)
    if (length(lines) == 0) {
      break
    }
    chrline <- lines[startsWith(lines, ">")]
    if(length(chrline)>0 | is.null(chrlines)){
      chrlines <- unique(c(chrlines, sub(".", "", chrline)))
    }
    cat(".")
    b <- b + 1
    if(b%%100 == 0)
      cat("\n")  
  }
  if(isOpen(con, rw = "")){
    close(con)
  }
  cat("\n")
  return(chrlines)
}

######################################
######## GET BSMAP IndexInterval  ####
getBSMAPIndexInterval <- function(mem_size){
  mem_size <- mem_size/1e+09
  if(mem_size <= 8){
    indexInterval <- 16
  }else if(mem_size <= 12){
    indexInterval <- 8
  }else if(mem_size <= 16){
    indexInterval <- 4
  }else if(mem_size <= 24){
    indexInterval <- 2
  }else{
    indexInterval <- 1
  }
  indexInterval <- paste0("-I ", indexInterval)
  return(indexInterval)
}
