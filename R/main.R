options(scipen=999, digits=22)
library("yaml")
library("rjson")
library("tools")
library("data.table")
library("ggplot2")
library("RColorBrewer")
library("methylKit")
library("GenomicRanges")
library("genomation")
library("stringr")
library("benchmarkme")
# library("BSgenome")
# library("Rsamtools")

source("./R/utils.R")
source("./R/main.utils.R")
source("./R/mainQC.R")

#########################################
CytoMethInfo <- function(){
  cat("############################################\n")
  cat("### CytoMeth ver 1.0.3  (April 1st 2023) ###\n")
  cat("############################################\n")
  cat("### Created by Michal Draminski, Agata Dziedzic, Rafal Guzik, Damian Loska, Bartosz Wojtas and Michal J. Dabrowski ###\n")
  cat("### Computational Biology Lab, Polish Academy of Science, Warsaw, Poland ###\n")
  cat("### Neurobiology Center, Nencki Institute of Experimental Biology, Warsaw, Poland ###\n\n")
}

#####################  
#' CytoMeth
#' 
#' Function runs CytoMeth processing.
#' 
#' @param config - configuration CytoMeth file 
#' @param input_file - FASTA R1 file in .fastq format or already mapped .bam file - if NULL then all all samples from /input/ directory will be processed
#' 
#' @return TRUE if succeed
#'
#' @examples
#' conf <- readConfig()
#' CytoMeth(conf)
#' 
#' @import yaml tools data.table rjson ggplot2
#' @export
#' 
#config <- conf
CytoMeth <- function(config, input_file = NULL){
  
  if(!file.exists(config$input_path)){
    stop("Input directory does not exist.")
  }
  
  if(!is.null(input_file)){
    CytoMethSingleSample(config, input_file)
  }else{
    input_files <- list.files(file.path(config$input_path), pattern = "\\.fastq|\\.bam")
    if(length(input_files) == 0 ){
      stop("Input directory does not contain any '*.fastq' or '*.bam' files. [lower case extension is required]")
    }
    input_files <- input_files[grepl("_R1.fastq|_r1.fastq|.bam",input_files)]
    if(length(input_files)==0){
      stop("Input directory does not contain any '*_R1.fastq'/'_R2.fastq or' '*.bam' files.")
    }
    cat(paste0("Input files #",length(input_files)," [", paste0(input_files, collapse = ", "),"]","\n"))
    for(i in 1:length(input_files)){
      sample_name <- gsub("_R1.fastq|_r1.fastq|_R2.fastq|_r2.fastq|.bam", "", input_files[i])
      cat(paste0("Running Processing on Sample: #",i, "/",length(input_files),"\n"))
      cat(paste0("Sample Name: ",sample_name,"\n"))
      input_file <- file.path(config$input_path, input_files[i])
      if(!CytoMethSingleSample(config, input_file))
        print(paste0("Error processing of sample: ", sample_name, "!"))
    }
  }
}

#########################################
# CytoMethSingleSample
#########################################
#config <- conf; config$clean_tmp_files <- F; input_file <-"~/workspace/CytoMeth/input/small_FAKE03_R1.fastq";
#config <- conf; config$clean_tmp_files <- F; input_file <-"~/workspace/CytoMeth/input/DA01_R1.fastq";
CytoMethSingleSample <- function(config, input_file){
  
  #show vignette
  CytoMethInfo()
  #check tools and reference files
  if(!checkRequiredTools(config)) return(F)
  if(!checkRequiredFiles(config)) return(F)
  if(!checkRequiredSettings(config)) return(F)
  
  config$myAppName <- "### CytoMeth ### "
  config_tools <- read.csv(file.path(config$tools_path, config$tools_config), stringsAsFactors = FALSE)
  sample_basename <- basename_sample(input_file)
  full_start_time <- Sys.time()
  
  if(tolower(file_ext(input_file)) == "fastq"){
    if(!endsWith(tolower(basename_noext(input_file)),"_r1")){
      print(paste0("Error! Input 'fastq' file does not contain required 'R1/r1' signature. File '*_R1.fastq' is expected. File: ", input_file))
      return(F)
    }
  }
  
  if(tolower(file_ext(input_file)) == "fastq"){
    config$file_r1 <- input_file
    config$file_r2 <- paste0(substr(input_file,1,nchar(input_file)-nchar("1.fastq")),'2.',file_ext(input_file))
    #check input files
    if(!file.exists(config$file_r1)){
      print(paste0("Error! File ", config$file_r1, " does not exist!"))
      return(F)
    }
    if(!file.exists(config$file_r2)){
      print(paste0("Error! File ", config$file_r2, " does not exist!"))
      return(F)
    }      
  }else if(tolower(file_ext(input_file)) == "bam"){
    config$file_bam <- input_file
    if(!file.exists(config$file_bam)){
      print(paste0("Error! File ", config$file_bam, " does not exist!"))
      return(F)
    }
  }else{
    print(paste0("Error! Input file is not in required format. Files: 'fastq' or 'bam' are required!"))
    return(F)
  }
  
  methyl_result_file <- NULL
  for(i in 1:length(config$meth_tool)){
    methyl_result_file <- c(methyl_result_file, paste0(file.path(config$results_path, config_tools[config_tools$proces==config$meth_tool[i],"temp_results_dirs"], sample_basename),".methylation_results.bed"))
  }
  
  if(!config$overwrite_results & all(file.exists(methyl_result_file))){
    cat('#############################################\n')
    cat(paste0("Sample '",sample_basename,"' is already processed. Skipping the file!\n"))
    cat('#############################################\n')
    ## TODO: add exit codes?
    return(T)
  }
  
  cat('#############################################\n')
  cat(paste0("Processing the sample - '", sample_basename, "'\n"))
  cat('#############################################\n')
  
  #create (if doesnt exist) full folders structure for temp and results files
  if(!createResultDirs(conf)) 
    return(F)
  
  if(tolower(file_ext(input_file)) == "fastq"){
    if(config$sqtk_run){
      cat(paste0("###### 0. Select a Subsample of Reads [seqtk] ######\n"))
      config_ret <- run_seqtk(config, config_tools)
      if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
    }
    
    cat(paste0("###### 1. Examine Sequence Read Quality Control [FastQC] ######\n"))
    config_ret <- run_FastQC(config, config_tools)
    if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
    
    cat(paste0("###### 2. Adapter trimming and Quality Filtering [Trimmomatic] ######\n"))
    config_ret <- run_Trimmomatic(config, config_tools)
    if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
    
    cat(paste0("###### 3. Mapping Reads [BSMAP] ######\n"))
    config_ret <- run_BSMAP(config, config_tools)
    if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
  }
  
  ###############################
  ## start here if input is BAM file
  
  ## Setup config$file_bam if it is not defined yet
  if(!"file_bam" %in% names(config) | tolower(file_ext(input_file)) == "bam"){
    config$file_bam <- file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"], basename(input_file))
    file.copy(from = input_file, to = config$file_bam)
  }
  
  cat(paste0("###### 4. Sorting and Removing Duplicates ######\n"))
  config_ret <- run_RemoveDuplicates(config, config_tools)
  if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
  
  cat(paste0("###### 5. Filter BAM file to keep only mapped and properly paired reads ######\n"))
  config_ret <- run_FilterBAM(config, config_tools)
  if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
  
  cat(paste0("###### 6. Clip Overlapping Reads [bamUtil] ######\n"))
  config_ret <- run_ClipOverlap(config, config_tools)
  if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
  
  cat(paste0("###### 7. Calculate Mapping Metrics [Picard] ######\n"))
  config_ret <- run_MappingMetrics(config, config_tools)
  if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
  
  cat(paste0("###### 8. Calculate Depth of Coverage [GATK] ######\n"))
  config_ret <- run_DepthOfCoverage(config, config_tools)
  if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
  
  methratio_time <- NA
  if('methratio' %in% config$meth_tool) {
    methratio_time <- Sys.time()
    config$meth_tool_current <- "methratio"
    cat(paste0("###### 9a. Determine Methylation Percentage [BSMAP/Methratio] ######\n"))
    config_ret <- run_Methratio(config, config_tools)
    if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
    config_ret <- methratio_to_bed(config, config_tools)
    if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
    config_ret <- run_PanelIntersect(config, config_tools)
    if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
    methratio_time <- format(Sys.time() - methratio_time, digits=3)
  }
  
  bssnper_time <- NA
  if('bssnper' %in% config$meth_tool){
    bssnper_time <- Sys.time()
    config$meth_tool_current <- "bssnper"
    cat(paste0("###### 9b. Determine Methylation Percentage and SNPs [BS-Snper] ######\n"))
    config$bssnper_intersect <- F
    config_ret <- run_BSSnper(config, config_tools)
    if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
    config_ret <- bssnper_to_bed(config, config_tools)
    if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
    config_ret <- run_PanelIntersect(config, config_tools)
    if(is.null(config_ret)) {return(FALSE)} else{config <- config_ret}
    bssnper_time <- format(Sys.time() - bssnper_time, digits=3)
  }
  
  #Remove Temp files
  if(config$clean_tmp_files){
    print("Removing temporary files...")
    files_removed <- sapply(unique(config$tmp_files), fileRemoveIfExists)
    runSystemCommand(config$myAppName, 'rm', 'tmp folder', paste0("rm -rf ", config$tmp_data_path, "/*"), FALSE)
    print("Done.")
  }
  
  cat(paste0("###### 10. Calculate final QC stats ######\n"))
  sampleQC <- getSampleQCSummary(sample_basename, config, result_format = 'bed')
  
  full_stop_time <- Sys.time()
  processing_time <- format(full_stop_time - full_start_time, digits=3)
  sampleQC$processing_time <- processing_time
  sampleQC$methratio_time <- methratio_time
  sampleQC$bssnper_time <- bssnper_time
  print(paste0(sample_basename, " - QC report: "))
  print(qc2dataframe(sampleQC))
  
  # Write QC report to yaml file
  sample_report_file <- file.path(config$results_path, "QC_report", paste0(sample_basename,"_QC_summary.yml"))
  print(paste0("Saving QC report file to: ", sample_report_file))
  yaml::write_yaml(lapply(sampleQC, as.character), sample_report_file)
  
  
  cat('#############################################\n')
  cat(paste0("Processing the sample - '", sample_basename, "' is finished. [", processing_time ,"]\n"))
  cat('#############################################\n\n')
  
  return(T)
}
