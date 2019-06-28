library("yaml")
library("rjson")
library("tools")
library("data.table")
library("ggplot2")
library("RColorBrewer")
library("methylKit")
library("GenomicRanges")
library("genomation")

source("./R/utils.R")
source("./R/mainQC.R")

#########################################
CytoMethInfo <- function(){
  cat("#######################################\n")
  cat("### CytoMeth ver 0.9.8 (27-06-2019) ###\n")
  cat("#######################################\n")
  cat("### Created by Michal Draminski, Agata Dziedzic, Rafal Guzik, Bartosz Wojtas and Michal J. Dabrowski ###\n")
  cat("### Computational Biology Lab, Polish Academy of Science, Warsaw, Poland ###\n")
  cat("### Neurobiology Center, Nencki Institute of Experimental Biology, Warsaw, Poland ###\n\n")
}

#####################  
#' CytoMeth
#' 
#' Function runs CytoMeth processing.
#' 
#' @param config - configuration CytoMeth file 
#' @param file_r1 - FASTA R1 file in .fastq format - if NULL then all all samples from /input/ directory will be processed
#' @param file_r2 - FASTA R2 file in .fastq format - if NULL then all all samples from /input/ directory will be processed
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
CytoMeth <- function(config, file_r1 = NULL, file_r2 = NULL){
  
  if(!file.exists(config$input_path)){
    stop("Input directory does not exist.")
  }
  
  if(!is.null(file_r1) | !is.null(file_r2))
    CytoMethSingleSample(config, file_r1, file_r2)
  
  input_files <- list.files(file.path(config$input_path))
  input_files <- input_files[file_ext(input_files) == "fastq"]
  
  if(length(input_files) == 0 ){
    stop("Input directory does not contain any '*.fastq' files. [fastq ext lower case is required]")
  }
  
  input_files <- gsub("_R1.fastq", "", input_files)
  input_files <- gsub("_R2.fastq", "", input_files)
  input_files <- gsub("_r1.fastq", "", input_files)
  input_files <- gsub("_r2.fastq", "", input_files)
  
  input_files <- unique(input_files)
  
  for(i in 1:length(input_files)){
    file_r1 <- file.path(config$input_path,paste0(input_files[i],"_R1.fastq"))
    file_r2 <- file.path(config$input_path,paste0(input_files[i],"_R2.fastq"))
    if(!CytoMethSingleSample(config, file_r1, file_r2))
      print(paste0("Error processing of sample: ",input_files[i], "!"))
  }
}

#########################################
# CytoMethSingleSample
#########################################
#config <- conf; file_r1 <-"/home/mdraminski/workspace/CytoMeth/input/small_FAKE03_R1.fastq"; file_r2 <- "/home/mdraminski/workspace/CytoMeth/input/small_FAKE03_R2.fastq";
CytoMethSingleSample <- function(config, file_r1, file_r2){
  #show vignette
  CytoMethInfo()
  #check tools and reference files
  if(!checkRequiredTools(config)) return(F)
  if(!checkRequiredFiles(config)) return(F)
  
  #create full folders structure for temp and results files
  #create if they are not exists
  if(!createResultDirs(conf)) return(F)
  
  config_tools <- read.csv(file.path(config$tools_path, config$tools_config), stringsAsFactors = FALSE)

  myAppName <- "### CytoMeth ### "  
  full_start_time <- Sys.time()
  
  if(!file.exists(file_r1)){
    print(paste0("Error! File ",file_r1," does not exist!"))
    return(F)
  }
  if(!file.exists(file_r2)){
    print(paste0("Error! File ",file_r2," does not exist!"))
    return(F)
  }
  
  if(tolower(file_ext(file_r1)) != "fastq" | tolower(file_ext(file_r2)) != "fastq"){
    print(paste0("Error! Input files does not contain required 'R1' or 'R2' signature. Files 'fastq' are expected."))
    return(F)
  }else{
    basename_r1 <- trimws(substr(basename(file_r1), 1, nchar(basename(file_r1)) - nchar(".fastq")) )
    basename_r2 <- trimws(substr(basename(file_r2), 1, nchar(basename(file_r2)) - nchar(".fastq")) )
    
    #define current sample name
    sample_basename <- substr(basename_r1, 1, nchar(basename_r1) - nchar("_R1"))
    
    cat('#############################################\n')
    cat(paste0("Processing the sample - ", sample_basename, "\n"))
    cat('#############################################\n')
    
    ####################
    ######  seqtk ######
    if(!config$sqtk_omit){
      seqtk_r1_file <- paste0(file.path(config$results_path, config_tools[config_tools$proces=="seqtk","temp_results_dirs"], sample_basename),"_R1.fastq")
      seqtk_r2_file <- paste0(file.path(config$results_path, config_tools[config_tools$proces=="seqtk","temp_results_dirs"], sample_basename),"_R2.fastq")
      
      if(config$overwrite_results | !(file.exists(seqtk_r1_file) &
                                      file.exists(seqtk_r2_file))){
        
        src_command <- paste0(file.path(config$anaconda_bin_path, config$seqtk), " sample -s 10000 ", file_r1, " ", config$sqtk_subset, " > ", seqtk_r1_file)
        runSystemCommand(myAppName, 'seqtk', 'subsample of reads R1', src_command, config$verbose)
        
        src_command <- paste0(file.path(config$anaconda_bin_path, config$seqtk), " sample -s 10000 ", file_r2, " ", config$sqtk_subset, " > ", seqtk_r2_file)
        runSystemCommand(myAppName, 'seqtk', 'subsample of reads R2', src_command, config$verbose)
      }else{
        skipProcess(myAppName, 'seqtk', 'subsample of reads',
                    file.path(config$results_path, config_tools[config_tools$proces=="seqtk","temp_results_dirs"],"/"))
      }  
      if(!checkIfFileExists(seqtk_r1_file)) return(F)
      if(!checkIfFileExists(seqtk_r2_file)) return(F)
      file_r1 <- seqtk_r1_file
      file_r2 <- seqtk_r2_file
    }
    
    #####################
    ######  fastqc ######
    fastqc_result_dir <- file.path(config$results_path, config_tools[config_tools$tool=="fastqc","temp_results_dirs"])
    fastqc_result_r1_file <- file.path(fastqc_result_dir, paste0(basename_r1,"_fastqc.html"))
    fastqc_result_r2_file <- file.path(fastqc_result_dir, paste0(basename_r2,"_fastqc.html"))
    
    if(config$overwrite_results | !(file.exists(fastqc_result_r1_file) &
                                    file.exists(fastqc_result_r2_file))){
      # src_command <- paste0(file.path(config$anaconda_bin_path, config$fastqc), " --nogroup ", 
      #                       seqtk_result_file,"_subset_R1.fastq ",
      #                       seqtk_result_file,"_subset_R2.fastq", 
      #                       " --outdir ", fastqc_result_dir)
    
      src_command <- paste0(file.path(config$anaconda_bin_path, config$fastqc), " --nogroup ", 
                            file_r1, " ", file_r2, 
                            " --outdir ", fastqc_result_dir)
      runSystemCommand(myAppName, 'fastqc', 'FastQC Report generation', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Trimmomatic', 'trimming', fastqc_result_dir)
    }
    if(!checkIfFileExists(fastqc_result_r1_file)) return(F)
    if(!checkIfFileExists(fastqc_result_r2_file)) return(F)

    ###############################
    ######  trimmomatic ######
    trimming_result_r1_file <- file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"], basename_r1)
    trimming_result_r2_file <- file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"], basename_r2)
    trimmomatic_out_logfile <- file.path(config$results_path,"logs",paste0(sample_basename,"_",config_tools[config_tools$tool=="trimmomatic","logfile"]))
    
    if(config$overwrite_results | !(file.exists(paste0(trimming_result_r1_file,"_trimmed.fq")) &
                                    file.exists(paste0(trimming_result_r2_file,"_trimmed.fq")) &
                                    file.exists(paste0(trimming_result_r1_file,"_unpaired.fq")) &
                                    file.exists(paste0(trimming_result_r2_file,"_unpaired.fq")) &
                                    file.exists(trimmomatic_out_logfile)
                                    )){
      src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem," -jar ",file.path(config$tools_path,config$trimmomatic),
                            " PE -threads ",config$threads," -phred33 ",file_r1," ", file_r2," ",
                            trimming_result_r1_file,"_trimmed.fq ", trimming_result_r1_file,"_unpaired.fq ",
                            trimming_result_r2_file,"_trimmed.fq ", trimming_result_r2_file,"_unpaired.fq ",
                            "ILLUMINACLIP:", file.path(config$tools_path, config$trimmomatic_adapter),
                            ":2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:",config$trimmomatic_MINLEN,
                            " 2> ",trimmomatic_out_logfile)
      runSystemCommand(myAppName, 'Trimmomatic', 'trimming', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Trimmomatic', 'trimming',
                   file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"],"/"))
    }
    if(!checkIfFileExists(paste0(trimming_result_r1_file,"_trimmed.fq"))) return(F)
    if(!checkIfFileExists(paste0(trimming_result_r2_file,"_trimmed.fq"))) return(F)
    if(!checkIfFileExists(paste0(trimming_result_r1_file,"_unpaired.fq"))) return(F)
    if(!checkIfFileExists(paste0(trimming_result_r2_file,"_unpaired.fq"))) return(F)
    if(!checkIfFileExists(trimmomatic_out_logfile)) return(F)
    
    #Create QC_raport
    sampleQC <- list()
    #Add to QC_raport
    sampleQC$Sample_ID <- sample_basename
    trim_log <- readLines(trimmomatic_out_logfile)
    if(!any(grepl("Completed successfully", trim_log))){
      print(paste0("Error! Trimming phase is not completed successfully! See log file: ", trimmomatic_out_logfile))
      print(trim_log)
      return(F)
    }
    trim_log <- strsplit(trim_log[startsWith(trim_log, "Input Read Pairs")], " ")
    sampleQC$Input_read_pairs <- as.numeric(trim_log[[1]][4])
    sampleQC$Read_Pairs_Surviving_trimming <- as.numeric(trim_log[[1]][7])
    sampleQC$Prc_Read_Pairs_Survaving_trimming <- as.numeric(substr((trim_log[[1]])[8], 2, nchar((trim_log[[1]][8]))-2))
    
    ###############################
    ####### BSMAP - mapping #######
    mapping_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"], sample_basename)
    mapping_dir <- file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"],"/")
    
    if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".sam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bsmap), " -r 0 -s 16 -n 1 ",
                            " -a ", trimming_result_r1_file,"_trimmed.fq",
                            " -b ", trimming_result_r2_file,"_trimmed.fq",
                            " -d ", file.path(config$ref_data_path, config$reference_sequence_file),
                            " -p ", min(config$threads, 8)," -o ", paste0(mapping_result_file, ".sam"))
      runSystemCommand(myAppName, 'BSMAP', 'mapping', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'BSMAP', 'mapping', mapping_dir)
    }
    
    if(!checkIfFileExists(paste0(mapping_result_file, ".sam"))) return(F)
    
    ###############################
    ##### picard - AddOrReplaceReadGroups
    if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".bam")))) {
      src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem, 
                            " -Dpicard.useLegacyParser=false ",
                            " -jar ", file.path(config$tools_path,config$picard), 
                            " AddOrReplaceReadGroups",
                            " -VALIDATION_STRINGENCY LENIENT ",
                            " -INPUT ", paste0(mapping_result_file, ".sam"),
                            " -OUTPUT ", paste0(mapping_result_file, ".bam"),
                            " -CREATE_INDEX true",
                            " -RGID SAMPLE -RGLB SAMPLE -RGPL illumina -RGSM SAMPLE -RGPU platform_unit")
      runSystemCommand(myAppName, 'Picard', 'AddOrReplaceReadGroups', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Picard', 'AddOrReplaceReadGroups', mapping_dir)
    }
    
    if(!checkIfFileExists(paste0(mapping_result_file, ".bam"))) return(F)
    
    ###############################
    #bamtools - split
    if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".TAG_ZS_+-.bam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bamtools), " split -tag ZS",
                            " -in ",paste0(mapping_result_file,".bam"))
      runSystemCommand(myAppName, 'Bamtools', 'split', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Bamtools', 'split', mapping_dir)
    }
    
    if(!checkIfFileExists(paste0(mapping_result_file, ".TAG_ZS_+-.bam"))) return(F)
    
    ###############################
    #bamtools - merge top
    if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".top.bam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " merge",
                            " -in ", paste0(mapping_result_file,".TAG_ZS_++.bam"),
                            " -in ", paste0(mapping_result_file,".TAG_ZS_+-.bam"),
                            " -out ", paste0(mapping_result_file,".top.bam"))
      runSystemCommand(myAppName, 'Bamtools', 'merge top', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Bamtools', 'merge top', mapping_dir)
    }
    
    if(!checkIfFileExists(paste0(mapping_result_file, ".top.bam"))) return(F)
    
    ###############################
    #bamtools - merge bottom
    if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".bottom.bam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bamtools), " merge",
                            " -in ", paste0(mapping_result_file,".TAG_ZS_-+.bam"),
                            " -in ", paste0(mapping_result_file,".TAG_ZS_--.bam"),
                            " -out ", paste0(mapping_result_file,".bottom.bam"))
      runSystemCommand(myAppName, 'Bamtools', 'merge bottom', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Bamtools', 'merge bottom', mapping_dir)
    }
    
    if(!checkIfFileExists(paste0(mapping_result_file, ".bottom.bam"))) return(F)
    
    ###############################
    #bamtools - sort top
    if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".top.bam.sorted")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " sort",
                            " -in ", paste0(mapping_result_file,".top.bam"), 
                            " -out ", paste0(mapping_result_file,".top.bam.sorted"))
      runSystemCommand(myAppName, 'Bamtools', 'sort top', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Bamtools', 'sort top', mapping_dir)
    }
    
    if(!checkIfFileExists(paste0(mapping_result_file, ".top.bam.sorted"))) return(F)
    
    ###############################
    #bamtools - sort bottom
    if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".bottom.bam.sorted")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " sort ",
                            " -in ", paste0(mapping_result_file,".bottom.bam"), 
                            " -out ", paste0(mapping_result_file,".bottom.bam.sorted"))
      runSystemCommand(myAppName, 'Bamtools', 'sort bottom', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Bamtools', 'sort bottom', mapping_dir)
    }
    
    if(!checkIfFileExists(paste0(mapping_result_file, ".bottom.bam.sorted"))) return(F)
    
    ###############################
    #picard - MarkDuplicates - top
    rmdups_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_top","temp_results_dirs"], sample_basename)
    rmdups_result_top_logfile <- file.path(config$results_path, "logs", paste0(sample_basename, "_", config_tools[config_tools$process=="mark_duplicates_top","logfile"]))
    rmdups_dir <- file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_top","temp_results_dirs"],"/")
      
    if(config$overwrite_results | !(file.exists(paste0(rmdups_result_file, ".top.rmdups.bam")) &
                                    file.exists(paste0(rmdups_result_file, ".top.rmdups_metrics.txt")) &
                                    file.exists(rmdups_result_top_logfile)
                                    )) {
      src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem,
                            " -Dpicard.useLegacyParser=false ",
                            " -XX:+UseG1GC -XX:MaxGCPauseMillis=100",
                            " -jar ",file.path(config$tools_path, config$picard), 
                            " MarkDuplicates ",
                            " -VALIDATION_STRINGENCY LENIENT",
                            " -INPUT ", paste0(mapping_result_file,".top.bam.sorted"),
                            " -OUTPUT ", paste0(rmdups_result_file,".top.rmdups.bam"),
                            " -METRICS_FILE ",paste0(rmdups_result_file,".top.rmdups_metrics.txt"),
                            " -REMOVE_DUPLICATES true -ASSUME_SORTED true -CREATE_INDEX true", 
                            " 2> ", rmdups_result_top_logfile)
      runSystemCommand(myAppName, 'Picard', 'MarkDuplicates - top', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Picard', 'MarkDuplicates - top', rmdups_dir)
    }
    
    if(!checkIfFileExists(paste0(rmdups_result_file, ".top.rmdups.bam"))) return(F)
    if(!checkIfFileExists(paste0(rmdups_result_file, ".top.rmdups_metrics.txt"))) return(F)
    if(!checkIfFileExists(rmdups_result_top_logfile)) return(F)
    
    ###############################
    #picard - MarkDuplicates - bottom
    rmdups_result_bottom_logfile <- file.path(config$results_path, "logs", paste0(sample_basename, "_", config_tools[config_tools$process=="mark_duplicates_bottom","logfile"]))
    
    if(config$overwrite_results | !(file.exists(paste0(rmdups_result_file, ".bottom.rmdups.bam")) &
                                    file.exists(paste0(rmdups_result_file, ".bottom.rmdups_metrics.txt")) &
                                    file.exists(rmdups_result_bottom_logfile)
                                    )) {
      src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem,
                            " -XX:+UseG1GC -XX:MaxGCPauseMillis=100",
                            " -Dpicard.useLegacyParser=false ",
                            " -jar ",file.path(config$tools_path,config$picard), 
                            " MarkDuplicates ",
                            " -VALIDATION_STRINGENCY LENIENT",
                            " -INPUT ", paste0(mapping_result_file,".bottom.bam.sorted"),
                            " -OUTPUT ", paste0(rmdups_result_file,".bottom.rmdups.bam"),
                            " -METRICS_FILE ",paste0(rmdups_result_file,".bottom.rmdups_metrics.txt"),
                            " -REMOVE_DUPLICATES true -ASSUME_SORTED true -CREATE_INDEX true", 
                            " 2> ", rmdups_result_bottom_logfile)
      runSystemCommand(myAppName, 'Picard', 'MarkDuplicates - bottom', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Picard', 'MarkDuplicates - bottom',
                   file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_bottom","temp_results_dirs"],"/"))
    }
    
    if(!checkIfFileExists(paste0(rmdups_result_file, ".bottom.rmdups.bam"))) return(F)
    if(!checkIfFileExists(paste0(rmdups_result_file, ".bottom.rmdups_metrics.txt"))) return(F)
    if(!checkIfFileExists(rmdups_result_bottom_logfile)) return(F)
    
    ###############################
    #Add to QC_raport
    top_dup <- readLines(file.path(paste0(rmdups_result_file,".top.rmdups_metrics.txt")))
    top_dup <-  strsplit(top_dup[startsWith(top_dup, "SAMPLE")], "\t")
    sampleQC$Prc_duplicated_reads_top <- as.numeric(top_dup[[1]][9])*100
    
    bottom_dup <- readLines(file.path(paste0(rmdups_result_file,".bottom.rmdups_metrics.txt")))
    bottom_dup <-  strsplit(bottom_dup[startsWith(bottom_dup, "SAMPLE")], "\t")
    sampleQC$Prc_duplicated_reads_bottom <- as.numeric(bottom_dup[[1]][9])*100
    
    ###############################
    #bamtools - merge
    if(config$overwrite_results | !(file.exists(paste0(rmdups_result_file,".rmdups.bam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " merge", 
                            " -in ", paste0(rmdups_result_file,".top.rmdups.bam"),
                            " -in ", paste0(rmdups_result_file,".bottom.rmdups.bam"),
                            " -out ", paste0(rmdups_result_file,".rmdups.bam"))
      runSystemCommand(myAppName, 'Bamtools', 'merge', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Bamtools', 'merge',
                   file.path(config$results_path, config_tools[config_tools$proces=="merge2_bam","temp_results_dirs"],"/"))
    }
    
    if(!checkIfFileExists(paste0(rmdups_result_file,".rmdups.bam"))) return(F)
    
    ###############################
    #bamtools - filter
    filtered_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="filter_bam","temp_results_dirs"], sample_basename)
    
    if(config$overwrite_results | !(file.exists(paste0(filtered_result_file,".filtered.bam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " filter -isMapped true -isPaired true -isProperPair true -forceCompression",
                            " -in ", paste0(rmdups_result_file,".rmdups.bam"), 
                            " -out ", paste0(filtered_result_file,".filtered.bam"))
      runSystemCommand(myAppName, 'Bamtools', 'filter', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Bamtools', 'filter',
                   file.path(config$results_path, config_tools[config_tools$proces=="filter_bam","temp_results_dirs"],"/"))
    }
    
    if(!checkIfFileExists(paste0(filtered_result_file,".filtered.bam"))) return(F)
    
    ###############################
    #Samtools - flagstat flagstat_filtered_bam
    flagstat_result_file <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename,"_",config_tools[config_tools$process=="flagstat_filtered_bam","logfile"]))
    
    if(config$overwrite_results | !(file.exists(paste0(flagstat_result_file)))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$samtools), " flagstat ",  
                            paste0(filtered_result_file,".filtered.bam"), 
                            " 1> ", flagstat_result_file)
      runSystemCommand(myAppName, 'Samtools', 'flagstat', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Samtools', 'flagstat', file.path(config$results_path,"logs","/"))
    }
    
    if(!checkIfFileExists(flagstat_result_file)) return(F)
    
    ###############################
    #bamUtil - clipOverlap
    clipping_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"], sample_basename)
    
    if(config$overwrite_results | !(file.exists(paste0(clipping_result_file,".clipped.bam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bamUtil), " clipOverlap --stats",
                            " --in ",  paste0(filtered_result_file,".filtered.bam"),
                            " --out ", paste0(clipping_result_file,".clipped.bam"), 
                            " 2> ", file.path(file.path(config$results_path,"logs"),paste0(sample_basename,"_",config_tools[config_tools$process=="clip_overlap","logfile"])))
      runSystemCommand(myAppName, 'BamUtil', 'clipOverlap', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'BamUtil', 'clipOverlap',
                   file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"],"/"))
    }
    
    if(!checkIfFileExists(paste0(clipping_result_file,".clipped.bam"))) return(F)
    
    ###############################
    #samtools - index
    if(config$overwrite_results | !(file.exists(paste0(clipping_result_file,".clipped.bam.bai")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$samtools), " index ", paste0(clipping_result_file,".clipped.bam"))
      runSystemCommand(myAppName, 'Samtools', 'index', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Samtools', 'index',
                   file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"],"/"))
    }
    
    if(!checkIfFileExists(paste0(clipping_result_file,".clipped.bam.bai"))) return(F)
    
    ###############################
    #picard - Basic Mapping Metrics
    basic_mapping_metrics_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="basic_mapping_metrics","temp_results_dirs"], sample_basename)
    
    if(config$overwrite_results | !(file.exists(paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt")))) {
      src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem, 
                            " -Dpicard.useLegacyParser=false ",
                            " -jar ", file.path(config$tools_path,config$picard), 
                            " CollectAlignmentSummaryMetrics ",
                            " -VALIDATION_STRINGENCY LENIENT",
                            " -METRIC_ACCUMULATION_LEVEL ALL_READS ",
                            " -INPUT ", paste0(rmdups_result_file,".rmdups.bam"), 
                            " -OUTPUT ", paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt"),
                            " -REFERENCE_SEQUENCE ", file.path(config$ref_data_path, config$reference_sequence_file))
      runSystemCommand(myAppName, 'Picard', 'Basic Mapping Metrics', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Picard', 'Basic Mapping Metrics',
                   file.path(config$results_path, config_tools[config_tools$proces=="basic_mapping_metrics","temp_results_dirs"],"/"))
    }
    
    if(!checkIfFileExists(paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt"))) return(F)
    
    ###############################
    #Add to QC raport
    metrics <- readLines(paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt"))
    metrics <- strsplit(metrics[!startsWith(metrics, "#")],'\t')
    sampleQC$Number_of_reads_after_removing_duplicates <- as.numeric(metrics[[5]][2])
    
    metrics <- readLines(flagstat_result_file)
    metrics <- strsplit(metrics[1]," ")
    sampleQC$Number_of_reads_after_filtering <- as.numeric(metrics[[1]][1])
    sampleQC$Prc_passed_filtering_step <- as.numeric(sampleQC$Number_of_reads_after_filtering)/as.numeric(sampleQC$Number_of_reads_after_removing_duplicates)*100
    
    ###############################
    #picard - Insert Size
    insert_size_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="insert_size","temp_results_dirs"], sample_basename)
    
    if(config$overwrite_results | !(file.exists(paste0(insert_size_result_file,"_insert_size_metrics.txt")))) {
      src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem, 
                            " -Dpicard.useLegacyParser=false ",
                            " -jar ", file.path(config$tools_path,config$picard), 
                            " CollectInsertSizeMetrics ",
                            " -VALIDATION_STRINGENCY LENIENT ",
                            " -Histogram_FILE ", paste0(insert_size_result_file,"_insert_size_plot.pdf"),
                            " -INPUT ", paste0(filtered_result_file,".filtered.bam"),
                            " -OUTPUT ", paste0(insert_size_result_file,"_insert_size_metrics.txt"))
      runSystemCommand(myAppName, 'Picard', 'Insert Size', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'Picard', 'Insert Size',
                   file.path(config$results_path, config_tools[config_tools$proces=="insert_size","temp_results_dirs"],"/"))
    }
    
    if(!checkIfFileExists(paste0(insert_size_result_file,"_insert_size_metrics.txt"))) return(F)
    
    ###############################
    #On-target reads (BEDTools - intersect)
    on_target_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="on_target_reads","temp_results_dirs"], sample_basename)
    
    if(config$overwrite_results | !(file.exists(paste0(on_target_result_file,"_on_target_reads")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bedtools), 
                            " intersect -bed -abam ",  paste0(rmdups_result_file, ".rmdups.bam"), 
                            " -b ", file.path(config$ref_data_path, config$intervals_file), 
                            " > ",paste0(on_target_result_file, "_on_target_reads"))
      runSystemCommand(myAppName, 'BEDTools', 'Intersect on_target_reads', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'BEDTools', 'Intersect on_target_reads',
                   file.path(config$results_path, config_tools[config_tools$proces=="on_target_reads","temp_results_dirs"],"/"))
    }
    
    if(!checkIfFileExists(paste0(on_target_result_file,"_on_target_reads"))) return(F)
    
    #Add to QC raport
    sampleQC$Number_of_on_target_reads <- getLinesNumber(paste0(on_target_result_file,"_on_target_reads"))
    sampleQC$Prc_of_on_target_reads <- sampleQC$Number_of_on_target_read/as.numeric(sampleQC$Number_of_reads_after_removing_duplicates)*100

    ###############################
    #gatk - DepthOfCoverage
    depth_of_cov_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="depth_of_coverage","temp_results_dirs"], sample_basename)
    gatk_depth_logfile <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename, "_", config_tools[config_tools$process=="depth_of_coverage","logfile"]))
    
    if(config$overwrite_results | !(file.exists(paste0(depth_of_cov_result_file,"_gatk_target_coverage")) &
                                    file.exists(paste0(depth_of_cov_result_file,"_gatk_target_coverage.sample_summary")) &
                                    file.exists(gatk_depth_logfile)
                                    )) {
      src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem,
                            " -jar ", file.path(config$tools_path, config$gatk),
                            " -T DepthOfCoverage -R ", file.path(config$ref_data_path, config$reference_sequence_file), 
                            " -I ", paste0(clipping_result_file,".clipped.bam"),
                            " -o ", paste0(depth_of_cov_result_file,"_gatk_target_coverage"),
                            " -L ", file.path(config$ref_data_path, config$intervals_file),
                            " -ct 1 -ct 10 -ct 20", 
                            " > ", gatk_depth_logfile)
      runSystemCommand(myAppName, 'GATK', 'depth_of_coverage', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'GATK', 'depth_of_coverage',
                   file.path(config$results_path, config_tools[config_tools$proces=="depth_of_coverage","temp_results_dirs"],"/"))
    }
    
    if(!checkIfFileExists(paste0(depth_of_cov_result_file,"_gatk_target_coverage"))) return(F)
    if(!checkIfFileExists(paste0(depth_of_cov_result_file,"_gatk_target_coverage.sample_summary"))) return(F)
    if(!checkIfFileExists(gatk_depth_logfile)) return(F)
    
    ###############################
    #Add to QC raport
    depth_summary  <- readLines(paste0(depth_of_cov_result_file,"_gatk_target_coverage.sample_summary"))
    depth_summary  <- strsplit(depth_summary, '\t')
    sampleQC$Mean_coverage <- as.numeric(depth_summary[[2]][3])

    print(qc2dataframe(sampleQC))
    
    ###############################
    #python - methratio (BSMAP) 
    methyl_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"], sample_basename)
    methratio_logfile <- file.path(file.path(config$results_path,"logs"),paste0(sample_basename,"_",config_tools[config_tools$process=="methratio","logfile"]))
    
    if(config$overwrite_results | !(file.exists(paste0(methyl_result_file,".methylation_results.txt")) &
                                    file.exists(methratio_logfile)
                                    )) {
      src_command <- paste0(config$python2, " ", file.path(config$tools_path, config$methratio), 
                            " -d ", file.path(config$ref_data_path, config$reference_sequence_file), 
                            #" -s ", config$anaconda_bin_path,
                            " -m 1 -z -i skip", 
                            " -o ", paste0(methyl_result_file,".methylation_results.txt")," ", paste0(clipping_result_file,".clipped.bam"),
                            " 2> ", methratio_logfile)
      runSystemCommand(myAppName, 'BSMAP', 'methratio', src_command, config$verbose)
    }else{
      skipProcess(myAppName, 'BSMAP', 'methratio',
                   file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"],"/"))
    }

    if(!checkIfFileExists(methratio_logfile)) return(F)
    if(!checkIfFileExists(paste0(methyl_result_file,".methylation_results.txt"))) {
      cat(paste0(readLines(methratio_logfile),"\n"))
      return(F)
    }

    ###############################
    #make bed file from methyl results
    print(paste0("File conversion..."))
    print(paste0("Reading the file: ", methyl_result_file,".methylation_results.txt"))
    methyl_result_data <- data.table::fread(paste0(methyl_result_file,".methylation_results.txt"), header = T, sep = '\t')
    methyl_result_data$end <- methyl_result_data$pos
    methyl_result_data <- methyl_result_data[,c("chr","pos","end","strand","context","ratio","eff_CT_count","C_count", "CT_count")]
    print(paste0("Number of non zero C_count: ", length(methyl_result_data$C_count[methyl_result_data$C_count > 0])))
    print(head(methyl_result_data[methyl_result_data$C_count > 0,]))
    print(paste0("Number of > 1 C_count: ", length(methyl_result_data$C_count[methyl_result_data$C_count > 1])))
    print(head(methyl_result_data[methyl_result_data$C_count > 1,]))
    if(nrow(methyl_result_data)==0){
      print(paste0("Error! File: ", paste0(methyl_result_file,".methylation_results.txt")," is empty!"))
      print(paste0("QC Report of the Sample: ",sample_basename))
      print(qc2dataframe(sampleQC))
      return(F)
    }
    print(paste0("Saving the file: ",methyl_result_file,".methylation_results.bed"))
    data.table::fwrite(methyl_result_data, file=paste0(methyl_result_file,".methylation_results.bed"), quote=FALSE, sep='\t', row.names = F)
    
    ###############################
    #BEDTools - intersect capture region
    
    if(config$overwrite_results | !(file.exists(paste0(methyl_result_file,".methylation_results.bed.panel")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bedtools), 
                            " intersect -bed -abam ",  paste0(methyl_result_file,".methylation_results.bed"), 
                            " -b ",file.path(config$ref_data_path, config$intervals_file),
                            " -u > ",paste0(methyl_result_file,".methylation_results.bed.panel"))
      runSystemCommand(myAppName, 'BEDTools', 'intersect - capture region', src_command, config$verbose)
      
    }else{
      skipProcess(myAppName, 'BEDTools', 'intersect - capture region',
                  file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"],"/"))
    }
    
    if(!checkIfFileExists(paste0(methyl_result_file,".methylation_results.bed.panel"))) return(F)
    
    ###############################
    # Calculate conversion efficiency
    methyl_res_panel_df <- data.table::fread(paste0(methyl_result_file,".methylation_results.bed.panel"), header = F, sep = '\t')
    
    colnames(methyl_res_panel_df) <- c("chr","start","end","strand","context","ratio","eff_CT_count","C_count", "CT_count")
    control <- methyl_res_panel_df[methyl_res_panel_df$chr=='NC_001416',]
    sampleQC$Number_of_Cs_in_control <- nrow(control)
    conversion_eff <- c((1-sum(control$C_count)/sum(control$CT_count))*100)
    sampleQC$Conversion_eff <- conversion_eff
    
    #Add to QC raport
    methyl_res_panel_df_no_control <- methyl_res_panel_df[methyl_res_panel_df$chr!='NC_001416',]
    sampleQC$Number_of_Cs_in_panel <- nrow(methyl_res_panel_df_no_control)
    
    methyl_res_panel_df_no_control_CpG <-methyl_res_panel_df_no_control[methyl_res_panel_df_no_control$context=='CG',]
    sampleQC$Number_of_Cs_in_panel_CpG <- nrow(methyl_res_panel_df_no_control_CpG)
    
    methyl_res_panel_df_no_control_nonCpG <-methyl_res_panel_df_no_control[methyl_res_panel_df_no_control$context!='CG',]
    sampleQC$Number_of_Cs_in_panel_non_CpG <- nrow(methyl_res_panel_df_no_control_nonCpG)
    
    methyl_res_panel_df_no_control_nonCpG_min10 <-methyl_res_panel_df_no_control_nonCpG[methyl_res_panel_df_no_control_nonCpG$eff_CT_count>=10,]
    sampleQC$Number_of_Cs_in_panel_non_CpG_cov_min10 <- nrow(methyl_res_panel_df_no_control_nonCpG_min10)
    
    methyl_res_panel_df_no_control_CpG_min10 <- methyl_res_panel_df_no_control_CpG[methyl_res_panel_df_no_control_CpG$eff_CT_count>=10,]
    sampleQC$Number_of_Cs_in_panel_CpG_cov_min10 <- nrow(methyl_res_panel_df_no_control_CpG_min10)
    
    methyl_res_panel_df_no_control_CpG_max9 <- methyl_res_panel_df_no_control_CpG[methyl_res_panel_df_no_control_CpG$eff_CT_count<10,]
    sampleQC$Number_of_Cs_in_panel_CpG_cov_max9 <- nrow(methyl_res_panel_df_no_control_CpG_max9)
    
    sampleQC$Prc_of_Cs_in_panel_CpG_cov_min10 <-  nrow(methyl_res_panel_df_no_control_CpG_min10)/sum(nrow(methyl_res_panel_df_no_control_CpG_max9)+nrow(methyl_res_panel_df_no_control_CpG_min10))*100
    sampleQC$Prc_of_Cs_in_panel_CpG_cov_max9 <- nrow(methyl_res_panel_df_no_control_CpG_max9)/sum(nrow(methyl_res_panel_df_no_control_CpG_max9)+nrow(methyl_res_panel_df_no_control_CpG_min10))*100
    
    
    #Write yaml raport
    yaml::write_yaml(sampleQC, file.path(config$results_path,"QC_report",paste0(sample_basename,"_QC_summary.yml")))
    # Prepare input for methylkit
    data.table::fwrite(methyl_res_panel_df_no_control_CpG_min10, file=paste0(methyl_result_file,".methylation_results.bed.panel.no_control.CpG_min10"), quote=FALSE, sep='\t', row.names = F)
    data.table::fwrite(methyl_res_panel_df_no_control_nonCpG_min10, file=paste0(methyl_result_file,".methylation_results.bed.panel.no_control.non_CpG_min10"), quote=FALSE, sep='\t', row.names = F)
    
    #Remove Temp files
    if(config$clean_tmp_files){
      print("Removing temporary files...")
      files2remove <- c(paste0(trimming_result_r1_file, c("_trimmed.fq","_unpaired.fq")),
                        paste0(trimming_result_r2_file, c("_trimmed.fq","_unpaired.fq")),
                        paste0(mapping_result_file, c(".sam",".bam",".top.bam",".bottom.bam",".TAG_ZS_++.bam",".TAG_ZS_+-.bam",".TAG_ZS_-+.bam",".TAG_ZS_--.bam",".top.bam.sorted",".bottom.bam.sorted")),
                        paste0(rmdups_result_file, c(".top.rmdups.bam",".bottom.rmdups.bam",".top.rmdups.bai",".bottom.rmdups.bai",".rmdups.bam")),
                        paste0(filtered_result_file,".filtered.bam"), paste0(on_target_result_file,"_on_target_reads"),
                        paste0(methyl_result_file, c(".methylation_results.txt",".methylation_results.bed")),
                        paste0(depth_of_cov_result_file, c("_gatk_target_coverage","_gatk_target_coverage.sample_cumulative_coverage_counts",
                                                           "_gatk_target_coverage.sample_interval_summary","_gatk_target_coverage.sample_cumulative_coverage_proportions",
                                                           "_gatk_target_coverage.sample_interval_statistics")))
      files_removed <- sapply(files2remove, file.remove)
      print("Done.")
    }
    
    full_stop_time <- Sys.time()
    cat('#############################################\n')
    cat(paste0("Processing the sample - '", sample_basename, "' is finished. [", format(full_stop_time - full_start_time, digits=3) ,"]\n"))
    cat('#############################################\n\n')
  }
  
  return(T)
}
