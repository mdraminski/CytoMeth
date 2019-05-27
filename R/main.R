library("yaml")
library("rjson")
library("tools")
library("data.table")
library("ggplot2")

source("./R/utils.R")
source("./R/mainQC.R")

#########################################
CytoMethInfo <- function(){
  cat("#######################################\n")
  cat("### CytoMeth ver 0.9.2 (19-05-2019) ###\n")
  cat("#######################################\n")
  cat("### Created by Michal J. Dabrowski, Agata Dziedzic and Michal Draminski ###\n")
  cat("### Computational Biology Lab - IPI PAN Warsaw, Poland ###\n\n")
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
  i=1
  for(i in 1:length(input_files)){
    file_r1 <- file.path(config$input_path,paste0(input_files[i],"_R1.fastq"))
    file_r2 <- file.path(config$input_path,paste0(input_files[i],"_R2.fastq"))
    CytoMethSingleSample(config, file_r1, file_r2)
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
  checkRequiredTools(config)
  checkRequiredFiles(config)
  
  #create full folders structure for temp and results files
  #create if they are not exists
  createResultDirs(conf)
  
  config_tools <- read.csv(file.path(config$tools_path, config$tools_config), stringsAsFactors = FALSE)

  myAppName <- "### CytoMeth ### "  
  full_start_time <- Sys.time()
  
  if(!file.exists(file_r1)){
    stop(paste0("File ",file_r1," does not exist!"))
  }
  if(!file.exists(file_r2)){
    stop(paste0("File ",file_r2," does not exist!"))
  }
  
  if(tolower(file_ext(file_r1)) != "fastq" | tolower(file_ext(file_r2)) != "fastq"){
    stop(paste0("Input files does not contain required 'R1' or 'R2' signature. Files 'fastq' are expected."))
  }else{
    basename_r1 <- trimws( substr( basename(file_r1), 1, nchar(basename(file_r1)) - nchar(".fastq")) )
    basename_r2 <- trimws( substr( basename(file_r2), 1, nchar(basename(file_r2)) - nchar(".fastq")) )
    
    #define current sample name
    sample_basename <- substr( basename_r1, 1, nchar(basename_r1) - nchar("_R1"))
    
    print('#############################################')
    print(paste0("Processing the sample - ", sample_basename, ""))
    print('#############################################')
    
    ###############################
    ######  trimmomatic ######
    trimming_result_r1_file <- file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"], basename_r1)
    trimming_result_r2_file <- file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"], basename_r2)
    trimmomatic_out_logfile <- file.path(config$results_path,"logs",paste0(sample_basename,"_",config_tools[config_tools$tool=="trimmomatic","logfile"]))
    
    if(config$overwrite_results || !(file.exists(paste0(trimming_result_r1_file,"_trimmed.fq")) &&
                                     file.exists(paste0(trimming_result_r2_file,"_trimmed.fq")) &&
                                     file.exists(paste0(trimming_result_r1_file,"_unpaired.fq")) &&
                                     file.exists(paste0(trimming_result_r2_file,"_unpaired.fq")))){
      
      src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem," -jar ",file.path(config$tools_path,config$trimmomatic),
                            " PE -threads ",config$threads," -phred33 ",file_r1," ", file_r2," ",
                            trimming_result_r1_file,"_trimmed.fq ", trimming_result_r1_file,"_unpaired.fq ",
                            trimming_result_r2_file,"_trimmed.fq ", trimming_result_r2_file,"_unpaired.fq ",
                            "ILLUMINACLIP:", file.path(config$tools_path, config$trimmomatic_adapter),":2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50",
                            " 2> ",trimmomatic_out_logfile)
      runSystemCommand(myAppName, 'Trimmomatic', 'trimming', src_command)
    }else{
      print(paste0(myName, "Process 'trimming' is skipped. Results files are present in directory: ",file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"],"/")))
    }
    checkIfFileExists(paste0(trimming_result_r1_file,"_trimmed.fq"))
    checkIfFileExists(paste0(trimming_result_r2_file,"_trimmed.fq"))
    checkIfFileExists(paste0(trimming_result_r1_file,"_unpaired.fq"))
    checkIfFileExists(paste0(trimming_result_r2_file,"_unpaired.fq"))
    checkIfFileExists(trimmomatic_out_logfile)
    
    #Add to QC_raport
    sample_qc_list <- list()
    sample_qc_list$Sample_ID <- sample_basename
    trim_log <- readLines(trimmomatic_out_logfile)
    trim_log <- strsplit(trim_log[startsWith(trim_log, "Input Read Pairs")], " ")
    sample_qc_list$Input_read_pairs <- as.numeric(trim_log[[1]][4])
    sample_qc_list$Read_Pairs_Surviving_trimming <- as.numeric(trim_log[[1]][7])
    sample_qc_list$Prc_Read_Pairs_Survaving_trimming <- as.numeric(substr((trim_log[[1]])[8], 2, nchar((trim_log[[1]][8]))-2))
    
    ###############################
    ####### BSMAP - mapping #######
    mapping_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"], sample_basename)
    if(config$overwrite_results || !(file.exists(paste0(mapping_result_file, ".sam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bsmap), " -r 0 -s 16 -n 1 ",
                            " -a ", trimming_result_r1_file,"_trimmed.fq",
                            " -b ", trimming_result_r2_file,"_trimmed.fq",
                            " -d ", file.path(config$ref_data_path, config$reference_sequence_file),
                            " -p ", min(config$threads, 8)," -o ", paste0(mapping_result_file, ".sam"))
      runSystemCommand(myAppName, 'BSMAP', 'mapping', src_command)
    }else{
      print(paste0(myName, "Process BSMAP 'mapping' is skipped. Results files are present in the directory: ", file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"],"/")))
    }
    
    checkIfFileExists(paste0(mapping_result_file, ".sam"))
    
    ###############################
    ##### picard - AddOrReplaceReadGroups
    if(config$overwrite_results || !(file.exists(paste0(mapping_result_file, ".bam")))) {
      src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem, 
                            " -Dpicard.useLegacyParser=false ",
                            " -jar ", file.path(config$tools_path,config$picard), 
                            " AddOrReplaceReadGroups",
                            " -VALIDATION_STRINGENCY LENIENT ",
                            " -INPUT ", paste0(mapping_result_file, ".sam"),
                            " -OUTPUT ", paste0(mapping_result_file, ".bam"),
                            " -RGID SAMPLE -RGLB SAMPLE -RGPL illumina -RGSM SAMPLE -RGPU platform_unit")
      runSystemCommand(myAppName, 'Picard', 'AddOrReplaceReadGroups', src_command)
    }else{
      print(paste0(myName, "Process 'AddOrReplaceReadGroups' is skipped. Results files are present in directory: ",file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"],"/")))
    }
    checkIfFileExists(paste0(mapping_result_file, ".bam"))
    
    ###############################
    #bamtools - split
    if(config$overwrite_results || !(file.exists(paste0(mapping_result_file, ".TAG_ZS_+-.bam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bamtools), " split -tag ZS",
                            " -in ",paste0(mapping_result_file,".bam"))
      runSystemCommand(myAppName, 'Bamtools', 'split', src_command)
    }else{
      print(paste0("Process 'split' is skipped. Results files are present in directory: ",file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"],"/")))
    }
    checkIfFileExists(paste0(mapping_result_file, ".TAG_ZS_+-.bam"))
    
    ###############################
    #bamtools - merge top
    if(config$overwrite_results || !(file.exists(paste0(mapping_result_file, ".top.bam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " merge",
                            " -in ", paste0(mapping_result_file,".TAG_ZS_++.bam"),
                            " -in ", paste0(mapping_result_file,".TAG_ZS_+-.bam"),
                            " -out ", paste0(mapping_result_file,".top.bam"))
      runSystemCommand(myAppName, 'Bamtools', 'merge top', src_command)
    }else{
      print(paste0("Process 'merge top' is skipped. Results files are present in directory: ",file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"],"/")))
    }
    checkIfFileExists(paste0(mapping_result_file, ".top.bam"))
    
    ###############################
    #bamtools - merge bottom
    if(config$overwrite_results || !(file.exists(paste0(mapping_result_file, ".bottom.bam")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bamtools), " merge",
                            " -in ", paste0(mapping_result_file,".TAG_ZS_-+.bam"),
                            " -in ", paste0(mapping_result_file,".TAG_ZS_--.bam"),
                            " -out ", paste0(mapping_result_file,".bottom.bam"))
      runSystemCommand(myAppName, 'Bamtools', 'merge bottom', src_command)
    }else{
      print(paste0("Process 'merge bottom' is skipped. Results files are present in directory: ",file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"],"/")))
    }
    checkIfFileExists(paste0(mapping_result_file, ".bottom.bam"))
    
    ###############################
    #bamtools - sort top
    if(config$overwrite_results || !(file.exists(paste0(mapping_result_file, ".top.bam.sorted")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " sort",
                            " -in ", paste0(mapping_result_file,".top.bam"), 
                            " -out ", paste0(mapping_result_file,".top.bam.sorted"))
      runSystemCommand(myAppName, 'Bamtools', 'sort top', src_command)
    }else{
      print(paste0("Process 'sort top' is skipped. Results files are present in directory: ",file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"],"/")))
    }
    checkIfFileExists(paste0(mapping_result_file, ".top.bam.sorted"))
    
    ###############################
    #bamtools - sort bottom
    if(config$overwrite_results || !(file.exists(paste0(mapping_result_file, ".bottom.bam.sorted")))) {
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " sort ",
                            " -in ", paste0(mapping_result_file,".bottom.bam"), 
                            " -out ", paste0(mapping_result_file,".bottom.bam.sorted"))
      runSystemCommand(myAppName, 'Bamtools', 'sort bottom', src_command)
    }else{
      print(paste0("Process 'sort bottom' is skipped. Results files are present in directory: ",file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"],"/")))
    }
    checkIfFileExists(paste0(mapping_result_file, ".bottom.bam.sorted"))
    
    ###############################
    #picard - MarkDuplicates - top
    rmdups_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_top","temp_results_dirs"], sample_basename)
    rmdups_result_top_logfile <- file.path(config$results_path, "logs", paste0(sample_basename, "_", config_tools[config_tools$process=="mark_duplicates_top","logfile"]))
    
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
    runSystemCommand(myAppName, 'Picard', 'MarkDuplicates - top', src_command)
    
    checkIfFileExists(paste0(rmdups_result_file, ".top.rmdups.bam"))
    checkIfFileExists(paste0(rmdups_result_file, ".top.rmdups_metrics.txt"))
    checkIfFileExists(rmdups_result_top_logfile)
    
    ###############################
    #picard - MarkDuplicates - bottom
    rmdups_result_bottom_logfile <- file.path(config$results_path, "logs", paste0(sample_basename, "_", config_tools[config_tools$process=="mark_duplicates_bottom","logfile"]))
    
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
    runSystemCommand(myAppName, 'Picard', 'MarkDuplicates - bottom', src_command)
    
    checkIfFileExists(paste0(rmdups_result_file, ".bottom.rmdups.bam"))
    checkIfFileExists(paste0(rmdups_result_file, ".bottom.rmdups_metrics.txt"))
    checkIfFileExists(rmdups_result_bottom_logfile)
    
    ###############################
    #Add to QC_raport
    top_dup <- readLines(file.path(paste0(rmdups_result_file,".top.rmdups_metrics.txt")))
    top_dup <-  strsplit(top_dup[startsWith(top_dup, "SAMPLE")], "\t")
    sample_qc_list$Prc_duplicated_reads_top <- as.numeric(top_dup[[1]][9])*100
    
    bottom_dup <- readLines(file.path(paste0(rmdups_result_file,".bottom.rmdups_metrics.txt")))
    bottom_dup <-  strsplit(bottom_dup[startsWith(bottom_dup, "SAMPLE")], "\t")
    sample_qc_list$Prc_duplicated_reads_bottom <- as.numeric(bottom_dup[[1]][9])*100
    
    ###############################
    #bamtools - merge
    src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " merge", 
                          " -in ", paste0(rmdups_result_file,".top.rmdups.bam"),
                          " -in ", paste0(rmdups_result_file,".bottom.rmdups.bam"),
                          " -out ", paste0(rmdups_result_file,".rmdups.bam"))
    runSystemCommand(myAppName, 'Bamtools', 'merge', src_command)
    
    checkIfFileExists(paste0(rmdups_result_file,".rmdups.bam"))
    
    ###############################
    #bamtools - filter
    filtered_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="filter_bam","temp_results_dirs"], sample_basename)
    src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " filter -isMapped true -isPaired true -isProperPair true -forceCompression",
                          " -in ", paste0(rmdups_result_file,".rmdups.bam"), 
                          " -out ", paste0(filtered_result_file,".filtered.bam"))
    runSystemCommand(myAppName, 'Bamtools', 'filter', src_command)
    
    checkIfFileExists(paste0(filtered_result_file,".filtered.bam"))
    
    ###############################
    #Samtools - flagstat flagstat_filtered_bam
    flagstat_result_file <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename,"_",config_tools[config_tools$process=="flagstat_filtered_bam","logfile"]))
    src_command <- paste0(file.path(config$anaconda_bin_path, config$samtools), " flagstat ",  
                          paste0(filtered_result_file,".filtered.bam"), 
                          " 1> ", flagstat_result_file)
    runSystemCommand(myAppName, 'Samtools', 'flagstat', src_command)
    
    checkIfFileExists(flagstat_result_file)
    
    ###############################
    #bamUtil - clipOverlap
    clipping_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"], sample_basename)
    src_command <- paste0(file.path(config$anaconda_bin_path, config$bamUtil), " clipOverlap --stats",
                          " --in ",  paste0(filtered_result_file,".filtered.bam"),
                          " --out ", paste0(clipping_result_file,".clipped.bam"), 
                          " 2> ", file.path(file.path(config$results_path,"logs"),paste0(sample_basename,"_",config_tools[config_tools$process=="clip_overlap","logfile"])))
    runSystemCommand(myAppName, 'BamUtil', 'clipOverlap', src_command)
    
    checkIfFileExists(paste0(clipping_result_file,".clipped.bam"))
    
    ###############################
    #samtools - index
    src_command <- paste0(file.path(config$anaconda_bin_path, config$samtools), " index ", paste0(clipping_result_file,".clipped.bam"))
    runSystemCommand(myAppName, 'Samtools', 'index', src_command)
    
    checkIfFileExists(paste0(clipping_result_file,".clipped.bam.bai"))
    
    ###############################
    #picard - Basic Mapping Metrics
    basic_mapping_metrics_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="basic_mapping_metrics","temp_results_dirs"], sample_basename)
    src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem, 
                          " -Dpicard.useLegacyParser=false ",
                          " -jar ", file.path(config$tools_path,config$picard), 
                          " CollectAlignmentSummaryMetrics ",
                          " -VALIDATION_STRINGENCY LENIENT",
                          " -METRIC_ACCUMULATION_LEVEL ALL_READS ",
                          " -INPUT ", paste0(rmdups_result_file,".rmdups.bam"), 
                          " -OUTPUT ", paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt"),
                          " -REFERENCE_SEQUENCE ", file.path(config$ref_data_path, config$reference_sequence_file))
    runSystemCommand(myAppName, 'Picard', 'Basic Mapping Metrics', src_command)
    
    checkIfFileExists(paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt"))
    
    ###############################
    #Add to QC raport
    metrics <- readLines(paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt"))
    metrics <- strsplit(metrics[!startsWith(metrics, "#")],'\t')
    sample_qc_list$Number_of_reads_after_removing_duplicates <- as.numeric(metrics[[5]][2])
    
    metrics <- readLines(flagstat_result_file)
    metrics <- strsplit(metrics[1]," ")
    sample_qc_list$Number_of_reads_after_filtering <- as.numeric(metrics[[1]][1])
    sample_qc_list$Prc_passed_filtering_step <- as.numeric(sample_qc_list$Number_of_reads_after_filtering)/as.numeric(sample_qc_list$Number_of_reads_after_removing_duplicates)*100
    
    ###############################
    #picard - Insert Size
    insert_size_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="insert_size","temp_results_dirs"], sample_basename)
    src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem, 
                          " -Dpicard.useLegacyParser=false ",
                          " -jar ", file.path(config$tools_path,config$picard), 
                          " CollectInsertSizeMetrics ",
                          " -VALIDATION_STRINGENCY LENIENT ",
                          " -Histogram_FILE ", paste0(insert_size_result_file,"_insert_size_plot.pdf"),
                          " -INPUT ", paste0(filtered_result_file,".filtered.bam"),
                          " -OUTPUT ", paste0(insert_size_result_file,"_insert_size_metrics.txt"))
    runSystemCommand(myAppName, 'Picard', 'Insert Size', src_command)
    
    checkIfFileExists(paste0(insert_size_result_file,"_insert_size_metrics.txt"))
    
    ###############################
    #On-target reads (BEDTools - intersect)
    on_target_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="on_target_reads","temp_results_dirs"], sample_basename)
    
    src_command <- paste0(file.path(config$anaconda_bin_path, config$bedtools), 
                          " intersect -nonamecheck -bed -abam ",  paste0(rmdups_result_file, ".rmdups.bam"),
                          #" intersect -bed -abam ",  paste0(rmdups_result_file, ".rmdups.bam"), 
                          " -b ", file.path(config$ref_data_path, config$intervals_file), 
                          " > ",paste0(on_target_result_file, "_on_target_reads"))
    runSystemCommand(myAppName, 'BEDTools', 'Intersect on_target_reads', src_command)
    
    checkIfFileExists(paste0(on_target_result_file,"_on_target_reads"))
    
    #Add to QC raport
    sample_qc_list$Number_of_on_target_reads <- getLinesNumber(paste0(on_target_result_file,"_on_target_reads"))
    sample_qc_list$Prc_of_on_target_reads <- sample_qc_list$Number_of_on_target_read/as.numeric(sample_qc_list$Number_of_reads_after_removing_duplicates)*100
    
    ###############################
    #gatk - DepthOfCoverage
    depth_of_cov_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="depth_of_coverage","temp_results_dirs"], sample_basename)
    gatk_depth_logfile <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename, "_", config_tools[config_tools$process=="depth_of_coverage","logfile"]))
    
    src_command <- paste0("java -Xms", config$java_mem," -Xmx", config$java_mem,
                          " -jar ", file.path(config$tools_path, config$gatk),
                          " -T DepthOfCoverage -R ", file.path(config$ref_data_path, config$reference_sequence_file), 
                          " -I ", paste0(clipping_result_file,".clipped.bam"),
                          " -o ", paste0(depth_of_cov_result_file,"_gatk_target_coverage"),
                          " -L ", file.path(config$ref_data_path, config$intervals_file),
                          " -ct 1 -ct 10 -ct 20", 
                          " > ", gatk_depth_logfile)
    runSystemCommand(myAppName, 'GATK', 'depth_of_coverage', src_command)
    
    checkIfFileExists(paste0(depth_of_cov_result_file,"_gatk_target_coverage"))
    checkIfFileExists(paste0(depth_of_cov_result_file,"_gatk_target_coverage.sample_summary"))
    checkIfFileExists(gatk_depth_logfile)
    
    ###############################
    #Add to QC raport
    depth_summary  <- readLines(paste0(depth_of_cov_result_file,"_gatk_target_coverage.sample_summary"))
    depth_summary  <- strsplit(depth_summary, '\t')
    sample_qc_list$Mean_coverage <- as.numeric(depth_summary[[2]][3])
    
    ###############################
    #python - methratio (BSMAP) 
    methyl_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"], sample_basename)
    methratio_logfile <- file.path(file.path(config$results_path,"logs"),paste0(sample_basename,"_",config_tools[config_tools$process=="methratio","logfile"]))
    src_command <- paste0("python ", file.path(config$anaconda_bin_path, config$methratio), 
                          " -d ", file.path(config$ref_data_path, config$reference_sequence_file), 
                          #" -s ", file.path(config$tools_path, config$samtools),
                          " -m 1 -z -i skip", 
                          " -o ", paste0(methyl_result_file,".methylation_results.txt")," ", paste0(clipping_result_file,".clipped.bam"),
                          " 2> ", methratio_logfile)
    runSystemCommand(myAppName, 'BSMAP', 'methratio', src_command)
    
    checkIfFileExists(paste0(methyl_result_file,".methylation_results.txt"))
    checkIfFileExists(methratio_logfile)
    
    ###############################
    #make bed file from methyl results
    print(paste0("File conversion..."))
    print(paste0("Reading the file: ",methyl_result_file,".methylation_results.txt"))
    methyl_res_df <- data.table::fread(paste0(methyl_result_file,".methylation_results.txt"), header = T, sep = '\t')
    methyl_res_df$end <- methyl_res_df$pos
    methyl_res_df_bed <- methyl_res_df[,c("chr","pos","end","strand","context","ratio","eff_CT_count","C_count", "CT_count")]
    print(paste0("Saving the file: ",methyl_result_file,".methylation_results.bed"))
    data.table::fwrite(methyl_res_df_bed, file=paste0(methyl_result_file,".methylation_results.bed"), quote=FALSE, sep='\t', row.names = F)
    
    ###############################
    #BEDTools - intersect capture region
    src_command <- paste0(file.path(config$anaconda_bin_path, config$bedtools), 
                          " intersect -bed -abam ",  paste0(methyl_result_file,".methylation_results.bed"), 
                          " -b ",file.path(config$ref_data_path, config$intervals_file),
                          " -u > ",paste0(methyl_result_file,".methylation_results.bed.panel"))
    runSystemCommand(myAppName, 'BEDTools', 'intersect - capture region', src_command)
    
    checkIfFileExists(paste0(methyl_result_file,".methylation_results.bed.panel"))
    
    ###############################
    # Calculate conversion efficiency
    methyl_res_panel_df <- data.table::fread(paste0(methyl_result_file,".methylation_results.bed.panel"), header = F, sep = '\t')
    
    colnames(methyl_res_panel_df) <- c("chr","start","end","strand","context","ratio","eff_CT_count","C_count", "CT_count")
    control <- methyl_res_panel_df[methyl_res_panel_df$chr=='NC_001416',]
    sample_qc_list$Number_of_Cs_in_control <- nrow(control)
    conversion_eff <- c((1-sum(control$C_count)/sum(control$CT_count))*100)
    sample_qc_list$Conversion_eff <- conversion_eff
    
    #Add to QC raport
    methyl_res_panel_df_no_control <- methyl_res_panel_df[methyl_res_panel_df$chr!='NC_001416',]
    sample_qc_list$Number_of_Cs_in_panel <- nrow(methyl_res_panel_df_no_control)
    
    methyl_res_panel_df_no_control_CpG <-methyl_res_panel_df_no_control[methyl_res_panel_df_no_control$context=='CG',]
    sample_qc_list$Number_of_Cs_in_panel_CpG <- nrow(methyl_res_panel_df_no_control_CpG)
    
    methyl_res_panel_df_no_control_nonCpG <-methyl_res_panel_df_no_control[methyl_res_panel_df_no_control$context!='CG',]
    sample_qc_list$Number_of_Cs_in_panel_non_CpG <- nrow(methyl_res_panel_df_no_control_nonCpG)
    
    methyl_res_panel_df_no_control_nonCpG_min10 <-methyl_res_panel_df_no_control_nonCpG[methyl_res_panel_df_no_control_nonCpG$eff_CT_count>=10,]
    sample_qc_list$Number_of_Cs_in_panel_non_CpG_cov_min10 <- nrow(methyl_res_panel_df_no_control_nonCpG_min10)
    
    methyl_res_panel_df_no_control_CpG_min10 <- methyl_res_panel_df_no_control_CpG[methyl_res_panel_df_no_control_CpG$eff_CT_count>=10,]
    sample_qc_list$Number_of_Cs_in_panel_CpG_cov_min10 <- nrow(methyl_res_panel_df_no_control_CpG_min10)
    
    methyl_res_panel_df_no_control_CpG_max9 <- methyl_res_panel_df_no_control_CpG[methyl_res_panel_df_no_control_CpG$eff_CT_count<10,]
    sample_qc_list$Number_of_Cs_in_panel_CpG_cov_max9 <- nrow(methyl_res_panel_df_no_control_CpG_max9)
    
    sample_qc_list$Prc_of_Cs_in_panel_CpG_cov_min10 <-  nrow(methyl_res_panel_df_no_control_CpG_min10)/sum(nrow(methyl_res_panel_df_no_control_CpG_max9)+nrow(methyl_res_panel_df_no_control_CpG_min10))*100
    sample_qc_list$Prc_of_Cs_in_panel_CpG_cov_max9 <- nrow(methyl_res_panel_df_no_control_CpG_max9)/sum(nrow(methyl_res_panel_df_no_control_CpG_max9)+nrow(methyl_res_panel_df_no_control_CpG_min10))*100
    
    
    #Write yaml raport
    yaml::write_yaml(sample_qc_list, file.path(config$results_path,"QC_report",paste0(sample_basename,"_QC_summary")))
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
    print('#############################################')
    print(paste0("Processing the sample - '", sample_basename, "' is finished. [", format(full_stop_time - full_start_time, digits=3) ,"]"))
    print('#############################################')
  }
  return(T)
}
