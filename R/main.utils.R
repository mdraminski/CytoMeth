#########################
######  run_seqtk  ######
#########################
run_seqtk <- function(config, config_tools){
  seqtk_r1_file <- file.path(config$results_path, config_tools[config_tools$proces=="seqtk","temp_results_dirs"], basename(config$file_r1))
  seqtk_r2_file <- file.path(config$results_path, config_tools[config_tools$proces=="seqtk","temp_results_dirs"], basename(config$file_r2))
  
  if(config$overwrite_results | !(file.exists(seqtk_r1_file) & file.exists(seqtk_r2_file))){
    src_command <- paste0(file.path(config$anaconda_bin_path, config$seqtk), " sample -s 10000 ", config$file_r1, " ", config$sqtk_subset, " > ", seqtk_r1_file)
    runSystemCommand(config$myAppName, 'seqtk', 'subsample of reads R1', src_command, config$verbose)
    
    src_command <- paste0(file.path(config$anaconda_bin_path, config$seqtk), " sample -s 10000 ", config$file_r2, " ", config$sqtk_subset, " > ", seqtk_r2_file)
    runSystemCommand(config$myAppName, 'seqtk', 'subsample of reads R2', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'seqtk', 'subsample of reads', dirname(seqtk_r1_file))
  }
  
  if(!checkIfFileExists(seqtk_r1_file)) return(NULL)
  if(!checkIfFileExists(seqtk_r2_file)) return(NULL)
  
  config$file_r1 <- seqtk_r1_file
  config$file_r2 <- seqtk_r2_file
  
  return(config)
}

#########################
######  run_FastQC  ######
#########################
run_FastQC <- function(config, config_tools){
  fastqc_result_dir <- file.path(config$results_path, config_tools[config_tools$tool=="fastqc","temp_results_dirs"])
  fastqc_result_r1_file <- file.path(fastqc_result_dir, paste0(basename_noext(config$file_r1),"_fastqc.html"))
  fastqc_result_r2_file <- file.path(fastqc_result_dir, paste0(basename_noext(config$file_r2),"_fastqc.html"))

  if(config$overwrite_results | !(file.exists(fastqc_result_r1_file) & file.exists(fastqc_result_r2_file))){
    # src_command <- paste0(file.path(config$anaconda_bin_path, config$fastqc), " --nogroup ", 
    #                       seqtk_result_file,"_subset_R1.fastq ",
    #                       seqtk_result_file,"_subset_R2.fastq", 
    #                       " --outdir ", fastqc_result_dir)
    
    src_command <- paste0(file.path(config$anaconda_bin_path, config$fastqc), " -t ", conf$threads, " --nogroup ", 
                          config$file_r1, " ", config$file_r2, " --outdir ", fastqc_result_dir)
    runSystemCommand(config$myAppName, 'fastqc', 'FastQC Report generation', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'FastQC', 'examine sequence read quality', fastqc_result_dir)
  }
  if(!checkIfFileExists(fastqc_result_r1_file)) return(NULL)
  if(!checkIfFileExists(fastqc_result_r2_file)) return(NULL)
  
  return(config)
}

###############################
######  run_Trimmomatic  ######
###############################
run_Trimmomatic <- function(config, config_tools){
  trimming_result_r1_file <- file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"], basename_noext(config$file_r1))
  trimming_result_r2_file <- file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"], basename_noext(config$file_r2))
  trimmomatic_out_logfile <- file.path(config$results_path,"logs",paste0(basename_sample(config$file_r1),"_",config_tools[config_tools$tool=="trimmomatic","logfile"]))
  
  if(config$overwrite_results | !(file.exists(paste0(trimming_result_r1_file,"_trimmed.fq")) & file.exists(paste0(trimming_result_r2_file,"_trimmed.fq")) &
                                  file.exists(paste0(trimming_result_r1_file,"_unpaired.fq")) & file.exists(paste0(trimming_result_r2_file,"_unpaired.fq")) &
                                  file.exists(trimmomatic_out_logfile)
  )){
    src_command <- paste0("java -Xms", f2si(config$memory)," -Xmx", f2si(config$memory)," -jar ", file.path(config$tools_path,config$trimmomatic),
                          " PE -threads ",config$threads," -phred33 ",config$file_r1," ", config$file_r2," ",
                          trimming_result_r1_file,"_trimmed.fq ", trimming_result_r1_file,"_unpaired.fq ",
                          trimming_result_r2_file,"_trimmed.fq ", trimming_result_r2_file,"_unpaired.fq ",
                          "ILLUMINACLIP:", file.path(config$tools_path, config$ref_data_trimmomatic_adapter),
                          ":2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:",config$trimmomatic_MINLEN,
                          " 2> ",trimmomatic_out_logfile)
    runSystemCommand(config$myAppName, 'Trimmomatic', 'trimming', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Trimmomatic', 'trimming', dirname(trimming_result_r1_file))
  }
  if(!checkIfFileExists(paste0(trimming_result_r1_file,"_trimmed.fq"))) return(NULL)
  if(!checkIfFileExists(paste0(trimming_result_r2_file,"_trimmed.fq"))) return(NULL)
  if(!checkIfFileExists(paste0(trimming_result_r1_file,"_unpaired.fq"))) return(NULL)
  if(!checkIfFileExists(paste0(trimming_result_r2_file,"_unpaired.fq"))) return(NULL)
  if(!checkIfFileExists(trimmomatic_out_logfile)) return(NULL)
  if(!checkLog(trimmomatic_out_logfile, "Completed successfully", 'Trimming')) return(NULL)
  
  config$tmp_files <- c(config$tmp_files, paste0(trimming_result_r1_file, c("_trimmed.fq","_unpaired.fq")), paste0(trimming_result_r2_file, c("_trimmed.fq","_unpaired.fq")))
  return(config)
}  

#########################
######  run_BSMAP  ######
#########################
run_BSMAP <- function(config, config_tools){
  trimming_result_r1_file <- file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"], basename_noext(config$file_r1))
  trimming_result_r2_file <- file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"], basename_noext(config$file_r2))
  mapping_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"], basename_sample(config$file_r1))
  mapping_dir <- dirname(mapping_result_file)
  
  if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".sam")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path, config$bsmap), " -r 0 -s 16 -n 1 ",
                          getBSMAPIndexInterval(config$memory),
                          " -a ", trimming_result_r1_file,"_trimmed.fq",
                          " -b ", trimming_result_r2_file,"_trimmed.fq",
                          " -d ", file.path(config$ref_data_path, config$ref_data_sequence_file),
                          " -p ", min(config$threads, 8)," -o ", paste0(mapping_result_file, ".sam"))
    runSystemCommand(config$myAppName, 'BSMAP', 'mapping', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'BSMAP', 'mapping', mapping_dir)
  }
  if(!checkIfFileExists(paste0(mapping_result_file, ".sam"))) return(NULL)
  
  # Picard - SAM -> BAM
  if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".bam")))) {
    src_command <- paste0("java -Xms", f2si(config$memory)," -Xmx", f2si(config$memory), 
                          config$picard_parser,
                          " -jar ", file.path(config$tools_path,config$picard), 
                          " AddOrReplaceReadGroups",
                          " VALIDATION_STRINGENCY=LENIENT ",
                          " INPUT=", paste0(mapping_result_file, ".sam"),
                          " OUTPUT=", paste0(mapping_result_file, ".bam"),
                          " CREATE_INDEX=true",
                          " RGID=SAMPLE RGLB=SAMPLE RGPL=illumina RGSM=SAMPLE RGPU=platform_unit")
    runSystemCommand(config$myAppName, 'Picard', 'sam to bam conversion', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Picard', 'sam to bam conversion', mapping_dir)
  }
  if(!checkIfFileExists(paste0(mapping_result_file, ".bam"))) return(NULL)
  
  config$file_bam <- paste0(mapping_result_file, '.bam')
  config$tmp_files <- c(config$tmp_files, paste0(mapping_result_file, c(".sam",".bam")))
  return(config)
}

####################################
######  run_RemoveDuplicates  ######
####################################
run_RemoveDuplicates <- function(config, config_tools){
  sample_basename <- basename_sample(config$file_bam)
  mapping_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"], basename_sample(config$file_bam))
  mapping_dir <- dirname(mapping_result_file)
  
  ######  bamtools - split
  if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".TAG_ZS_+-.bam")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path, config$bamtools), " split -tag ZS",
                          " -in ", config$file_bam)
    runSystemCommand(config$myAppName, 'Bamtools', 'split', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Bamtools', 'split', mapping_dir)
  }
  if(!checkIfFileExists(paste0(mapping_result_file, ".TAG_ZS_++.bam"))) return(NULL)
  if(!checkIfFileExists(paste0(mapping_result_file, ".TAG_ZS_+-.bam"))) return(NULL)
  if(!checkIfFileExists(paste0(mapping_result_file, ".TAG_ZS_--.bam"))) return(NULL)
  if(!checkIfFileExists(paste0(mapping_result_file, ".TAG_ZS_-+.bam"))) return(NULL)
  
  ######  bamtools - merge top
  if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".top.bam")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " merge",
                          " -in ", paste0(mapping_result_file,".TAG_ZS_++.bam"),
                          " -in ", paste0(mapping_result_file,".TAG_ZS_+-.bam"),
                          " -out ", paste0(mapping_result_file,".top.bam"))
    runSystemCommand(config$myAppName, 'Bamtools', 'merge top', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Bamtools', 'merge top', mapping_dir)
  }
  if(!checkIfFileExists(paste0(mapping_result_file, ".top.bam"))) return(NULL)
  
  ######  bamtools - merge bottom
  if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".bottom.bam")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path, config$bamtools), " merge",
                          " -in ", paste0(mapping_result_file,".TAG_ZS_-+.bam"),
                          " -in ", paste0(mapping_result_file,".TAG_ZS_--.bam"),
                          " -out ", paste0(mapping_result_file,".bottom.bam"))
    runSystemCommand(config$myAppName, 'Bamtools', 'merge bottom', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Bamtools', 'merge bottom', mapping_dir)
  }
  if(!checkIfFileExists(paste0(mapping_result_file, ".bottom.bam"))) return(NULL)
  
  ######  bamtools - sort top
  if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".top.bam.sorted")))) {
    #bamtools deprecated    
    #src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " sort ",
    #                      " -in ", paste0(mapping_result_file,".top.bam"), 
    #                      " -out ", paste0(mapping_result_file,".top.bam.sorted"))
    #samtools multithreaded
    src_command <- paste0(file.path(config$anaconda_bin_path,config$samtools), " sort",
                          " -@ ", config$threads, " -m 4G",
                          " -o ", paste0(mapping_result_file,".top.bam.sorted"),
                          " ", paste0(mapping_result_file,".top.bam"))
    runSystemCommand(config$myAppName, 'samtools', 'sort top', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'samtools', 'sort top', mapping_dir)
  }
  if(!checkIfFileExists(paste0(mapping_result_file, ".top.bam.sorted"))) return(NULL)
  
  ######  bamtools - sort bottom
  if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".bottom.bam.sorted")))) {
    #bamtools deprecated    
    #src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " sort ",
    #                      " -in ", paste0(mapping_result_file,".bottom.bam"), 
    #                      " -out ", paste0(mapping_result_file,".bottom.bam.sorted"))
    #samtools multithreaded
    src_command <- paste0(file.path(config$anaconda_bin_path,config$samtools), " sort",
                          " -@ ", config$threads, " -m 4G",
                          " -o ", paste0(mapping_result_file,".bottom.bam.sorted"),
                          " ", paste0(mapping_result_file,".bottom.bam"))
    runSystemCommand(config$myAppName, 'samtools', 'sort bottom', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'samtools', 'sort bottom', mapping_dir)
  }
  if(!checkIfFileExists(paste0(mapping_result_file, ".bottom.bam.sorted"))) return(NULL)
  
  ######  picard - MarkDuplicates - top
  ###### TODO: samtools markdup instead?
  rmdups_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_top","temp_results_dirs"], sample_basename)
  rmdups_logfile <- file.path(config$results_path, "logs", paste0(sample_basename, "_", config_tools[config_tools$process=="mark_duplicates_top","logfile"]))
  rmdups_dir <- dirname(rmdups_result_file)
  if(config$overwrite_results | !(file.exists(paste0(rmdups_result_file, ".top.rmdups.bam")) 
                                  & file.exists(paste0(rmdups_result_file, ".top.rmdups_metrics.txt")) 
                                  & file.exists(rmdups_logfile))) {
    src_command <- paste0("java -Xms", f2si(config$memory)," -Xmx", f2si(config$memory),
                          config$picard_parser,
                          " -XX:+UseG1GC -XX:MaxGCPauseMillis=100",
                          " -jar ",file.path(config$tools_path, config$picard), 
                          " MarkDuplicates ",
                          " VALIDATION_STRINGENCY=LENIENT",
                          " INPUT=", paste0(mapping_result_file,".top.bam.sorted"),
                          " OUTPUT=", paste0(rmdups_result_file,".top.rmdups.bam"),
                          " METRICS_FILE=",paste0(rmdups_result_file,".top.rmdups_metrics.txt"),
                          " REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true",
                          " 2> ", rmdups_logfile)
    runSystemCommand(config$myAppName, 'Picard', 'MarkDuplicates - top', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Picard', 'MarkDuplicates - top', rmdups_dir)
  }
  if(!checkIfFileExists(paste0(rmdups_result_file, ".top.rmdups.bam"))) return(NULL)
  if(!checkIfFileExists(paste0(rmdups_result_file, ".top.rmdups_metrics.txt"))) return(NULL)
  if(!checkIfFileExists(rmdups_logfile)) return(NULL)
  if(!checkLog(rmdups_logfile, "MarkDuplicates done. Elapsed time:", 'MarkDuplicates - top')) return(NULL)
  
  ######  picard - MarkDuplicates - bottom
  rmdups_logfile <- file.path(config$results_path, "logs", paste0(sample_basename, "_", config_tools[config_tools$process=="mark_duplicates_bottom","logfile"]))
  if(config$overwrite_results | !(file.exists(paste0(rmdups_result_file, ".bottom.rmdups.bam")) &
                                  file.exists(paste0(rmdups_result_file, ".bottom.rmdups_metrics.txt")) &
                                  file.exists(rmdups_logfile)
  )) {
    src_command <- paste0("java -Xms", f2si(config$memory)," -Xmx", f2si(config$memory),
                          " -XX:+UseG1GC -XX:MaxGCPauseMillis=100",
                          config$picard_parser,
                          " -jar ",file.path(config$tools_path,config$picard), 
                          " MarkDuplicates ",
                          " VALIDATION_STRINGENCY=LENIENT",
                          " INPUT=", paste0(mapping_result_file,".bottom.bam.sorted"),
                          " OUTPUT=", paste0(rmdups_result_file,".bottom.rmdups.bam"),
                          " METRICS_FILE=",paste0(rmdups_result_file,".bottom.rmdups_metrics.txt"),
                          " REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true",
                          " 2> ", rmdups_logfile)
    runSystemCommand(config$myAppName, 'Picard', 'MarkDuplicates - bottom', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Picard', 'MarkDuplicates - bottom', rmdups_dir)
  }
  if(!checkIfFileExists(paste0(rmdups_result_file, ".bottom.rmdups.bam"))) return(NULL)
  if(!checkIfFileExists(paste0(rmdups_result_file, ".bottom.rmdups_metrics.txt"))) return(NULL)
  if(!checkIfFileExists(rmdups_logfile)) return(NULL)
  if(!checkLog(rmdups_logfile, "MarkDuplicates done. Elapsed time:", 'MarkDuplicates - bottom')) return(NULL)
  
  ######  bamtools - merge
  if(config$overwrite_results | !(file.exists(paste0(rmdups_result_file,".rmdups.bam")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " merge", 
                          " -in ", paste0(rmdups_result_file,".top.rmdups.bam"),
                          " -in ", paste0(rmdups_result_file,".bottom.rmdups.bam"),
                          " -out ", paste0(rmdups_result_file,".rmdups.bam"))
    runSystemCommand(config$myAppName, 'Bamtools', 'merge', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Bamtools', 'merge', rmdups_dir)
  }
  if(!checkIfFileExists(paste0(rmdups_result_file,".rmdups.bam"))) return(NULL)
  
  config$tmp_files <- c(config$tmp_files, 
                        paste0(mapping_result_file, c(".bam",".top.bam",".bottom.bam",".TAG_ZS_++.bam",".TAG_ZS_+-.bam",".TAG_ZS_-+.bam",".TAG_ZS_--.bam",".top.bam.sorted",".bottom.bam.sorted")),
                        paste0(rmdups_result_file, c(".top.rmdups.bam",".bottom.rmdups.bam",".top.rmdups.bai",".bottom.rmdups.bai",".rmdups.bam")))
  return(config)
}

#############################
######  run_FilterBAM  ######
#############################
run_FilterBAM <- function(config, config_tools){
  sample_basename <- basename_sample(config$file_bam)
  rmdups_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_top","temp_results_dirs"], sample_basename)
  
  ######  bamtools - filter
  filtered_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="filter_bam","temp_results_dirs"], sample_basename)
  if(config$overwrite_results | !(file.exists(paste0(filtered_result_file,".filtered.bam")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " filter -isMapped true -isPaired true -isProperPair true -forceCompression",
                          " -in ", paste0(rmdups_result_file,".rmdups.bam"), 
                          " -out ", paste0(filtered_result_file,".filtered.bam"))
    runSystemCommand(config$myAppName, 'Bamtools', 'filter', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Bamtools', 'filter', dirname(filtered_result_file))
  }
  if(!checkIfFileExists(paste0(filtered_result_file,".filtered.bam"))) return(NULL)
  
  ######  Samtools - flagstat flagstat_filtered_bam
  flagstat_result_file <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename,"_",config_tools[config_tools$process=="flagstat_filtered_bam","logfile"]))
  if(config$overwrite_results | !(file.exists(paste0(flagstat_result_file)))) {
    src_command <- paste0(file.path(config$anaconda_bin_path, config$samtools), " flagstat ", 
                          " -@ ", config$threads," ",
                          paste0(filtered_result_file,".filtered.bam"), 
                          " 1> ", flagstat_result_file)
    runSystemCommand(config$myAppName, 'Samtools', 'flagstat', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Samtools', 'flagstat', dirname(flagstat_result_file))
  }
  if(!checkIfFileExists(flagstat_result_file)) return(NULL)
  
  config$tmp_files <- c(config$tmp_files, paste0(filtered_result_file,".filtered.bam"))
  return(config)
}

###############################
######  run_ClipOverlap  ######
###############################
run_ClipOverlap <- function(config, config_tools){
  sample_basename <- basename_sample(config$file_bam)
  filtered_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="filter_bam","temp_results_dirs"], sample_basename)
  
  ######  bamUtil - clipOverlap
  clipping_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"], sample_basename)
  clipOverlap_logfile <- file.path(file.path(config$results_path,"logs"),paste0(sample_basename,"_",config_tools[config_tools$process=="clip_overlap","logfile"]))
  if(config$overwrite_results | !(file.exists(paste0(clipping_result_file,".clipped.bam")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path, config$bamUtil), " clipOverlap --stats",
                          " --in ",  paste0(filtered_result_file,".filtered.bam"),
                          " --out ", paste0(clipping_result_file,".clipped.bam"), 
                          " 2> ", clipOverlap_logfile)
    runSystemCommand(config$myAppName, 'BamUtil', 'clipOverlap', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'BamUtil', 'clipOverlap', dirname(clipping_result_file))
  }
  
  if(!checkIfFileExists(paste0(clipping_result_file,".clipped.bam"))) return(NULL)
  if(!checkLog(clipOverlap_logfile, "Completed ClipOverlap Successfully.", 'BamUtil - clipOverlap')) return(NULL)
  
  ######  samtools - index
  if(config$overwrite_results | !(file.exists(paste0(clipping_result_file,".clipped.bam.bai")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path, config$samtools), 
                          " index ", " -@ ", config$threads," ", 
                          paste0(clipping_result_file,".clipped.bam"))
    runSystemCommand(config$myAppName, 'Samtools', 'index', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Samtools', 'index', dirname(clipping_result_file))
  }
  
  if(!checkIfFileExists(paste0(clipping_result_file,".clipped.bam.bai"))) return(NULL)
  
  if(config$remove_clipped_bam){
    config$tmp_files <- c(config$tmp_files, paste0(clipping_result_file,".clipped.bam"), paste0(clipping_result_file,".clipped.bam.bai"))
  }
  
  return(config)
}

##################################
######  run_MappingMetrics  ######
##################################
run_MappingMetrics <- function(config, config_tools){
  sample_basename <- basename_sample(config$file_bam)
  filtered_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="filter_bam","temp_results_dirs"], sample_basename)
  ### TODO: ?? small_FAKE03.rmdups.bam not small_FAKE03.clipped.bam
  rmdups_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_top","temp_results_dirs"], sample_basename)
  ######  picard - Basic Mapping Metrics
  basic_mapping_metrics_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="basic_mapping_metrics","temp_results_dirs"], sample_basename)
  if(config$overwrite_results | !(file.exists(paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt")))) {
    src_command <- paste0("java -Xms", f2si(config$memory)," -Xmx", f2si(config$memory), 
                          config$picard_parser,
                          " -jar ", file.path(config$tools_path,config$picard), 
                          " CollectAlignmentSummaryMetrics ",
                          " VALIDATION_STRINGENCY=LENIENT",
                          " METRIC_ACCUMULATION_LEVEL=ALL_READS ",
                          " INPUT=", paste0(rmdups_result_file,".rmdups.bam"), ##TODO: ?? clipped.bam
                          " OUTPUT=", paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt"),
                          " REFERENCE_SEQUENCE=", file.path(config$ref_data_path, config$ref_data_sequence_file))
    runSystemCommand(config$myAppName, 'Picard', 'Basic Mapping Metrics', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Picard', 'Basic Mapping Metrics', dirname(basic_mapping_metrics_result_file))
  }
  if(!checkIfFileExists(paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt"))) return(NULL)
  
  ######  picard - Insert Size Metrics
  insert_size_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="insert_size","temp_results_dirs"], sample_basename)
  if(config$overwrite_results | !(file.exists(paste0(insert_size_result_file,"_insert_size_metrics.txt")))) {
    src_command <- paste0("java -Xms", f2si(config$memory)," -Xmx", f2si(config$memory), 
                          config$picard_parser,
                          " -jar ", file.path(config$tools_path, config$picard), 
                          " CollectInsertSizeMetrics ",
                          " VALIDATION_STRINGENCY=LENIENT ",
                          " Histogram_FILE=", paste0(insert_size_result_file,"_insert_size_plot.pdf"),
                          " INPUT=", paste0(filtered_result_file,".filtered.bam"),
                          " OUTPUT=", paste0(insert_size_result_file,"_insert_size_metrics.txt"))
    runSystemCommand(config$myAppName, 'Picard', 'Insert Size', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Picard', 'Insert Size', dirname(insert_size_result_file))
  }
  if(!checkIfFileExists(paste0(insert_size_result_file,"_insert_size_metrics.txt"))) return(NULL)
  
  
  ######  On-target reads (BEDTools - intersect)
  on_target_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="on_target_reads","temp_results_dirs"], sample_basename)
  if(config$overwrite_results | !(file.exists(paste0(on_target_result_file,"_on_target_reads")))) {
    ## If the panel file is defined use the standard path - intersect
    if(str_trim(config$ref_data_intervals_file) != ''){
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bedtools), 
                            " intersect -bed -abam ",  paste0(rmdups_result_file, ".rmdups.bam"), 
                            " -b ", file.path(config$ref_data_path, config$ref_data_intervals_file), 
                            " > ",paste0(on_target_result_file, "_on_target_reads"))
      ### TODO: ??? small_FAKE03.rmdups.bam? Shouldn't be .clipped.bam?
      runSystemCommand(config$myAppName, 'BEDTools', 'Intersect on_target_reads', src_command, config$verbose)  
    }else{
      ## If the panel file is not defined use convert bam2bed - (bedtools bamtobed -i reads.bam)
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bedtools), 
                            " bamtobed -i ",  paste0(rmdups_result_file, ".rmdups.bam"),
                            " > ",paste0(on_target_result_file, "_on_target_reads"))
      runSystemCommand(config$myAppName, 'BEDTools', 'bamtobed on_target_reads', src_command, config$verbose)
    }
  }else{
    skipProcess(config$myAppName, 'BEDTools', 'Intersect on_target_reads', dirname(on_target_result_file))
  }
  if(!checkIfFileExists(paste0(on_target_result_file,"_on_target_reads"))) return(NULL)
  
  #update and save result SAMPLE_on_target_reads.txt file
  on_target_result_file_txt <- paste0(on_target_result_file,"_on_target_reads.txt")
  if(config$overwrite_results | !(file.exists(on_target_result_file_txt))){
    number_of_on_target_reads <- getLinesNumber(paste0(on_target_result_file,"_on_target_reads"))
    on_target_params <- list(number_of_on_target_reads = as.character(number_of_on_target_reads))
    yaml::write_yaml(on_target_params, on_target_result_file_txt)
  }
  
  config$tmp_files <- c(config$tmp_files, paste0(on_target_result_file,"_on_target_reads"))
  return(config)
}

##################################
######  run_DepthOfCoverage  #####
##################################
run_DepthOfCoverage <- function(config, config_tools){
  sample_basename <- basename_sample(config$file_bam)
  clipping_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"], sample_basename)
  
  ######  gatk - DepthOfCoverage
  depth_of_cov_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="depth_of_coverage","temp_results_dirs"], sample_basename)
  gatk_depth_logfile <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename, "_", config_tools[config_tools$process=="depth_of_coverage","logfile"]))
  if(config$overwrite_results | !(file.exists(paste0(depth_of_cov_result_file,"_gatk_target_coverage")) &
                                  file.exists(paste0(depth_of_cov_result_file,"_gatk_target_coverage.sample_summary")) & file.exists(gatk_depth_logfile))){
    panel_command <- ""
    if(str_trim(config$ref_data_intervals_file != "")){
      panel_command <- paste0(" -L ", file.path(config$ref_data_path, config$ref_data_intervals_file))
    }
    src_command <- paste0("java -Xms", f2si(config$memory)," -Xmx", f2si(config$memory),
                          " -jar ", file.path(config$tools_path, config$gatk),
                          " -T DepthOfCoverage -R ", file.path(config$ref_data_path, config$ref_data_sequence_file), 
                          " -I ", paste0(clipping_result_file,".clipped.bam"),
                          " -o ", paste0(depth_of_cov_result_file,"_gatk_target_coverage"),
                          panel_command,
                          " -ct 1 -ct 10 -ct 20", 
                          " > ", gatk_depth_logfile)
    runSystemCommand(config$myAppName, 'GATK', 'depth_of_coverage', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'GATK', 'depth_of_coverage', dirname(depth_of_cov_result_file))
  }
  if(!checkIfFileExists(paste0(depth_of_cov_result_file,"_gatk_target_coverage"))) return(NULL)
  if(!checkIfFileExists(paste0(depth_of_cov_result_file,"_gatk_target_coverage.sample_summary"))) return(NULL)
  if(!checkIfFileExists(gatk_depth_logfile)) return(NULL)
  
  config$tmp_files <- c(config$tmp_files, paste0(depth_of_cov_result_file, c("_gatk_target_coverage","_gatk_target_coverage.sample_cumulative_coverage_counts",
                                                                             "_gatk_target_coverage.sample_interval_summary","_gatk_target_coverage.sample_cumulative_coverage_proportions",
                                                                             "_gatk_target_coverage.sample_interval_statistics")))
  return(config)
}

##################################
######     run_Methratio     #####
##################################
run_Methratio <- function(config, config_tools){
  start_time_meth <- Sys.time()
  chr_output_files <- NULL
  sample_basename <- basename_sample(config$file_bam)
  clipping_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"], sample_basename)
  clipping_result_file_bam <- paste0(clipping_result_file,".clipped.bam")
  clipping_result_file_bed <- paste0(clipping_result_file,".clipped.bed")
  clipping_result_file_sam <- stringr::str_replace(paste0(clipping_result_file,".clipped.sam"),"clipped","tmp")
  
  ######  python - methratio (BSMAP)
  methyl_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"], sample_basename)
  methyl_result_file_txt <- paste0(methyl_result_file,".methylation_results.txt")
  methyl_result_file_bed <- paste0(methyl_result_file,".methylation_results.rough.bed")
  methratio_logfile <- file.path(file.path(config$results_path,"logs"),paste0(sample_basename,"_",config_tools[config_tools$process=="methratio","logfile"]))
  dir.create(dirname(methyl_result_file_txt), showWarnings = FALSE)
  
  if(config$overwrite_results | !(file.exists(methyl_result_file_txt) & file.exists(methratio_logfile))) {
    if(tolower(config$meth_processing) == tolower("allCHR")){
      ######  Methratio ALL CHR PROCESSING
      src_command <- paste0(config$python2, " ", file.path(config$tools_path, config$methratio), 
                            " -d ", file.path(config$ref_data_path, config$ref_data_sequence_file), 
                            #" -s ", config$anaconda_bin_path,
                            " -m ", config$min_depth," -z -i skip", 
                            " -o ", methyl_result_file_txt," ", clipping_result_file_bam,
                            " 2> ", methratio_logfile)
      runSystemCommand(myAppName, 'BSMAP', 'methratio', src_command, config$verbose)
    }else if(startsWith(tolower(config$meth_processing),tolower("batchCHR"))){
      ######  Methratio CHR BATCH PROCESSING - SPLIT INPUT FILE (HDD)
      print(paste0("### Methratio CHR BATCH PROCESSING ###"))
      runSystemCommand(config$myAppName, 'rm', 'tmp folder', paste0("rm -rf ", config$tmp_data_path, "/*"), FALSE)
      print(paste0("Conversion of input clipped.bam to clipped.sam [bamToSam]..."))
      src_command <- paste0(file.path(config$anaconda_bin_path, config$samtools), " view ", " -@ ", config$threads," ", clipping_result_file_bam," > ",clipping_result_file_sam)
      runSystemCommand(config$myAppName, 'samtools', 'bamView', src_command, config$verbose)
      print(paste0("Splitting clipped sam file..."))
      src_command <- paste0("cd ", dirname(clipping_result_file_sam), "; awk '{print>$3}' ", basename(clipping_result_file_sam), ";")
      runSystemCommand(config$myAppName, 'awk', 'splitclippedbed', src_command, config$verbose)
      file.remove(clipping_result_file_sam)
      inputFiles <- list.files(dirname(clipping_result_file_sam),full.names = T)
      inputFiles <- sapply(inputFiles, function(x) {file.rename(x,paste0(x,".bsam"))})
      inputFiles <- list.files(dirname(clipping_result_file_sam),full.names = T)
      sample_chr <- stringr::str_replace(basename(inputFiles),".bsam","")
      #common processing for chr batch processing
      chr_output_files <- list()
      logAppend <- ">"
      for(i in 1:length(sample_chr)){
        print(paste0("#####################################################"))
        print(paste0("##### Running methratio for chr: '",sample_chr[i],"' [",i,"/",length(sample_chr),"] #####"))
        current_chr <- sample_chr[i]
        current_input <- inputFiles[i]
        current_output <- stringr::str_replace(methyl_result_file_txt, ".methylation_results.txt", paste0(".methylation_results.",sample_chr[i],".txt"))
        chr_output_files[[i]] <- current_output
        src_command <- paste0(config$python2, " ", file.path(config$tools_path, config$methratio), 
                              " -c ", current_chr," -d ", file.path(config$ref_data_path, config$ref_data_sequence_file), 
                              #" -s ", config$anaconda_bin_path,
                              " -m ", config$min_depth, " -z -i skip", 
                              " -o ", current_output, " ", current_input,
                              " 2", logAppend, methratio_logfile)
        #" 1", logAppend, methratio_logfile, " 2>&1")
        runSystemCommand(config$myAppName, 'BSMAP', 'methratio', src_command, config$verbose)
        logAppend <- ">>"
        if(file.exists(current_output))
          file.remove(current_input)
      }
      chr_output_files <- unlist(chr_output_files)
      if(!all(sapply(c(chr_output_files), checkIfFileExists))){
        cat(paste0(readLines(methratio_logfile),"\n"))
        return(NULL)
      }
      methyl_result_data <- list()
      ######  Make one txt file out of separated chr methyl results
      print(paste0("Methylation files binding ..."))
      for(i in 1:length(chr_output_files)){
        print(paste0("Reading the file: ", chr_output_files[i]))
        methyl_result_data[[i]] <- data.table::fread(chr_output_files[i], header = T, sep = '\t')
      }
      methyl_result_data <- data.table::rbindlist(methyl_result_data)
      data.table::fwrite(methyl_result_data, methyl_result_file_txt, sep = '\t')
      print("Removing temporary Methylation (CHR) files...")
      files_removed <- sapply(unique(chr_output_files), fileRemoveIfExists)
    }
  }else{
    skipProcess(config$myAppName, 'BSMAP', 'methratio', dirname(methyl_result_file_txt))
  }
  
  if(!checkIfFileExists(methratio_logfile)) return(NULL)
  if(!checkIfFileExists(methyl_result_file_txt)) {
    cat(paste0(readLines(methratio_logfile),"\n"))
    return(NULL)
  }
  
  stop_time_meth <- Sys.time()
  cat(paste0(config$myAppName, ": Methratio execution is finished. [", format(stop_time_meth - start_time_meth, digits=3) ,"]\n"))
  
  config$tmp_files <- c(config$tmp_files, methyl_result_file_txt, clipping_result_file_bed, clipping_result_file_sam, chr_output_files)
  return(config)
}

############################################
######  bssnper calculate methylation  #####
############################################
methratio_to_bed <- function(config, config_tools){
  start_time_meth <- Sys.time()
  
  sample_basename <- basename_sample(config$file_bam)
  methyl_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"], sample_basename)
  methyl_result_file_txt <- paste0(methyl_result_file,".methylation_results.txt")
  methyl_result_file_bed <- paste0(methyl_result_file,".methylation_results.rough.bed")

  if(config$overwrite_results | !(file.exists(methyl_result_file_bed))) {
    ######  Make bed file from methyl results
    print(paste0("File conversion txt -> bed. Reading the file: ", methyl_result_file_txt))
    methyl_result_data <- data.table::fread(methyl_result_file_txt, header = T, sep = '\t')
    methyl_result_data$start <- methyl_result_data$pos-1  ## setting the bed's start
    methyl_result_data$end <- methyl_result_data$pos
    methyl_result_data <- methyl_result_data[,c("chr","start","end","context","ratio","strand", "eff_CT_count", "C_count", "CT_count", "pos")]
    names(methyl_result_data) <- getMethylDataHeader(version = 2, size = 10)
    
    #head(methyl_result_data)
    #summary(methyl_result_data)
    if(nrow(methyl_result_data) == 0){
      print(paste0("Error! methyl_result_data is empty!"))
      print(paste0("QC Report of the Sample: ", sample_basename))
      sampleQC <- getSampleQCSummary(sample_basename, config)
      print(qc2dataframe(sampleQC))
      return(NULL)
    }else{
      methyl_result_info(methyl_result_data)
    }
  
    ## for testing 
    ## ref_fasta <- file.path(config$ref_data_path, config$ref_data_sequence_file)
    ## checkLetters(methyl_result_data,ref_fasta) 
    print(paste0("Saving the file: ", methyl_result_file_bed))
    data.table::fwrite(methyl_result_data, file = methyl_result_file_bed, quote=FALSE, sep='\t', row.names = F, col.names = F)
  } else{
    skipProcess(config$myAppName, 'methratio', 'txt2bed', dirname(methyl_result_file_bed))
  }
  
  if(!checkIfFileExists(methyl_result_file_bed)) {
    return(NULL)
  }
  
  stop_time_meth <- Sys.time()
  cat(paste0(config$myAppName, ": Conversion txt2bed is finished. [", format(stop_time_meth - start_time_meth, digits=3) ,"]\n"))
  
  config$tmp_files <- c(config$tmp_files)
  return(config)
}

##################################
######  run_PanelIntersect  #####
##################################
run_PanelIntersect <- function(config, config_tools){
  start_time_meth <- Sys.time()
  sample_basename <- basename_sample(config$file_bam)
  ######  python - methratio (BSMAP)
  methyl_result_file <- file.path(config$results_path, config_tools[config_tools$proces==config$meth_tool_current,"temp_results_dirs"], sample_basename)
  methyl_result_file_rough <- paste0(methyl_result_file,".methylation_results.rough.bed")
  methyl_result_file_bed <- paste0(methyl_result_file,".methylation_results.bed")
  ######  VCF SNP files
  vcf_in                <- paste0(methyl_result_file, '.sorted.vcf')
  vcf_out               <- paste0(methyl_result_file, '.snp.vcf')
  
  ######  BEDTools - intersect capture region (only if file is defined)
  if(str_trim(config$ref_data_intervals_file) != ''){
    print(paste0("Intersection of meth data with panel intervals..."))
    if(config$overwrite_results | !(file.exists(methyl_result_file_bed))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bedtools), 
                            " intersect -bed -abam ",  methyl_result_file_rough, 
                            " -b ", file.path(config$ref_data_path, config$ref_data_intervals_file),
                            " -u > ", methyl_result_file_bed)
      runSystemCommand(config$myAppName, 'BEDTools', 'intersect - meth data capture region', src_command, config$verbose)
      
      if(file.exists(vcf_in)){
        src_command <- paste0("bedtools intersect -header -a ", vcf_in, " -b ",file.path(config$ref_data_path, config$ref_data_intervals_file), " > ", vcf_out)
        runSystemCommand(config$myAppName, 'BEDTools', 'intersect - snp data capture region', src_command, config$verbose)
      }
    }else{
      skipProcess(config$myAppName, 'BEDTools', 'intersect - capture region', dirname(methyl_result_file_bed))
    }
    
  }else{
    file.rename(methyl_result_file_rough, methyl_result_file_bed)
    if(file.exists(vcf_in))
      file.rename(vcf_in, vcf_out)
  }
  if(!checkIfFileExists(methyl_result_file_bed)) return(NULL)
  
  print(paste0("Saving methylation to RDS file: ", gsub(".bed",".rds", methyl_result_file_bed)))
  methyl_result_data <- data.table::fread(methyl_result_file_bed, header = F, sep = '\t')
  names(methyl_result_data) <- getMethylDataHeader(version = 2, size = 10)
  saveRDS(methyl_result_data, str_replace(methyl_result_file_bed, "methylation_results.bed", "methylation_results.rds"))
  
  stop_time_meth <- Sys.time()
  cat(paste0(config$myAppName, ": Panel Intersect (",config$meth_tool_current,") execution is finished. [", format(stop_time_meth - start_time_meth, digits=3) ,"]\n"))
  
  config$tmp_files <- c(config$tmp_files, methyl_result_file_rough, vcf_in)
  return(config)
}

#config_tools[config_tools$proces=="bssnper","temp_results_dirs"] <- "methyl_results/bssnper3"
##########################
######  run_BSSnper  #####
##########################
run_BSSnper <- function(config, config_tools){
  start_time_meth <- Sys.time()
  sample_basename       <- basename_sample(config$file_bam)
  bssnper_result_file   <- file.path(config$results_path, config_tools[config_tools$proces=="bssnper","temp_results_dirs"], sample_basename)
  bssnper_log_file_log  <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename, "_", config_tools[config_tools$process=="bssnper","logfile"]))
  dir.create(dirname(bssnper_result_file), showWarnings = FALSE)
  
  vcf_out                <- paste0(bssnper_result_file, '.sorted.vcf')
  vcf_rough_out          <- paste0(bssnper_result_file, '.rough.vcf')
  bssnper_candidate_out  <- paste0(bssnper_result_file, '.snp.candidate.out')
  meth_cg_out            <- paste0(bssnper_result_file, '.meth.cg')
  meth_chg_out           <- paste0(bssnper_result_file, '.meth.chg')
  meth_chh_out           <- paste0(bssnper_result_file, '.meth.chh')
  
  clipping_result_file     <- file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"], sample_basename)
  clipping_result_file_bam <- paste0(clipping_result_file, '.clipped.bam')
  
  if(config$overwrite_results | !(file.exists(vcf_out) & file.exists(meth_cg_out) & file.exists(meth_chg_out) & file.exists(meth_chh_out))){
    if(tolower(config$meth_processing) == tolower("allCHR")){
      ######  BS-Snper ALL CHR PROCESSING
      src_command <- command_BSSnper(config, clipping_result_file_bam, bssnper_candidate_out, meth_cg_out, meth_chg_out, meth_chh_out, vcf_out, log_append=">", bssnper_log_file_log)
      runSystemCommand(config$myAppName, 'BS-Snper', 'BS-Snper calling', src_command, config$verbose)
    }else if(startsWith(tolower(config$meth_processing),tolower("batchCHR"))){
      ######  BS-Snper CHR BATCH PROCESSING - SPLIT INPUT FILE (HDD)
      print(paste0("### BS-Snper CHR BATCH PROCESSING ###"))
      runSystemCommand(config$myAppName, 'rm', 'tmp folder', paste0("rm -rf ", config$tmp_data_path, "/*"), FALSE)
      #bamtools split -in file.bam -reference
      print(paste0("Splitting clipped bam file..."))
      src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools)," split -in ", clipping_result_file_bam," -reference -refPrefix SPLIT.")
      runSystemCommand(config$myAppName, 'bamtools', 'split', src_command, config$verbose)
      src_command <- paste0("mv ", clipping_result_file,".clipped.SPLIT.*.bam ", config$tmp_data_path,"/")
      runSystemCommand(config$myAppName, 'mv', 'split files', src_command, config$verbose)
      inputFiles <- list.files(config$tmp_data_path, full.names = T)
      sample_chr <- stringr::str_replace(basename(inputFiles), paste0(sample_basename,".clipped.SPLIT."),"")
      sample_chr <- stringr::str_replace(sample_chr,".bam","")
      #common processing for chr batch processing
      chr_output_files <- list()
      logAppend <- ">"
      #i=1
      for(i in 1:length(sample_chr)){
        print(paste0("#####################################################"))
        print(paste0("##### Running BS-Snper for chr: '",sample_chr[i],"' [",i,"/",length(sample_chr),"] #####"))
        current_chr <- sample_chr[i]
        current_input <- inputFiles[i]
        current_output <- add_file_prefix(bssnper_candidate_out, current_chr)
        chr_output_files[[i]] <- current_output
        src_command <- command_BSSnper(config, current_input, current_output, 
                                       add_file_prefix(meth_cg_out, current_chr), add_file_prefix(meth_chg_out, current_chr), add_file_prefix(meth_chh_out, current_chr),
                                       add_file_prefix(vcf_rough_out, current_chr),
                                       logAppend, bssnper_log_file_log)
        runSystemCommand(config$myAppName, 'BS-Snper', 'BS-Snper calling', src_command, config$verbose)
        logAppend <- ">>"
        if(file.exists(current_output))
          file.remove(current_input)
      }
      
      chr_output_files <- unlist(chr_output_files)
      if(!all(sapply(c(chr_output_files), checkIfFileExists))){
        cat(paste0(readLines(bssnper_log_file_log),"\n"))
        return(NULL)
      }
      
      print(paste0("Wrapping BS-Snper results..."))
      wrap_files(config, meth_cg_out, remove_input_files=T, ext = "cg", header_lines = 1)
      wrap_files(config, meth_chg_out, remove_input_files=T, ext = "chg", header_lines = 1)
      wrap_files(config, meth_chh_out, remove_input_files=T, ext = "chh", header_lines = 1)
      wrap_files(config, bssnper_candidate_out, remove_input_files=T, ext = "out", header_lines = 1)
      wrap_files(config, vcf_rough_out, remove_input_files=T, ext = "vcf", header_lines = "#")
      
      print(paste0("Sorting vcf file..."))
      ref_file_fai <- paste0(config$ref_data_path,"/",config$ref_data_sequence_file,".fai")
      if(file.exists(ref_file_fai)){
        src_command <- paste0("bedtools sort -header -i ", vcf_rough_out, " -faidx ", ref_file_fai, " > ", vcf_out)
      }else{
        src_command <- paste0("bedtools sort -header -i ", vcf_rough_out, " > ", vcf_out)
      }
      runSystemCommand(config$myAppName, 'bedtools', 'vcf sort', src_command, config$verbose)
      fileRemoveIfExists(vcf_rough_out)
    }
  }else{
    skipProcess(config$myAppName, 'BS-Snper', 'BS-Snper calling', dirname(bssnper_result_file))
  }
  
  if(!checkIfFileExists(vcf_out))               return(NULL)
  if(!checkIfFileExists(bssnper_candidate_out)) return(NULL)
  if(!checkIfFileExists(meth_cg_out))           return(NULL)
  if(!checkIfFileExists(meth_chg_out))          return(NULL)
  if(!checkIfFileExists(meth_chh_out))          return(NULL)
  
  stop_time_meth <- Sys.time()
  cat(paste0(config$myAppName, ": BS-Snper execution is finished. [", format(stop_time_meth - start_time_meth, digits=3) ,"]\n"))
  
  config$tmp_files <- c(config$tmp_files, bssnper_candidate_out, meth_cg_out, meth_chg_out, meth_chh_out)
  
  return(config)
}

############################################
######  command_BSSnper  #####
############################################
command_BSSnper <- function(config, in_file, out_file, meth_cg_out, meth_chg_out, meth_chh_out, vcf_out, log_append = ">", log_out){
  panel_command <- ""
  if(config$bssnper_intersect==T & str_trim(config$ref_data_intervals_file)!=""){
    panel_command <- paste0("--regions-file ", file.path(config$ref_data_path, config$ref_data_intervals_file))
  }
  src_command <- paste0("perl "      , file.path(config$tools_path, config$bssnper),
                        " ", in_file, 
                        " --fa "     , file.path(config$ref_data_path, config$ref_data_sequence_file),
                        " --output " , out_file,
                        " --methcg " , meth_cg_out,
                        " --methchg ", meth_chg_out,
                        " --methchh ", meth_chh_out,
                        " --minhetfreq 0.1 ",
                        " --minhomfreq 0.85 ",
                        " --minquali 15 ",
                        " --mincover ", config$min_depth, 
                        " --maxcover 1000 ",
                        " --minread2 2 ",
                        " --errorate 0.02 ",
                        " --mapvalue 20 ",
                        panel_command,
                        " > ", vcf_out,
                        " 2", log_append, " ", log_out)
  return(src_command)
}

##################################
######  run_BSSnperOld  #############
##################################
run_BSSnperOld <- function(config, config_tools){
  
  start_time_meth <- Sys.time()
  sample_basename       <- basename_sample(config$file_bam)
  bssnper_result_file   <- file.path(config$results_path, config_tools[config_tools$proces=="bssnper","temp_results_dirs"], sample_basename)
  bssnper_log_file_log  <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename, "_", config_tools[config_tools$process=="bssnper","logfile"]))
  dir.create(dirname(bssnper_result_file), showWarnings = FALSE)
  
  bssnper_candidate_out  <- paste0(bssnper_result_file, '.snp.candidate.out')
  meth_cg_out            <- paste0(bssnper_result_file, '.meth.cg' )
  meth_chg_out           <- paste0(bssnper_result_file, '.meth.chg' )
  meth_chh_out           <- paste0(bssnper_result_file, '.meth.chh' )
  vcf_out                <- paste0(bssnper_result_file, '.bssnper.vcf' )
  
  clipping_result_file     <- file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"], sample_basename)
  clipping_result_file_bam <- paste0(clipping_result_file, '.clipped.bam')
  
  if(config$overwrite_results | !(file.exists(vcf_out) & file.exists(meth_cg_out) & file.exists(meth_chg_out) & file.exists(meth_chh_out))){
    panel_command <- ""
    # if(str_trim(config$ref_data_intervals_file != "")){
    #   panel_command <- paste0("--regions-file ", file.path(config$ref_data_path, config$ref_data_intervals_file))
    # }
    src_command <- paste0("perl "      , file.path(config$tools_path, config$bssnper),
                          " ", clipping_result_file_bam, 
                          " --fa "     , file.path(config$ref_data_path, config$ref_data_sequence_file),
                          " --output " , bssnper_candidate_out,
                          " --methcg " , meth_cg_out,
                          " --methchg ", meth_chg_out,
                          " --methchh ", meth_chh_out,
                          " --minhetfreq 0.1",
                          " --minhomfreq 0.85",
                          " --minquali 15",
                          " --mincover", config$min_depth, 
                          " --maxcover 1000",
                          " --minread2 2",
                          " --errorate 0.02",
                          " --mapvalue 20",
                          panel_command,
                          " > ", vcf_out,
                          " 2> ", bssnper_log_file_log
    )
    runSystemCommand(config$myAppName, 'BS-Snper', 'BS-Snper calling', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'BS-Snper', 'BS-Snper calling', dirname(bssnper_result_file))
  }

  if(!checkIfFileExists(bssnper_candidate_out)) return(NULL)
  if(!checkIfFileExists(meth_cg_out))           return(NULL)
  if(!checkIfFileExists(meth_chg_out))         return(NULL)
  if(!checkIfFileExists(meth_chh_out))         return(NULL)
  
  stop_time_meth <- Sys.time()
  cat(paste0(config$myAppName, ": BSsnper execution is finished. [", format(stop_time_meth - start_time_meth, digits=3) ,"]\n"))
  
  config$tmp_files <- c(config$tmp_files, bssnper_candidate_out, meth_cg_out, meth_chg_out, meth_chh_out)
  
  return(config)
}

############################################
######  bssnper calculate methylation  #####
############################################
bssnper_to_bed <- function(config, config_tools){
  
  start_time_meth <- Sys.time()
  sample_basename <- basename_sample(config$file_bam)
  bssnper_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="bssnper","temp_results_dirs"], sample_basename)
  bssnper_result_file_bed <- paste0(bssnper_result_file, ".methylation_results.rough.bed")
  
  if ( !(config$overwrite_results | !file.exists(bssnper_result_file_bed) )) {
    skipProcess(config$myAppName, 'BS-Snper', 'Calculate Methylation ', dirname(bssnper_result_file_bed))
  }else{
    ref_fasta <- file.path(config$ref_data_path, config$ref_data_sequence_file)
    
    meth_data <- list()
    meth_data$f1 <- data.table::fread(paste0(bssnper_result_file, '.meth.cg' ), na.strings='.',  header = T, sep = '\t')
    meth_data$f2 <- data.table::fread(paste0(bssnper_result_file, '.meth.chg'), na.strings='.',  header = T, sep = '\t')
    meth_data$f3 <- data.table::fread(paste0(bssnper_result_file, '.meth.chh'), na.strings='.',  header = T, sep = '\t')
    meth_data <- data.table::rbindlist(meth_data)
    meth_data        <- meth_data[,-c('Watson-QUAL', 'Crick-QUAL')]
    names(meth_data) <- c('chr', 'pos', 'context', 'Wmeth', 'Wcov', 'Cmeth', 'Ccov')
    
    meth_data$bv     <- NA
    meth_data$start  <- meth_data$pos
    meth_data$end    <- meth_data$pos
    meth_data$strand <- "."
    meth_data$cov    <- NA
    meth_data$met    <- NA
    meth_data <- meth_data[,c("chr","start","end","context", "bv", "strand","Wmeth","Wcov","Cmeth","Ccov","pos","met","cov")]

    ##### Porcess CG, CHG, CHH contexts
    cg_w <- cg_c <-  meth_data[meth_data$context=="CG"]
    chg_w  <- meth_data[meth_data$context=="CHG" & is.na(meth_data$Cmeth)]
    chg_c  <- meth_data[meth_data$context=="CHG" & is.na(meth_data$Wmeth)]
    chg_wc1 <- chg_wc2  <- meth_data[meth_data$context=="CHG" & !is.na(meth_data$Cmeth) & !is.na(meth_data$Wmeth)]
    chh_w  <- meth_data[meth_data$context=="CHH" & is.na(meth_data$Cmeth)]
    chh_c  <- meth_data[meth_data$context=="CHH" & is.na(meth_data$Wmeth)]
    rm(meth_data)
    
    #### GC PROCESSING
    cat("##processing CG context\n")
    ##GC Watson
    cg_w$start  <- cg_w$pos-1
    cg_w$bv     <- cg_w$Wmeth / cg_w$Wcov
    cg_w        <- cg_w[!is.na(cg_w$bv)]
    cg_w$strand <- "+"
    cg_w$cov    <- cg_w$Wcov
    cg_w$met    <- cg_w$Wmeth
    
    ## TODO: here we can do some automatic test(??)
    #checkLetters(cg_w,ref_fasta) 
    
    ##GC Crick
    # cg_c <- bssnper_to_bval(meth_data, "c","-", 0, "CG", "+", 1)
    cg_c$end    <-  cg_c$pos+1
    cg_c$bv     <-  cg_c$Cmeth / cg_c$Ccov
    cg_c        <-  cg_c[!is.na(cg_c$bv)]
    cg_c$strand <- "-"
    cg_c$cov    <- cg_c$Ccov
    cg_c$met    <- cg_c$Cmeth
    #checkLetters(cg_c, ref_fasta)
    
    ######## CHG PROCESSING
    cat("##processing CHG context\n")
    ##CHG Watson
    chg_w$start  <- chg_w$pos-1
    chg_w$bv     <- chg_w$Wmeth / chg_w$Wcov
    chg_w$strand <- "+"
    chg_w$cov    <- chg_w$Wcov
    chg_w$met    <- chg_w$Wmeth
    #checkLetters(chg_w, ref_fasta)
    
    ##CHG CRICK
    chg_c$start  <- chg_c$pos-1
    chg_c$bv     <- chg_c$Cmeth / chg_c$Ccov
    chg_c$strand <- "-"
    chg_c$cov    <- chg_c$Ccov
    chg_c$met    <- chg_c$Cmeth
    #checkLetters(chg_c,ref_fasta)
    
    ##CHG Watson-Crick
    chg_wc1$start  <- chg_wc1$pos-1
    chg_wc1$bv     <- chg_wc1$Wmeth / chg_wc1$Wcov
    chg_wc1$strand <- "+"
    chg_wc1$cov    <- chg_wc1$Wcov
    chg_wc1$met    <- chg_wc1$Wmeth
    #checkLetters(chg_wc1,ref_fasta)
    
    ##CHG Crick-Watson 
    chg_wc2$start  <- chg_wc2$pos+1
    chg_wc2$end    <- chg_wc2$pos+2
    chg_wc2$bv     <- chg_wc2$Cmeth / chg_wc2$Ccov
    chg_wc2$strand <- "-"
    chg_wc2$cov    <- chg_wc2$Ccov
    chg_wc2$met    <- chg_wc2$Cmeth
    #checkLetters(chg_wc2,ref_fasta)
    
    ######## CHH PROCESSING
    cat("##processing CHH context\n")
    ##CHH Watson
    chh_w$start  <- chh_w$pos-1
    chh_w$bv     <- chh_w$Wmeth / chh_w$Wcov
    chh_w$strand <- "+"
    chh_w$cov    <- chh_w$Wcov
    chh_w$met    <- chh_w$Wmeth
    #checkLetters(chh_w,ref_fasta)
    
    ##CHH Crick
    chh_c$start   <- chh_c$pos-1
    chh_c$bv      <- chh_c$Cmeth / chh_c$Ccov
    chh_c$strand  <- "-"
    chh_c$cov     <- chh_c$Ccov
    chh_c$met     <- chh_c$Cmeth
    #checkLetters(chh_c,ref_fasta)
    
    #### methylation result: bssnper
    meth_data <- data.table::rbindlist(list(cg_w, cg_c, chg_w, chg_c, chg_wc1, chg_wc2, chh_w, chh_c))[,-c("Wmeth", "Wcov", "Cmeth", "Ccov")]
    rm(cg_w, cg_c, chg_w, chg_c, chg_wc1, chg_wc2, chh_w, chh_c)
    meth_data <- meth_data[order(chr, start,end),]
    meth_data$bv <- round(meth_data$bv,3)
    
    ### make the names the same like in Methratio output
    meth_data$pos <- meth_data$end ## make sure it's the same
    meth_data$coverage <- meth_data$cov
    meth_data     <- meth_data[,c("chr","start","end","context", "bv", "strand", "coverage","met","cov","pos" )]
    names(meth_data) <-  getMethylDataHeader(version = 2, size = 10)

    print(paste0("Saving the file: ", bssnper_result_file_bed))
    data.table::fwrite(meth_data, file = bssnper_result_file_bed, quote=FALSE, sep='\t', row.names = F, col.names = F)

    stop_time_meth <- Sys.time()
    cat(paste0(config$myAppName, ": bssnper_to_bed execution is finished. [", format(stop_time_meth - start_time_meth, digits=3) ,"]\n"))
    
    config$tmp_files <- c(config$tmp_files, bssnper_result_file_bed)
  }
  
  return(config)
}

#############################################
######  compare BS-Snper and Methratio #######
#############################################
## this function is for testing purposes, it is not part of the main pipeline
 
# sample_basename <- basename_sample(config$file_bam)
# methyl_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"], sample_basename)
# bssnper_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="bssnper","temp_results_dirs"], sample_basename)
# meth_data_src  <- readRDS(paste0(methyl_result_file,".methylation_results.rds"))
# attr(meth_data_src, "name") <- "methratio"
# meth_data_dst  <- readRDS(paste0(bssnper_result_file, ".methylation_results.rds"))
# attr(meth_data_dst, "name") <- "bssnper"
# 
# compareMethData(meth_data_src, meth_data_dst)
 
compareMethData <- function(meth_data_src, meth_data_dst){

  meth_data_src <- setDT(meth_data_src)[order(chr, start, end),][,c("numCs","numTs","posCs"):=NULL]
  meth_data_dst <- setDT(meth_data_dst)[order(chr, start, end),][,c("numCs","numTs","posCs"):=NULL]

  chr_cnt_src <- data.table(table(meth_data_src$chr,meth_data_src$context))
  #names(chr_cnt_src) <- c("chr", "context", attr(meth_data_src, "name"))
  names(chr_cnt_src) <- c("chr", "context", "src")
  chr_cnt_dst <- data.table(table(meth_data_dst$chr,meth_data_dst$context))
  #names(chr_cnt_dst) <- c("chr", "context", attr(meth_data_dst, "name"))
  names(chr_cnt_dst) <- c("chr", "context", "dst")
  
  chr_cnt <- dplyr::full_join(chr_cnt_src, chr_cnt_dst, by=c("chr","context"))
  chr_cnt$diff <- chr_cnt[,4]-chr_cnt[,3]
  chr_cnt <- chr_cnt[order(chr,context)]
  chr_cnt <- split(chr_cnt, chr_cnt$context)
  
  #chr_cnt
  
  names(meth_data_src) <- paste0("src_",names(meth_data_src))
  names(meth_data_dst) <- paste0("dst_",names(meth_data_dst))
  meth_data_full <- dplyr::full_join(meth_data_src, meth_data_dst, by=c('src_chr'='dst_chr','src_start'='dst_start','src_end'='dst_end','src_strand'='dst_strand'))
  #meth_data_full

  
  meth_compare <- data.table(metric="meth_sum", value = nrow(meth_data_full))
  meth_compare <- rbind(meth_compare, data.table(metric="meth_common", value = sum(!is.na(meth_data_full$src_betaVal) & !is.na(meth_data_full$dst_betaVal)) ))
  #null in dst then unique src
  meth_compare <- rbind(meth_compare, data.table(metric="meth_unique_src", value = sum(is.na(meth_data_full$dst_betaVal)) ))
  meth_compare <- rbind(meth_compare, data.table(metric="meth_unique_dst", value = sum(is.na(meth_data_full$src_betaVal)) ))
  meth_compare$pct <- 100 * meth_compare$value/meth_compare$value[meth_compare$metric=="meth_sum"]
  meth_compare <- rbind(meth_compare, data.table(metric="meth_common_pearson", value = cor(meth_data_full$src_betaVal, meth_data_full$dst_betaVal, use = "complete.obs"), pct = NA ))
  meth_data_full$beta_diff <- meth_data_full$src_betaVal - meth_data_full$dst_betaVal
  beta_diff <- meth_data_full$beta_diff[!is.na(meth_data_full$beta_diff)]
  meth_compare <- rbind(meth_compare, data.table(metric="meth_diff_0", value = sum(abs(beta_diff)==0), pct = 100* sum(abs(beta_diff)==0)/meth_compare$value[meth_compare$metric=="meth_common"] ))
  meth_compare <- rbind(meth_compare, data.table(metric="meth_diff_less_0.1", value = sum(abs(beta_diff)<0.1), pct = 100* sum(abs(beta_diff)<0.1)/meth_compare$value[meth_compare$metric=="meth_common"] ))
  meth_compare <- rbind(meth_compare, data.table(metric="meth_diff_less_0.5", value = sum(abs(beta_diff)<0.5), pct = 100* sum(abs(beta_diff)<0.5)/meth_compare$value[meth_compare$metric=="meth_common"] ))
  meth_compare <- rbind(meth_compare, data.table(metric="meth_diff_less_0.9", value = sum(abs(beta_diff)<0.9), pct = 100* sum(abs(beta_diff)<0.9)/meth_compare$value[meth_compare$metric=="meth_common"] ))
  #meth_compare

  meth_data_diff <- meth_data_full[is.na(meth_data_full$beta_diff) | meth_data_full$beta_diff>=0.1]

  gg_beta_diff <- ggplot(dplyr::sample_n(data.table(beta_diff=beta_diff), size=min(c(100000, length(beta_diff))) ) , aes(x=beta_diff)) + 
    geom_histogram(aes(y=..density..), breaks=seq(-1,1,0.1), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666") 

  retlist <- list(chr_cnt = chr_cnt, meth_compare = meth_compare, meth_data_diff = meth_data_diff, gg_beta_diff = gg_beta_diff)
  
  return(retlist)
}


