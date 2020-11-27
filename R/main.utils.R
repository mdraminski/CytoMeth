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
    skipProcess(config$myAppName, 'seqtk', 'subsample of reads',
                file.path(config$results_path, config_tools[config_tools$proces=="seqtk","temp_results_dirs"],"/"))
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
    
    src_command <- paste0(file.path(config$anaconda_bin_path, config$fastqc), " --nogroup ", 
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
    skipProcess(config$myAppName, 'Trimmomatic', 'trimming',
                file.path(config$results_path, config_tools[config_tools$proces=="trimming","temp_results_dirs"],"/"))
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
  mapping_dir <- file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"],"/")
  
  if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".sam")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path,config$bsmap), " -r 0 -s 16 -n 1 ",
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
  mapping_dir <- file.path(config$results_path, config_tools[config_tools$proces=="mapping","temp_results_dirs"],"/")
  
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
    src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " sort",
                          " -in ", paste0(mapping_result_file,".top.bam"), 
                          " -out ", paste0(mapping_result_file,".top.bam.sorted"))
    runSystemCommand(config$myAppName, 'Bamtools', 'sort top', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Bamtools', 'sort top', mapping_dir)
  }
  if(!checkIfFileExists(paste0(mapping_result_file, ".top.bam.sorted"))) return(NULL)
  
  ######  bamtools - sort bottom
  if(config$overwrite_results | !(file.exists(paste0(mapping_result_file, ".bottom.bam.sorted")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path,config$bamtools), " sort ",
                          " -in ", paste0(mapping_result_file,".bottom.bam"), 
                          " -out ", paste0(mapping_result_file,".bottom.bam.sorted"))
    runSystemCommand(config$myAppName, 'Bamtools', 'sort bottom', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Bamtools', 'sort bottom', mapping_dir)
  }
  if(!checkIfFileExists(paste0(mapping_result_file, ".bottom.bam.sorted"))) return(NULL)
  
  ######  picard - MarkDuplicates - top
  rmdups_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_top","temp_results_dirs"], sample_basename)
  rmdups_logfile <- file.path(config$results_path, "logs", paste0(sample_basename, "_", config_tools[config_tools$process=="mark_duplicates_top","logfile"]))
  rmdups_dir <- file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_top","temp_results_dirs"])
  if(config$overwrite_results | !(file.exists(paste0(rmdups_result_file, ".top.rmdups.bam")) & file.exists(paste0(rmdups_result_file, ".top.rmdups_metrics.txt")) &
                                  file.exists(rmdups_logfile)
  )) {
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
    skipProcess(config$myAppName, 'Bamtools', 'merge',
                file.path(config$results_path, config_tools[config_tools$proces=="merge2_bam","temp_results_dirs"],"/"))
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
    skipProcess(config$myAppName, 'Bamtools', 'filter',
                file.path(config$results_path, config_tools[config_tools$proces=="filter_bam","temp_results_dirs"],"/"))
  }
  if(!checkIfFileExists(paste0(filtered_result_file,".filtered.bam"))) return(NULL)
  
  ######  Samtools - flagstat flagstat_filtered_bam
  flagstat_result_file <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename,"_",config_tools[config_tools$process=="flagstat_filtered_bam","logfile"]))
  if(config$overwrite_results | !(file.exists(paste0(flagstat_result_file)))) {
    src_command <- paste0(file.path(config$anaconda_bin_path, config$samtools), " flagstat ",  
                          paste0(filtered_result_file,".filtered.bam"), 
                          " 1> ", flagstat_result_file)
    runSystemCommand(config$myAppName, 'Samtools', 'flagstat', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Samtools', 'flagstat', file.path(config$results_path,"logs","/"))
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
    skipProcess(config$myAppName, 'BamUtil', 'clipOverlap',
                file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"],"/"))
  }
  
  if(!checkIfFileExists(paste0(clipping_result_file,".clipped.bam"))) return(NULL)
  if(!checkLog(clipOverlap_logfile, "Completed ClipOverlap Successfully.", 'BamUtil - clipOverlap')) return(NULL)

  ######  samtools - index
  if(config$overwrite_results | !(file.exists(paste0(clipping_result_file,".clipped.bam.bai")))) {
    src_command <- paste0(file.path(config$anaconda_bin_path, config$samtools), " index ", paste0(clipping_result_file,".clipped.bam"))
    runSystemCommand(config$myAppName, 'Samtools', 'index', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Samtools', 'index',
                file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"],"/"))
  }
  
  if(!checkIfFileExists(paste0(clipping_result_file,".clipped.bam.bai"))) return(NULL)
    
  return(config)
}

##################################
######  run_MappingMetrics  ######
##################################
run_MappingMetrics <- function(config, config_tools){
  sample_basename <- basename_sample(config$file_bam)
  filtered_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="filter_bam","temp_results_dirs"], sample_basename)
  rmdups_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_top","temp_results_dirs"], sample_basename) ### TODO ?? small_FAKE03.rmdups.bam not small_FAKE03.clipped.bam

  ######  picard - Basic Mapping Metrics
  basic_mapping_metrics_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="basic_mapping_metrics","temp_results_dirs"], sample_basename)
  if(config$overwrite_results | !(file.exists(paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt")))) {
    src_command <- paste0("java -Xms", f2si(config$memory)," -Xmx", f2si(config$memory), 
                          config$picard_parser,
                          " -jar ", file.path(config$tools_path,config$picard), 
                          " CollectAlignmentSummaryMetrics ",
                          " VALIDATION_STRINGENCY=LENIENT",
                          " METRIC_ACCUMULATION_LEVEL=ALL_READS ",
                          " INPUT=", paste0(rmdups_result_file,".rmdups.bam"), 
                          " OUTPUT=", paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt"),
                          " REFERENCE_SEQUENCE=", file.path(config$ref_data_path, config$ref_data_sequence_file))
    runSystemCommand(config$myAppName, 'Picard', 'Basic Mapping Metrics', src_command, config$verbose)
  }else{
    skipProcess(config$myAppName, 'Picard', 'Basic Mapping Metrics',
                file.path(config$results_path, config_tools[config_tools$proces=="basic_mapping_metrics","temp_results_dirs"],"/"))
  }
  if(!checkIfFileExists(paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt"))) return(NULL)

  ######  picard - Insert Size Metrics
  insert_size_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="insert_size","temp_results_dirs"], sample_basename)  ### TODO ??? small_FAKE03.filtered.bam? Shouldn't be .clipped.bam?
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
    skipProcess(config$myAppName, 'Picard', 'Insert Size',
                file.path(config$results_path, config_tools[config_tools$proces=="insert_size","temp_results_dirs"],"/"))
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
      runSystemCommand(config$myAppName, 'BEDTools', 'Intersect on_target_reads', src_command, config$verbose)  ### TODO ??? small_FAKE03.rmdups.bam? Shouldn't be .clipped.bam?
    }else{
      ## If the panel file is not defined use convert bam2bed - (bedtools bamtobed -i reads.bam)
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bedtools), 
                            " bamtobed -i ",  paste0(rmdups_result_file, ".rmdups.bam"),
                            " > ",paste0(on_target_result_file, "_on_target_reads"))
      runSystemCommand(config$myAppName, 'BEDTools', 'bamtobed on_target_reads', src_command, config$verbose)
    }
  }else{
    skipProcess(config$myAppName, 'BEDTools', 'Intersect on_target_reads',
                file.path(config$results_path, config_tools[config_tools$proces=="on_target_reads","temp_results_dirs"],"/"))
  }
  if(!checkIfFileExists(paste0(on_target_result_file,"_on_target_reads"))) return(NULL)

  #update and save result SAMPLE_on_target_reads.txt file
  number_of_on_target_reads <- getLinesNumber(paste0(on_target_result_file,"_on_target_reads"))
  on_target_params <- list(number_of_on_target_reads = as.character(number_of_on_target_reads))
  yaml::write_yaml(on_target_params, paste0(on_target_result_file,"_on_target_reads.txt"))

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
    skipProcess(config$myAppName, 'GATK', 'depth_of_coverage',
                file.path(config$results_path, config_tools[config_tools$proces=="depth_of_coverage","temp_results_dirs"],"/"))
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
  methyl_result_file <- paste0(file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"], sample_basename),".methylation_results.txt")
  methratio_logfile <- file.path(file.path(config$results_path,"logs"),paste0(sample_basename,"_",config_tools[config_tools$process=="methratio","logfile"]))
  
  if(config$overwrite_results | !(file.exists(methyl_result_file) & file.exists(methratio_logfile))) {
    ######  Methratio ALL CHR PROCESSING
    if(tolower(config$methratio_processing) == tolower("allCHR")){
      src_command <- paste0(config$python2, " ", file.path(config$tools_path, config$methratio), 
                            " -d ", file.path(config$ref_data_path, config$ref_data_sequence_file), 
                            #" -s ", config$anaconda_bin_path,
                            " -m 1 -z -i skip", 
                            " -o ", methyl_result_file," ", clipping_result_file_bam,
                            " 2> ", methratio_logfile)
      runSystemCommand(myAppName, 'BSMAP', 'methratio', src_command, config$verbose)
    }else if(startsWith(tolower(config$methratio_processing),tolower("batchCHR"))){
      ######  Methratio CHR BATCH PROCESSING - SPLIT INPUT FILE (HDD)
      if(endsWith(tolower(config$methratio_processing),tolower("HDD"))){
        batchHDD <- TRUE
        print(paste0("Conversion of input clipped.bam to clipped.sam [bamToSam]..."))
        src_command <- paste0(file.path(config$anaconda_bin_path,"samtools")," view ", clipping_result_file_bam," > ",clipping_result_file_sam)
        runSystemCommand(config$myAppName, 'samtools', 'bamView', src_command, config$verbose)
        print(paste0("Splitting clipped sam file..."))
        src_command <- paste0("cd ", dirname(clipping_result_file_sam), "; awk '{print>$3}' ", basename(clipping_result_file_sam), ";")
        runSystemCommand(config$myAppName, 'awk', 'splitclippedbed', src_command, config$verbose)
        file.remove(clipping_result_file_sam)
        inputFiles <- list.files(dirname(clipping_result_file_sam),full.names = T)
        inputFiles <- sapply(inputFiles, function(x) {file.rename(x,paste0(x,".bsam"))})
        inputFiles <- list.files(dirname(clipping_result_file_sam),full.names = T)
        sample_chr <- stringr::str_replace(basename(inputFiles),".bsam","")
      }else{
      ######  Methratio CHR BATCH PROCESSING - NO SPLITTING
        batchHDD <- FALSE
        print(paste0("Conversion of input clipped.bam to clipped.bed [bamToBed]..."))
        src_command <- paste0(file.path(config$anaconda_bin_path,"bamToBed")," -i ", clipping_result_file_bam," > ", clipping_result_file_bed)
        runSystemCommand(config$myAppName, 'bedtools', 'bamToBed', src_command, config$verbose)
        print(paste0("Reading sample chromosomes..."))
        sample_chr <- getSampleChr(clipping_result_file_bed)
        file.remove(clipping_result_file_bed)
        print(head(sample_chr, 12))
        inputFiles <- rep(clipping_result_file_bam, length(sample_chr))
      }
      #common processing for chr batch processing
      chr_output_files <- list()
      logAppend <- ">"
      for(i in 1:length(sample_chr)){
        print(paste0("#####################################################"))
        print(paste0("##### Running methratio for chr: '",sample_chr[i],"' [",i,"/",length(sample_chr),"] #####"))
        current_chr <- sample_chr[i]
        current_input <- inputFiles[i]
        current_output <- stringr::str_replace(methyl_result_file, ".methylation_results.txt", paste0(".methylation_results.",sample_chr[i],".txt"))
        chr_output_files[[i]] <- current_output
        src_command <- paste0(config$python2, " ", file.path(config$tools_path, config$methratio), 
                              " -c ", current_chr," -d ", file.path(config$ref_data_path, config$ref_data_sequence_file), 
                              #" -s ", config$anaconda_bin_path,
                              " -m 1 -z -i skip", 
                              " -o ", current_output, " ", current_input,
                              " 2", logAppend, methratio_logfile)
                              #" 1", logAppend, methratio_logfile, " 2>&1")
        runSystemCommand(config$myAppName, 'BSMAP', 'methratio', src_command, config$verbose)
        logAppend <- ">>"
        if(batchHDD & file.exists(current_output))
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
      data.table::fwrite(methyl_result_data, methyl_result_file, sep = '\t')
    }
  }else{
    skipProcess(config$myAppName, 'BSMAP', 'methratio',
                file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"],"/"))
  }
  
  if(!checkIfFileExists(methratio_logfile)) return(NULL)
  if(!checkIfFileExists(methyl_result_file)) {
    cat(paste0(readLines(methratio_logfile),"\n"))
    return(NULL)
  }
  stop_time_meth <- Sys.time()
  cat(paste0(config$myAppName, ": All methratio Calculations are finished. [", format(stop_time_meth - start_time_meth, digits=3) ,"]\n"))
  
  config$tmp_files <- c(config$tmp_files, clipping_result_file_bed, clipping_result_file_sam, chr_output_files)
  return(config)
}

##################################
######  run_CalcMethylation  #####
##################################
run_CalcMethylation <- function(config, config_tools){

  sample_basename <- basename_sample(config$file_bam)
  ######  python - methratio (BSMAP)
  methyl_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"], sample_basename)
  methyl_result_file_txt <- paste0(methyl_result_file,".methylation_results.txt")
  methyl_result_file <- paste0(methyl_result_file,".methylation_results.bed")
  
  ######  Make bed file from methyl results
  print(paste0("File conversion txt -> bed..."))
  print(paste0("Reading the file: ", methyl_result_file_txt))
  methyl_result_data <- data.table::fread(methyl_result_file_txt, header = T, sep = '\t')
  
  methyl_result_data$end <- methyl_result_data$pos
  methyl_result_data <- methyl_result_data[,c("chr","pos","end","context","ratio","strand","eff_CT_count","C_count", "CT_count")]
  names(methyl_result_data) <- getMethylDataHeader(version = 2, size = 9)
  #head(methyl_result_data)
  #summary(methyl_result_data)
  if(nrow(methyl_result_data) == 0){
    print(paste0("Error! methyl_result_data is empty!"))
    print(paste0("QC Report of the Sample: ", sample_basename))
    sampleQC <- getSampleQCSummary(sample_basename, config)
    print(qc2dataframe(sampleQC))
    return(NULL)
  }else{
    print(paste0("Number of numCs > 0: ", length(methyl_result_data$numCs[methyl_result_data$numCs > 0])))
    print(head(methyl_result_data[methyl_result_data$numCs > 0,],5))
    print(paste0("Number of numCs > 1: ", length(methyl_result_data$numCs[methyl_result_data$numCs > 1])))
    print(head(methyl_result_data[methyl_result_data$numCs > 1,],5))
  }
  methyl_result_prime_file <- str_replace(methyl_result_file,".methylation_results.bed",".methylation_results.prime.bed")
  print(paste0("Saving the file: ",methyl_result_prime_file))
  data.table::fwrite(methyl_result_data, file = methyl_result_prime_file, quote=FALSE, sep='\t', row.names = F, col.names = F)
  
  ######  BEDTools - intersect capture region (only if file is defined)
  if(str_trim(config$ref_data_intervals_file) != ''){
    if(config$overwrite_results | !(file.exists(methyl_result_file))) {
      src_command <- paste0(file.path(config$anaconda_bin_path, config$bedtools), 
                            " intersect -bed -abam ",  methyl_result_prime_file, 
                            " -b ", file.path(config$ref_data_path, config$ref_data_intervals_file),
                            " -u > ", methyl_result_file)
      runSystemCommand(config$myAppName, 'BEDTools', 'intersect - capture region', src_command, config$verbose)
    }else{
      skipProcess(config$myAppName, 'BEDTools', 'intersect - capture region',
                  file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"],"/"))
    }
    if(!checkIfFileExists(methyl_result_file)) return(NULL)
  }else{
    file.rename(methyl_result_prime_file, methyl_result_file)
  }

  print(paste0("Final conversion of methylations..."))
  #finish the methylation file and save rds version
  methyl_result_data <- readMethResult(methyl_result_file, version = 2)
  methyl_result_data$posCs <- methyl_result_data$start 
  # lets make end = start
  methyl_result_data$end <- methyl_result_data$start
  # for positive strand end++ (end = start + 1)
  methyl_result_data$end[methyl_result_data$strand == '+'] <- methyl_result_data$end[methyl_result_data$strand == '+'] + 1
  # for negative strand start-- (end = start - 1)
  methyl_result_data$start[methyl_result_data$strand == '-'] <- methyl_result_data$start[methyl_result_data$strand == '-'] - 1
  
  print(paste0("Saving final methylation files: ",methyl_result_file, " and ", gsub(".bed",".rds", methyl_result_file)))
  #save bed file
  data.table::fwrite(methyl_result_data, methyl_result_file, quote=FALSE, sep='\t', row.names = F, col.names = F, scipen=50)
  #save rds file
  saveRDS(methyl_result_data, str_replace(methyl_result_file, "methylation_results.bed","methylation_results.rds"))
  
  config$tmp_files <- c(config$tmp_files, methyl_result_file_txt, methyl_result_prime_file)
  return(config)
}



##################################
######  run_BSsnper  #############
##################################
run_BSsnper <- function(config, config_tools){
  
  meth_cg_suf  <- '.meth.cg'
  meth_chg_suf <- '.meth.chg'
  meth_chh_suf <- '.meth.chh'
  vcf_suf      <- '.bssnper.vcf'
  
  sample_basename     <- basename_sample(config$file_bam)
  bssnper_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="bssnper","temp_results_dirs"], sample_basename)
  bssnper_log_file    <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename, "_", config_tools[config_tools$process=="bssnper","logfile"]))
  bssnper_candidate_snps_file <- paste0(bssnper_result_file, '.snp.candidate.out')
  clipping_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="clip_overlap","temp_results_dirs"], sample_basename)
  
  if(config$overwrite_results | 
     !(file.exists(paste0(bssnper_result_file, vcf_suf )) 
       & file.exists(paste0(bssnper_result_file, meth_cg_suf ))
       & file.exists(paste0(bssnper_result_file, meth_chg_suf)) 
       & file.exists(paste0(bssnper_result_file, meth_chh_suf)) 
      )
     ){
    
    panel_command <- ""
    
    if(str_trim(config$ref_data_intervals_file != "")){
      panel_command <- paste0("--regions-file ", file.path(config$ref_data_path, config$ref_data_intervals_file))
    }
    
    src_command <- paste0("perl "      , file.path(config$tools_path, config$bssnper),
                          " --fa "     , file.path(config$ref_data_path, config$ref_data_sequence_file),
                          " --input "  , paste0(clipping_result_file, '.clipped.bam'), 
                          " --output " , bssnper_candidate_snps_file,
                          " --methcg " , paste0(bssnper_result_file, meth_cg_suf ),
                          " --methchg ", paste0(bssnper_result_file, meth_chg_suf),
                          " --methchh ", paste0(bssnper_result_file, meth_chh_suf),
                          " --minhetfreq 0.1 ",
                          " --minhomfreq 0.85 ",
                          " --minquali 15 ",
                          " --mincover 10 ",
                          " --maxcover 1000 ",
                          " --minread2 2 ",
                          " --errorate 0.02 ",
                          " --mapvalue 20 ",
                          panel_command,
                          " > ", paste0(bssnper_result_file, vcf_suf),
                          " 2> ", bssnper_log_file
                          )
    
    runSystemCommand(config$myAppName, 'BS-Snper', 'BS-Snper calling', src_command, config$verbose)
    
  }else{
    skipProcess(config$myAppName, 'BS-Snper', 'BS-Snper calling',
                file.path(config$results_path, config_tools[config_tools$proces=="bssnper","temp_results_dirs"],"/"))
  }
  
  config$tmp_files <- c(config$tmp_files, bssnper_candidate_snps_file)
  return(config)
}










