################################
###### get and save SAMPLE_QC_summary ####
#config <- conf; save = F; result_format = 'bed'; sample_basename <- "DA01"; sample_basename <- "small_FAKE01"; 
getSampleQCSummary <- function(sample_basename, config, result_format = c('bed','rds')){
  result_format <- checkFormat(result_format, supported = c('bed','rds'))
  config_tools <- read.csv(file.path(config$tools_path, config$tools_config), stringsAsFactors = FALSE)
  sampleQC <- list()
  sampleQC$Sample_ID <- sample_basename
  
  ### trimmomatic
  trimmomatic_out_logfile <- file.path(config$results_path,"logs",paste0(sample_basename,"_",config_tools[config_tools$tool=="trimmomatic","logfile"]))
  if(file.exists(trimmomatic_out_logfile)){
    trim_log <- readLines(trimmomatic_out_logfile)
    trim_log <- str_split(trim_log[startsWith(trim_log, "Input Read Pairs")], " ")
    sampleQC$Input_read_pairs <- as.double(trim_log[[1]][4])
    sampleQC$Read_Pairs_Surviving_trimming <- as.double(trim_log[[1]][7])
    sampleQC$Prc_Read_Pairs_Surviving_trimming <- 100 * sampleQC$Read_Pairs_Surviving_trimming/sampleQC$Input_read_pairs
  }else{
    print(paste0("File: ", trimmomatic_out_logfile, " does not exist!"))
  }
  
  ### top_rmdups
  rmdups_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="mark_duplicates_top", "temp_results_dirs"], sample_basename)
  top_rmdups_result_file <- file.path(paste0(rmdups_result_file,".top.rmdups_metrics.txt"))
  
  if(file.exists(top_rmdups_result_file)){
    top_dup <- readLines(top_rmdups_result_file)
    top_dup <- makeDataFrame(data = unlist(str_split(top_dup[startsWith(top_dup, "SAMPLE")], "\t")), 
                             header = tolower(unlist(str_split(top_dup[startsWith(top_dup, "LIBRARY")], "\t"))))
    sampleQC$Prc_duplicated_reads_top <- as.numeric(top_dup[,"percent_duplication"])*100
  }else{
    print(paste0("File: ", top_rmdups_result_file, " does not exist!"))
  }
  
  ### bottom_rmdups
  bottom_rmdups_result_file <- file.path(paste0(rmdups_result_file,".bottom.rmdups_metrics.txt"))
  if(file.exists(bottom_rmdups_result_file)){
    bottom_dup <- readLines(bottom_rmdups_result_file)
    bottom_dup <- makeDataFrame(data = unlist(str_split(bottom_dup[startsWith(bottom_dup, "SAMPLE")], "\t")), 
                                header = tolower(unlist(str_split(bottom_dup[startsWith(bottom_dup, "LIBRARY")], "\t"))))
    sampleQC$Prc_duplicated_reads_bottom <- as.numeric(bottom_dup[,"percent_duplication"])*100
  }else{
    print(paste0("File: ", bottom_rmdups_result_file, " does not exist!"))
  }
  
  ### basic_mapping
  basic_mapping_metrics_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="basic_mapping_metrics","temp_results_dirs"], sample_basename)
  basic_mapping_metrics_result_file <- paste0(basic_mapping_metrics_result_file,"_basic_mapping_metrics.txt")
  if(file.exists(basic_mapping_metrics_result_file)){
    metrics <- readLines(basic_mapping_metrics_result_file)
    headerStartLine <- which(startsWith(metrics,"CATEGORY"))
    metrics_df <- makeDataFrame(data = unlist(str_split(metrics[headerStartLine + 1], "\t")), 
                                header = tolower(unlist(str_split(metrics[headerStartLine], "\t"))))
    metrics_df <- rbind(metrics_df, unlist(str_split(metrics[headerStartLine + 2], "\t")))
    metrics_df <- rbind(metrics_df, unlist(str_split(metrics[headerStartLine + 3], "\t")))
    sampleQC$Number_of_reads_after_removing_duplicates <- metrics_df[metrics_df$category == "PAIR", "pf_reads"]
  }else{
    print(paste0("File: ", basic_mapping_metrics_result_file, " does not exist!"))
  }
  
  ### flagstat
  flagstat_result_file <- file.path(file.path(config$results_path,"logs"), paste0(sample_basename,"_",config_tools[config_tools$process=="flagstat_filtered_bam","logfile"]))
  if(file.exists(flagstat_result_file)){
    metrics <- readLines(flagstat_result_file)
    metrics <- unlist(str_split(metrics[1]," "))
    sampleQC$Number_of_reads_after_filtering <- as.numeric(metrics[1])
    sampleQC$Prc_passed_filtering_step <- 100 * as.numeric(sampleQC$Number_of_reads_after_filtering)/as.numeric(sampleQC$Number_of_reads_after_removing_duplicates)
  }else{
    print(paste0("File: ", flagstat_result_file, " does not exist!"))
  }
  
  ### on_target_reads
  on_target_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="on_target_reads","temp_results_dirs"], sample_basename)
  on_target_result_file <- paste0(on_target_result_file,"_on_target_reads.txt")
  if(file.exists(on_target_result_file)){
    sampleQC$Number_of_on_target_reads <- as.numeric(yaml::yaml.load_file(on_target_result_file)[['number_of_on_target_reads']])
    sampleQC$Prc_of_on_target_reads <- 100 * sampleQC$Number_of_on_target_read/as.numeric(sampleQC$Number_of_reads_after_removing_duplicates)
  }else{
    print(paste0("File: ", on_target_result_file, " does not exist!"))
  }
  
  ### depth_of_cov
  depth_of_cov_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="depth_of_coverage","temp_results_dirs"], sample_basename)
  depth_of_cov_result_file <- paste0(depth_of_cov_result_file,"_gatk_target_coverage.sample_summary")
  if(file.exists(depth_of_cov_result_file)){
    depth_summary <- readLines(depth_of_cov_result_file)
    depth_summary <- makeDataFrame(data = unlist(str_split(depth_summary[2], "\t")), 
                                   header = tolower(unlist(str_split(depth_summary[1], "\t"))))
    sampleQC$Mean_coverage <- as.numeric(depth_summary[1, "mean"])
  }else{
    print(paste0("File: ", on_target_result_file, " does not exist!"))
  }
  
  ### Calculate conversion efficiency
  methyl_result_file <- file.path(config$results_path, config_tools[config_tools$proces=="methratio","temp_results_dirs"], sample_basename)
  methyl_result_data <- as.data.frame(readMethResult(paste0(methyl_result_file,".methylation_results.", result_format),version=2))
  #head(methyl_result_data)
  if(!is.null(methyl_result_data)){
    control <- methyl_result_data[methyl_result_data$chr == config$ref_control_sequence_name,]
    sampleQC$Number_of_Cs_in_control <- nrow(control)
    conversion_eff <- c(100 * (1-sum(control$numCs)/sum(control$numTs)))
    sampleQC$Conversion_eff <- conversion_eff
    
    #Add to QC Sample report
    methyl_res_panel_df_no_control <- methyl_result_data[methyl_result_data$chr != config$ref_control_sequence_name,]
    sampleQC$Number_of_Cs_in_panel <- nrow(methyl_res_panel_df_no_control)
    
    methyl_res_panel_df_no_control_CpG <- methyl_res_panel_df_no_control[methyl_res_panel_df_no_control$context %in% c('CG'),]
    sampleQC$Number_of_Cs_in_panel_CpG <- nrow(methyl_res_panel_df_no_control_CpG)
    
    methyl_res_panel_df_no_control_nonCpG <- methyl_res_panel_df_no_control[methyl_res_panel_df_no_control$context %in% c('CHG','CHH'),]
    sampleQC$Number_of_Cs_in_panel_non_CpG <- nrow(methyl_res_panel_df_no_control_nonCpG)
    
    sampleQC$Number_of_Cs_in_panel_CpG_cov_min10 <- nrow(filterMethResult(methyl_result_data, config$ref_control_sequence_name, context = c('CG'), min_coverage = 10))
    sampleQC$Number_of_Cs_in_panel_non_CpG_cov_min10 <- nrow(filterMethResult(methyl_result_data, config$ref_control_sequence_name, context = c('CHG','CHH'), min_coverage = 10))

    #methyl_res_panel_df_no_control_CpG_max9 <- methyl_res_panel_df_no_control_CpG[methyl_res_panel_df_no_control_CpG$coverage<10,]
    sampleQC$Number_of_Cs_in_panel_CpG_cov_max9 <- sum(methyl_res_panel_df_no_control_CpG$coverage<10)
    
    sampleQC$Prc_of_Cs_in_panel_CpG_cov_min10 <- 100 * (sampleQC$Number_of_Cs_in_panel_CpG_cov_min10/sampleQC$Number_of_Cs_in_panel_CpG)
    sampleQC$Prc_of_Cs_in_panel_CpG_cov_max9 <- 100 * (sampleQC$Number_of_Cs_in_panel_CpG_cov_max9/sampleQC$Number_of_Cs_in_panel_CpG)
  }
  
  return(sampleQC)
}

####################################################
###### Table of QC common for all present samples ####
createAllQCSummary <- function(config){
  config_tools <- read.csv(file.path(config$tools_path, config$tools_config), stringsAsFactors = FALSE)
  samples <- list.files(file.path(config$results_path,"logs"),full.names = F)
  samples <- samples[endsWith(samples, config_tools[config_tools$tool=="trimmomatic","logfile"])]
  samples <- str_replace_all(samples, pattern = paste0("_",config_tools[config_tools$tool=="trimmomatic","logfile"]),replace = "")
  
  for(s in samples){
    print(paste0("Processing sample: ", s))
          sampleQC <- getSampleQCSummary(s, config)
  }
  return(T)
}

####################################################
###### Table of QC common for all present samples ####
readAllQCSummary <- function(config, save = T){
  
  input_dir <- file.path(config$results_path, "QC_report")
  input_files <- list.files(input_dir, pattern = "\\_QC_summary.yml$", full.names = T)
  
  if(length(input_files) == 0)
    stop(paste0("There is no sample result files to create the report. Directory", input_dir))
  
  qc_summary <- lapply(input_files, function(x) as.data.frame(yaml::read_yaml(x), stringsAsFactors = F))
  qc_summary <- data.frame(rbindlist(qc_summary))
  
  no_numeric_cols <- c("Sample_ID", "processing_time")
  for(curcol in names(qc_summary)){
    if(!curcol %in% no_numeric_cols){
      qc_summary[,curcol] <- as.numeric(qc_summary[,curcol])
    }
  }
  
  if(save){
    fwrite(qc_summary, file = file.path(config$results_path,"QC_report/","SummaryQC.csv"))
  }
  return(qc_summary)
}

####################################################
#### PLOT: Number of sites in CpG context covered by more/less than 10
plotSitesCpG <- function(cov_summary_data, config, min_coverage = 10, pal = brewer.pal(8, "Dark2"), share = F, save = T, fontsize = 10){
  cov <- min_coverage
  if(!cov %in% cov_summary_data$min_coverage){
    warning(paste0('Coverage: ',cov, ' is not present in input cov_summary_data.'))
    return(NULL)
  }

  cov_summary_data <- cov_summary_data[cov_summary_data$context == 'CpG' & (cov_summary_data$min_coverage == cov | cov_summary_data$max_coverage == cov),]
  cov_summary_data$coverage <- paste0("cov>=",cov)
  cov_summary_data$coverage[!is.na(cov_summary_data$max_coverage) & cov_summary_data$max_coverage == cov] <- paste0("cov<",cov)
  cov_summary_data$coverage <- factor(cov_summary_data$coverage, levels = c(paste0("cov<",cov), paste0("cov>=",cov)))
  
  gg <- ggplot(cov_summary_data, aes(fill=coverage, y=value, x=SampleID))
  if(share){
    gg <- gg + geom_bar(stat="identity", position="fill")
  }else{
    gg <- gg + geom_bar(stat="identity")
  }
  
  gg <- gg + theme_minimal() +
    theme_light() + scale_fill_manual(values=pal[1:2]) +
    theme(text = element_text(size = fontsize)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = fontsize)) +
    theme(axis.text.y = element_text(hjust = 1, size = fontsize)) +
    theme(legend.text = element_text(size = fontsize)) + 
    ylab("Number of sites") +
    xlab("Sample ID") +
    ggtitle(paste0("Number of sites in CpG context \ncovered by more/less than ", cov))
  
  if(save){
    ggsave(filename = file.path(config$results_path,"QC_report/",paste0("SummarySitesCovBy",cov,'.',config$plot_format)), plot = gg, device=config$plot_format)
  }
  
  return(gg)
}

####################################################
#### PLOT: Number of sites covered by minimum 10 reads \nseparately for sites in CpG and non-CpG context
plotSitesNonCpG <- function(cov_summary_data, config, min_coverage = 10, pal = brewer.pal(8, "Dark2"), share = F, save = T, fontsize = 10){
  cov <- min_coverage
  if(!cov %in% cov_summary_data$min_coverage){
    warning(paste0('Coverage: ',cov, ' is not present in input cov_summary_data.'))
    return(NULL)
  }
  
  cov_summary_data <- cov_summary_data[cov_summary_data$min_coverage == cov,]
  cov_summary_data$context <- factor(cov_summary_data$context, levels = c("CpG","non-CpG"))
  
  gg <- ggplot(cov_summary_data, aes(fill=context, y=value, x=SampleID))
  if(share){
    gg <- gg + geom_bar(stat="identity", position="fill")
  }else
    gg <- gg + geom_bar(stat="identity")
  
  gg <- gg + theme_minimal() +
    theme_light() + scale_fill_manual(values=pal[3:4]) + 
    theme(text = element_text(size = fontsize)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = fontsize)) +
    theme(axis.text.y = element_text(hjust = 1, size = fontsize)) +
    theme(legend.text = element_text(size = fontsize)) + 
    ylab("Number of sites") +
    xlab("Sample ID") +
    ggtitle(paste0("Number of cytosines covered by at least ",cov," reads \nin the CpG and non-CpG context"))
  if(save){
    ggsave(filename=file.path(config$results_path,"QC_report/",paste0("SummarySitesCovBy",cov,"CpGnonCpG.",config$plot_format)), plot=gg, device=config$plot_format)
  }
  
  return(gg)
}

####################################################
###### Plot Basic MethylKit stats for one sample #######
plotSingleMethStats <- function(single_meth_data){
  par(mfrow=c(1,2))
  getMethylationStats(single_meth_data, plot=T, both.strands=F, labels=F)
  grid()
  getCoverageStats(single_meth_data, plot=T, both.strands=F, labels=F)
  grid()
  par(mfrow=c(1,1))
  
  return(1)
}

####################################################
###### Generate and save Basic MethylKit stats for all samples #######
plotMethStats <- function(meth_data, config, pal = brewer.pal(8, "Dark2"), save = T){
  for(i in 1:length(meth_data)){
    openPlotFile(file.path(config$results_path, 'QC_report', paste0(meth_data[[i]]@sample.id,"_histCpGStats.", config$plot_format)))
    plotSingleMethStats(meth_data[[i]])
    dev.off()
  }
  return(length(meth_data))
}

####################################################
###### PLOT: Generate and save Coverage boxplot for all samples (by default log10 coverage)
plotMethStatsSummary <- function(meth_data, config, log = TRUE, pal = brewer.pal(8, "Dark2"), save = T, fontsize = 10){
  
  sampleCov <- lapply(meth_data, function(x){
    data.frame(coverage = x$coverage, sample_id = x@sample.id)
  })
  
  sampleCov <- rbindlist(sampleCov)
  sampleCov$sample_id <- factor(sampleCov$sample_id)
  
  xlab <- "Read coverage per base"
  outfile_sufix <- ""
  if(log){
    sampleCov$coverage = log10(sampleCov$coverage)
    xlab <- "Read coverage per base [log10]"
    outfile_sufix <- "Log10"
  }
  gg <- ggplot(sampleCov, aes(x=sample_id, y=coverage)) + 
    geom_boxplot(outlier.colour="darkred", outlier.shape=1, outlier.size=0.5, fill = pal[1], color="black") +
    theme_minimal() +
    theme(text = element_text(size = fontsize)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = fontsize)) +
    theme(axis.text.y = element_text(hjust = 1, size = fontsize)) +
    theme(legend.text=element_text(size = fontsize)) + 
    #theme(plot.title = element_text(hjust = 0.5, size = 14)) +
    ggtitle("CpG coverage for the experiment") +
    scale_y_continuous(name = xlab) +
    scale_x_discrete(name = "Sample ID") 
  
  if(save){
    plotfile <- file.path(config$results_path,'QC_report',paste0('SummaryCpGCoverage',outfile_sufix,".",config$plot_format))
    ggsave(plotfile, gg)
  }
  
  return(gg)
}

####################################################
### PLOT: Beta Values Boxplot
#meth_data <- methData
#config <- conf
plotBetaValuesSummary <- function(meth_data, config, sample_size = 100000, pal = brewer.pal(8, "Dark2"), save = T, fontsize = 10){
  
  data_context <- attr(meth_data, 'context')
  plot_title <- paste0("Beta values of ", data_context, " sites across samples")
  head(meth_data)
  percentage_meth <- lapply(meth_data, function(x){x$numCs/x$coverage})
  methyl_levels <- list()
  i = 1
  for (i in 1:length(meth_data)){
    methyl_levels[[i]] <- data.frame(sample.id = meth_data[[i]]@sample.id, 
                                     beta_values = as.numeric(percentage_meth[[i]]),
                                     stringsAsFactors = F)
    if(sample_size > 0 & !is.na(sample_size) & sample_size < nrow(methyl_levels[[i]])){
      methyl_levels[[i]] <- (methyl_levels[[i]])[sample(nrow(methyl_levels[[i]]), sample_size), ]
    }
  }
  methyl_levels <- rbindlist(methyl_levels)
  methyl_levels <- methyl_levels[order(methyl_levels$sample.id),]
  methyl_levels <- data.frame(methyl_levels)
  methyl_levels$sample.id <- factor(methyl_levels$sample.id, levels = sort(unique(methyl_levels$sample.id)))

  gg <- ggplot(methyl_levels, aes(x=sample.id, y=beta_values)) +
    geom_boxplot(outlier.colour="red", outlier.shape=42, outlier.size=3, notch=FALSE, fill=pal[1], color="black") +
    theme_minimal() +
    theme(text = element_text(size = fontsize)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = fontsize)) +
    theme(axis.text.y = element_text(hjust = 1, size = fontsize)) +
    theme(legend.position="none") +  
    ggtitle(plot_title) +
    xlab("Sample ID") + 
    ylab("Beta value")
  
  if(save){
    plotfile <- file.path(config$results_path,'QC_report', paste0('BetaValuesSummary', data_context, '.', config$plot_format))
    ggsave(plotfile, gg)
  }
  
  return(gg)
}
####################################################
### PLOT: Methylation Levels
plotMethLevels <- function(meth_data, config, breaks = c(0,0.1,0.2,0.4,0.6,0.8,0.9,1), pal = brewer.pal(8, "Dark2"), share = F, save = T, fontsize = 10){
  
  percentage_meth <- lapply(meth_data, function(x){x$numCs/x$coverage } )
  methyl_levels <- list()
  for (i in 1:length(meth_data)){
    cnt <- table(cut(percentage_meth[[i]], breaks = breaks, right = F, include.lowest = T))
    methyl_levels[[i]] <- data.frame(sample.id = meth_data[[i]]@sample.id, 
                                     frequency = as.numeric(cnt),
                                     methyl_level = paste0(gsub("]",'',gsub(")",'',gsub(",",'-',gsub("\\[",'',names(cnt))))),"%"),
                                     stringsAsFactors = F)
  }
  
  methyl_levels <- rbindlist(methyl_levels)
  methyl_levels <- methyl_levels[order(methyl_levels$sample.id),]
  methyl_levels <- data.frame(methyl_levels)
  class(methyl_levels$sample.id)
  methyl_levels$sample.id <- factor(methyl_levels$sample.id, levels = sort(unique(methyl_levels$sample.id)))
  methyl_levels$methyl_level <- factor(methyl_levels$methyl_level, levels = rev(sort(unique(methyl_levels$methyl_level))))
  
  gg <- ggplot(methyl_levels, aes(fill=methyl_level, y=frequency, x=sample.id))
  if(share){
    gg <- gg + geom_bar(stat="identity", position="fill")
  }else
    gg <- gg + geom_bar(stat="identity")
  
  gg <- gg + theme_minimal() +
    scale_fill_manual(values = pal) +
    theme(text = element_text(size = fontsize)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = fontsize)) +
    theme(axis.text.y = element_text(hjust = 1, size = fontsize)) +
    theme(legend.text=element_text(size = fontsize)) + 
    ggtitle("Frequency of CpGs according to \nspecific methylation level ranges") +
    xlab("Sample ID") + 
    ylab("Frequency")
  
  if(save){
    plotfile <- file.path(config$results_path,'QC_report',paste0('SummaryMethylationLevel.',config$plot_format))
    ggsave(plotfile, gg)
  }
  
  return(gg)
}

####################################################
###### PLOT: Hyper/Hypo CpG genom & islands annotation for each sample.
#meth_data <- methData
#config <- conf
plotCpGAnnotation <- function(meth_data, hypo_hyper_def = c(0.2,0.8), config, pal = brewer.pal(8, "Dark2"), share = F, save = T, fontsize = 10){
  
  cpg.gene <- NULL
  cpg.island <- NULL
  if(str_trim(config$ref_data_CpGGenomAnnotation) != ""){
    cpg.gene <- readTranscriptFeatures(file.path(config$ref_data_path, config$ref_data_CpGGenomAnnotation))
  }
  if(str_trim(config$ref_data_CpgIslandAnnotation) != ""){
    cpg.island <- readFeatureFlank(file.path(config$ref_data_path, config$ref_data_CpgIslandAnnotation), feature.flank.name=c("CpGi","shores"))
  }
  percentage_meth <- lapply(meth_data, function(x){ x$numCs/x$coverage } )
  
  meth_data_hypo <- list()
  meth_data_hyper <- list()
  
  for (i in 1:(length(meth_data))){
    meth_data_hypo[[meth_data[[i]]@sample.id]] <- meth_data[[i]][percentage_meth[[i]] < hypo_hyper_def[1]]
    meth_data_hyper[[meth_data[[i]]@sample.id]] <- meth_data[[i]][percentage_meth[[i]] >= hypo_hyper_def[2]]
  }
  
  plotList <- list()
  if(!is.null(cpg.gene)){
    plotList[[length(plotList)+1]] <- plotCpGGenomAnnotation(meth_data_hyper, cpg.gene, config, pal, subtitle = "Hypermethylated", share, save, fontsize)
    plotList[[length(plotList)+1]] <- plotCpGGenomAnnotation(meth_data_hypo, cpg.gene, config, pal, subtitle = "Hypomethylated", share, save, fontsize)
  }
  if(!is.null(cpg.island)){
    plotList[[length(plotList)+1]] <- plotCpGIslandsAnnotation(meth_data_hyper, cpg.island, config, pal, subtitle = "Hypermethylated", share, save, fontsize)
    plotList[[length(plotList)+1]] <- plotCpGIslandsAnnotation(meth_data_hypo, cpg.island, config, pal, subtitle = "Hypomethylated", share, save, fontsize)
  }
  
  return(plotList)
}

####################################################
###### PLOT: Hyper/Hypo CpG genom annotation for each sample.
plotCpGGenomAnnotation <- function(meth_data, gene_annot_data, config, pal = brewer.pal(8, "Dark2"), subtitle = "", share = F, save = T, fontsize = 10){
  
  annot_summary <- lapply(meth_data, function(x){ 
    annot <- annotateWithGeneParts(as(x,"GRanges"), gene_annot_data) 
    d <- data.frame(reshape2::melt(t(getTargetAnnotationStats(annot, percentage=TRUE,precedence=F))))
    d$Var1 <- x@sample.id
    names(d) <- c("Sample_Id", "Gene_Part", "Frequency")
    return(d)
  })
  
  annot_summary <- rbindlist(annot_summary)
  annot_summary <- annot_summary[order(annot_summary$Sample_Id),]
  annot_summary$Sample_Id <- factor(annot_summary$Sample_Id, levels = sort(unique(annot_summary$Sample_Id)))
  
  gg <- ggplot(annot_summary, aes(fill=Gene_Part, y=Frequency, x=Sample_Id))
  if(share){
    gg <- gg + geom_bar(stat="identity", position="fill")
  }else{
    gg <- gg + geom_bar(stat="identity")
  }
  gg <- gg + theme_minimal() +
    scale_fill_manual(values = pal) +
    theme(text = element_text(size = fontsize)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = fontsize)) +
    theme(axis.text.y = element_text(hjust = 1, size = fontsize)) +
    theme(legend.text=element_text(size = fontsize)) + 
    ggtitle(paste(subtitle,"CpGs annotated to \ngenomic regions")) +
    xlab("Sample ID") + 
    ylab("Frequency")
  
  if(save){
    ggsave(filename=file.path(config$results_path,"QC_report/",paste0("Summary",subtitle,"CpGGenomAnnotation.",config$plot_format)), plot=gg, device=config$plot_format)
  }
  
  return(gg)
}

####################################################
###### PLOT: Hyper/Hypo CpG islands annotation for each sample.
plotCpGIslandsAnnotation <- function(meth_data, cpg_annot_data, config, pal = brewer.pal(8, "Dark2"), subtitle = "", share = F, save = T, fontsize = 10){
  
  annot_summary <- lapply(meth_data, function(x){ 
    annot <- annotateWithFeatureFlank(as(x,"GRanges"), cpg_annot_data$CpGi, cpg_annot_data$shores, feature.name="CpGi", flank.name="shores")
    d <- data.frame(reshape2::melt(t(getTargetAnnotationStats(annot, percentage=TRUE,precedence=F))))
    d$Var1 <- x@sample.id
    names(d) <- c("Sample_Id", "Region", "Frequency")
    return(d)
  })
  
  annot_summary <- rbindlist(annot_summary)
  annot_summary <- annot_summary[order(annot_summary$Sample_Id),]
  annot_summary$Sample_Id <- factor(annot_summary$Sample_Id, levels = sort(unique(annot_summary$Sample_Id)))
  
  gg <- ggplot(annot_summary, aes(fill=Region, y=Frequency, x=Sample_Id))
  if(share){
    gg <- gg + geom_bar(stat="identity", position="fill")
  }else{
    gg <- gg + geom_bar(stat="identity")
  }
  gg <- gg + theme_minimal() +
    scale_fill_manual(values = pal) +
    theme(text = element_text(size = fontsize)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = fontsize)) +
    theme(axis.text.y = element_text(hjust = 1, size = fontsize)) +
    theme(legend.text=element_text(size = fontsize)) + 
    ggtitle(paste(subtitle,"CpGs annotated to \nCpG-islands")) +
    xlab("Sample ID") + 
    ylab("Frequency")
  
  if(save){
    ggsave(filename=file.path(config$results_path,"QC_report/",paste0("Summary",subtitle,"CpGIslandsAnnotation.",config$plot_format)), plot=gg, device=config$plot_format)
  }
  
  return(gg)
}

###################################################
##### Number of Common CpG Shared Between Samples
plotCntCommonCpG <- function(meth_data, config, pal = brewer.pal(8, "Dark2"), save = T, fontsize = 10){
  
  commonCpg <- lapply(meth_data, function(x){
    data.frame(cpg = paste0(x$chr, "_", x$start))
  })
  commonCpg <- as.data.frame(table(rbindlist(commonCpg)))
  
  gg <- ggplot(commonCpg, aes(x=Freq)) + 
    geom_bar(fill = pal[1]) +
    theme_minimal() +
    scale_x_continuous(breaks = round(seq(min(1), max(length(meth_data), by = 1),1))) +
    theme(text = element_text(size = fontsize)) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = fontsize)) +
    theme(axis.text.y = element_text(hjust = 1, size = fontsize)) +
    theme(legend.text=element_text(size = fontsize)) + 
    ggtitle("Common CpG shared between samples") +
    xlab("Number of samples") + 
    ylab("Number of common CpG")
  
  if(save){
    ggsave(filename=file.path(config$results_path,"QC_report/",paste0("SummaryCntCommonCpG.",config$plot_format)), plot=gg, device=config$plot_format)
  }
  
  return(gg)
}

