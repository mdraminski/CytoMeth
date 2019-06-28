####################################################
###### Table of QC common for all present samples ####
getQCSummary <- function(config, save = T){
  
  input_dir <- file.path(config$results_path, "QC_report")
  input_files <- list.files(input_dir, pattern = "\\_QC_summary.yml$", full.names = T)
  
  if(length(input_files) == 0)
    stop(paste0("There is no sample result files to create the report. Directory", input_dir))
  
  qc_summary <- lapply(input_files, function(x) as.data.frame(yaml::read_yaml(x)))
  qc_summary <- data.frame(rbindlist(qc_summary))
  
  if(save){
    fwrite(qc_summary, file = file.path(config$results_path,"QC_report/","SummaryQC.csv"))
  }
  return(qc_summary)
}

####################################################
#### PLOT: Number of sites in CpG context covered by more/less than 10
plotSitesCovBy10 <- function(qc_summary, config, pal = brewer.pal(8, "Dark2"), share = F, save = T){
  
  qc_summary_coverage <- rbind(data.frame(SampleID = qc_summary$Sample_ID, Number = qc_summary$Number_of_Cs_in_panel_CpG_cov_min10, Coverage = "cov>=10", stringsAsFactors = F),
                               data.frame(SampleID = qc_summary$Sample_ID, Number = qc_summary$Number_of_Cs_in_panel_CpG_cov_max9, Coverage = "cov<10", stringsAsFactors = F))
  qc_summary_coverage$Coverage <- factor(qc_summary_coverage$Coverage, levels = c("cov<10","cov>=10"))
  
  gg <- ggplot(qc_summary_coverage, aes(fill=Coverage, y=Number, x=SampleID))
  if(share){
    gg <- gg + geom_bar(stat="identity", position="fill")
  }else
    gg <- gg + geom_bar(stat="identity")
  
  gg <- gg + theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    theme(legend.text=element_text(size=15)) + 
    theme(text = element_text(size=15)) +
    theme_light() + scale_fill_manual(values=pal[1:2]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Number of sites") +
    xlab("Sample ID") +
    ggtitle("Number of sites in CpG context covered by more/less than 10.")
  
  if(save){
    ggsave(filename = file.path(config$results_path,"QC_report/",paste0("SummarySitesCovBy10.",config$plot_format)), plot = gg, device = config$plot_format)
  }
  
  return(gg)
}

####################################################
#### PLOT: Number of sites covered by minimum 10 reads \nseparately for sites in CpG and non-CpG context
plotSitesCovBy10CpGnonCpG <- function(qc_summary, config, pal = brewer.pal(8, "Dark2"), share = F, save = T){
  
  qc_summary_cpg <- rbind(data.frame(SampleID = qc_summary$Sample_ID, Number = qc_summary$Number_of_Cs_in_panel_CpG_cov_min10, Context = "CpG", stringsAsFactors = F),
                          data.frame(SampleID = qc_summary$Sample_ID, Number = qc_summary$Number_of_Cs_in_panel_non_CpG_cov_min10, Context = "non-CpG", stringsAsFactors = F))
  
  gg <- ggplot(qc_summary_cpg, aes(fill=Context, y=Number, x=SampleID))
  if(share){
    gg <- gg + geom_bar(stat="identity", position="fill")
  }else
    gg <- gg + geom_bar(stat="identity")
  
  gg <- gg + theme_minimal() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(legend.text=element_text(size=15)) + 
    theme(text = element_text(size=15)) +
    theme_light() + scale_fill_manual(values=pal[3:4]) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Number of sites") +
    xlab("Sample ID") +
    ggtitle("Number of sites covered by minimum 10 reads \nseparately for sites in CpG and non-CpG context.")
  if(save){
    ggsave(filename=file.path(config$results_path,"QC_report/",paste0("SummarySitesCovBy10CpGnonCpG.",config$plot_format)), plot=gg, device = config$plot_format)
  }
  
  return(gg)
}

##############################################
######## Read data using MethylKit ############
readMethData <- function(config){
  
  result_dir <- file.path(config$results_path, "methyl_results")
  result_files <- list.files(result_dir, pattern = "\\.CpG_min10", full.names = T)
  sample_id <- lapply(strsplit(basename(result_files), '\\.'), function(x){x[1]})
  
  methData <- methRead(as.list(result_files),
                       sample.id = sample_id,
                       treatment = c(rep(0, length(sample_id))), # 0/1 control/test samples
                       assembly = "hg38",
                       header = TRUE,
                       context = "CpG",
                       resolution = "base",
                       pipeline = list(fraction=TRUE, chr.col=1, start.col=2, end.col=3, coverage.col=7, strand.col=4, freqC.col=6)
  )
  return(methData)
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
    openPlotFile(file.path(config$results_path,'QC_report',paste0(meth_data[[i]]@sample.id,"_histCpGStats.",config$plot_format)))
    plotSingleMethStats(meth_data[[i]])
    dev.off()
  }
  return(length(meth_data))
}

####################################################
###### Generate and save log10 Coverage boxplot for all samples
plotMethStatsSummary <- function(meth_data, config, pal = brewer.pal(8, "Dark2"), save = T){
  
  sampleCov <- lapply(meth_data, function(x){
    data.frame(coverage = log10(x$coverage), sample_id = x@sample.id)
  })
  sampleCov <- rbindlist(sampleCov)
  sampleCov$sample_id <- factor(sampleCov$sample_id)
  
  gg <- ggplot(sampleCov, aes(x=sample_id, y=coverage)) + 
    geom_boxplot(outlier.colour="darkred", outlier.shape=1, outlier.size=0.5, fill = pal[1], color="black") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("CpG Coverage for the experiment") +
    scale_y_continuous(name = "log 10 of read Coverage per base") +
    scale_x_discrete(name = "Sample Id") 
  
  if(save){
    plotfile <- file.path(config$results_path,'QC_report',paste0('SummaryCpGCoverage.',config$plot_format))
    ggsave(plotfile, gg)
  }
  
  return(gg)
}

########## create % info ###########
### Methylation Levels
plotMethLevels <- function(meth_data, config, breaks = c(0,10,20,40,60,80,90,100), pal = brewer.pal(8, "Dark2"), share = F, save = T){
  
  percentage_meth <- lapply(meth_data, function(x){x$numCs/x$coverage*100 } )
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
    geom_bar(stat="identity", position="fill") + 
    scale_fill_manual(values=pal) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.text=element_text(size=15)) + 
    theme(text = element_text(size=15)) + 
    ggtitle("Methylation Levels") +
    xlab("Sample Id") + 
    ylab("Frequency")
  
  if(save){
    plotfile <- file.path(config$results_path,'QC_report',paste0('SummaryMethylationLevel.',config$plot_format))
    ggsave(plotfile, gg)
  }
  
  return(gg)
}
