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
    ggsave(filename = file.path(config$results_path,"QC_report/",paste0("SummarySitesCovBy10.",config$plot_format)), plot = gg, device=config$plot_format)
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
    ggsave(filename=file.path(config$results_path,"QC_report/",paste0("SummarySitesCovBy10CpGnonCpG.",config$plot_format)), plot=gg, device=config$plot_format)
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
###### PLOT: Generate and save log10 Coverage boxplot for all samples
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

####################################################
### PLOT: Beta Values Boxplot
plotBetaValuesSummary <- function(meth_data, config, pal = brewer.pal(8, "Dark2"), save = T){
  
  percentage_meth <- lapply(meth_data, function(x){x$numCs/x$coverage } )
  methyl_levels <- list()
  for (i in 1:length(meth_data)){
    methyl_levels[[i]] <- data.frame(sample.id = meth_data[[i]]@sample.id, 
                                     beta_values = as.numeric(percentage_meth[[i]]),
                                     stringsAsFactors = F)
  }
  methyl_levels <- rbindlist(methyl_levels)
  methyl_levels <- methyl_levels[order(methyl_levels$sample.id),]
  methyl_levels <- data.frame(methyl_levels)
  methyl_levels$sample.id <- factor(methyl_levels$sample.id, levels = sort(unique(methyl_levels$sample.id)))
  
  gg <- ggplot(methyl_levels, aes(x=sample.id, y=beta_values)) +
    geom_boxplot(outlier.colour="red", outlier.shape=42, outlier.size=3, notch=FALSE, fill=pal[1], color="black") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(text = element_text(size=15)) + 
    theme(legend.position="none") +  
    ggtitle("Methylation over samples ") +
    xlab("Sample Id") + 
    ylab("Beta Value")
  
  if(save){
    plotfile <- file.path(config$results_path,'QC_report',paste0('BetaValuesSummary.',config$plot_format))
    ggsave(plotfile, gg)
  }
  
  return(gg)
}
####################################################
### PLOT: Methylation Levels
plotMethLevels <- function(meth_data, config, breaks = c(0,0.1,0.2,0.4,0.6,0.8,0.9,1), pal = brewer.pal(8, "Dark2"), share = F, save = T){
  
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

####################################################
###### PLOT: Hyper/Hypo CpG genom & islands annotation for each sample.
#meth_data <- methData
#config <- conf
plotCpGAnnotation <- function(meth_data, hypo_hyper_def = c(0.2,0.8), config, pal = brewer.pal(8, "Dark2"), share = F, save = T){
  
  cpg.gene <- readTranscriptFeatures(file.path(config$ref_data_path, config$ref_data_CpGGenomAnnotation))
  cpg.island <- readFeatureFlank(file.path(config$ref_data_path, config$ref_data_CpgIslandAnnotation), feature.flank.name=c("CpGi","shores"))
  percentage_meth <- lapply(meth_data, function(x){ x$numCs/x$coverage } )
  
  meth_data_hypo <- list()
  meth_data_hyper <- list()
  
  for (i in 1:(length(meth_data))){
    meth_data_hypo[[meth_data[[i]]@sample.id]] <- meth_data[[i]][percentage_meth[[i]] < hypo_hyper_def[1]]
    meth_data_hyper[[meth_data[[i]]@sample.id]] <- meth_data[[i]][percentage_meth[[i]] >= hypo_hyper_def[2]]
  }
  
  plotList <- list(
    plotCpGGenomAnnotation(meth_data_hyper, cpg.gene, config, pal, subtitle = "Hypermethylated", share, save),
    plotCpGGenomAnnotation(meth_data_hypo, cpg.gene, config, pal, subtitle = "Hypomethylated", share, save),
    plotCpGIslandsAnnotation(meth_data_hyper, cpg.island, config, pal, subtitle = "Hypermethylated", share, save),
    plotCpGIslandsAnnotation(meth_data_hypo, cpg.island, config, pal, subtitle = "Hypomethylated", share, save))
  return(plotList)
}

####################################################
###### PLOT: Hyper/Hypo CpG genom annotation for each sample.
plotCpGGenomAnnotation <- function(meth_data, gene_annot_data, config, pal = brewer.pal(8, "Dark2"), subtitle = "", share = F, save = T){
  
  annot_summary <- lapply(meth_data, function(x){ 
    annot <- annotateWithGeneParts(as(x,"GRanges"), gene_annot_data) 
    d <- data.frame(melt(t(getTargetAnnotationStats(annot, percentage=TRUE,precedence=F))))
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.text=element_text(size=15)) + 
    theme(text = element_text(size=15)) +
    ggtitle(paste(subtitle,"CpG Genom Annotation")) +
    xlab("Sample Id") + 
    ylab("Frequency")
  
  if(save){
    ggsave(filename=file.path(config$results_path,"QC_report/",paste0("Summary",subtitle,"CpGGenomAnnotation.",config$plot_format)), plot=gg, device=config$plot_format)
  }
  
  return(gg)
}

####################################################
###### PLOT: Hyper/Hypo CpG islands annotation for each sample.
plotCpGIslandsAnnotation <- function(meth_data, cpg_annot_data, config, pal = brewer.pal(8, "Dark2"), subtitle = "", share = F, save = T){
  
  annot_summary <- lapply(meth_data, function(x){ 
    annot <- annotateWithFeatureFlank(as(x,"GRanges"), cpg_annot_data$CpGi, cpg_annot_data$shores, feature.name="CpGi", flank.name="shores")
    d <- data.frame(melt(t(getTargetAnnotationStats(annot, percentage=TRUE,precedence=F))))
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.text=element_text(size=15)) + 
    theme(text = element_text(size=15)) +
    ggtitle(paste(subtitle,"CpG Islands Annotation")) +
    xlab("Sample Id") + 
    ylab("Frequency")
  
  if(save){
    ggsave(filename=file.path(config$results_path,"QC_report/",paste0("Summary",subtitle,"CpGIslandsAnnotation.",config$plot_format)), plot=gg, device=config$plot_format)
  }
  
  return(gg)
}

###################################################
##### Number of Common CpG Shared Between Samples
plotCntCommonCpG <- function(meth_data, config, pal = brewer.pal(8, "Dark2"), save = T){
  
  commonCpg <- lapply(meth_data, function(x){
    data.frame(cpg = paste0(x$chr, "_", x$start))
  })
  commonCpg <- as.data.frame(table(rbindlist(commonCpg)))
  
  gg <- ggplot(commonCpg, aes(x=Freq)) + 
    geom_bar(fill = pal[1]) +
    theme_minimal() +
    scale_x_continuous(breaks = round(seq(min(1), max(length(meth_data), by = 1),1))) +
    theme(legend.text=element_text(size=15)) + 
    theme(text = element_text(size=15)) + 
    ggtitle("Common CpG Shared Between Samples") +
    xlab("Number of samples") + 
    ylab("Number of common CpG")
  
  if(save){
    ggsave(filename=file.path(config$results_path,"QC_report/",paste0("SummaryCntCommonCpG.",config$plot_format)), plot=gg, device=config$plot_format)
  }
  
  return(gg)
}

