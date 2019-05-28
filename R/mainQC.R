#####################  
#' CytoMethQCReport
#' 
#' Function runs CytoMeth Quality Control processing.
#' 
#' @param config - configuration CytoMeth file 
#' @param qc_colors
#' 
#' @return TRUE if succeed
#'
#' @examples
#' conf <- readConfig()
#' CytoMethQCReport(conf)
#' 
#' @import yaml data.table ggplot2
#' @export
#' 
#config <- conf
CytoMethQCReport <- function(config, qc_colors = c("azure1","darkolivegreen3","antiquewhite2","firebrick")){
  
  file_type <- config$plot_format
  ###### Table of QC common for all present samples ####
  input_dir <- file.path(config$results_path, "QC_report")
  input_files <- list.files(input_dir, pattern = "\\_QC_summary.yml$", full.names = T)
  
  if(length(input_files) == 0)
    stop(paste0("There is no sample result files to create the report. Directory", input_dir))
  
  qc_summary <- lapply(input_files, function(x) as.data.frame(yaml::read_yaml(x)))
  qc_summary <- data.frame(rbindlist(qc_summary))
  
  fwrite(qc_summary, file = file.path(config$results_path,"QC_report/","QC_report.csv"))
  
  ######## barplot ggplot of covarage ###########
  qc_summary_coverage <- rbind(data.frame(SampleID = qc_summary$Sample_ID, Frequency = qc_summary$Number_of_Cs_in_panel_CpG_cov_min10, Coverage = "cov>=10", stringsAsFactors = F),
                               data.frame(SampleID = qc_summary$Sample_ID, Frequency = qc_summary$Number_of_Cs_in_panel_CpG_cov_max9, Coverage = "cov<10", stringsAsFactors = F))
  qc_summary_coverage$Coverage <- factor(qc_summary_coverage$Coverage, levels = c("cov<10","cov>=10"))
  
  gg_coverage <- ggplot(qc_summary_coverage, aes(fill=Coverage, y=Frequency, x=SampleID)) + 
    geom_bar(stat="identity") + 
    #geom_bar(stat="identity", position="fill") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    theme(legend.text=element_text(size=15)) + 
    theme(text = element_text(size=15)) +
    theme_light() + scale_fill_manual(values=qc_colors[1:2])
  #gg_coverage
  ggsave(filename = file.path(config$results_path,"QC_report/",paste0("Coverage_plot.",file_type)), plot = gg_coverage, device = file_type)
  
  ######## barplot ggplot of CG vs not-CG ###########
  qc_summary_cpg <- rbind(data.frame(SampleID = qc_summary$Sample_ID, Frequency = qc_summary$Number_of_Cs_in_panel_CpG_cov_min10, Context = "CpG", stringsAsFactors = F),
                          data.frame(SampleID = qc_summary$Sample_ID, Frequency = qc_summary$Number_of_Cs_in_panel_non_CpG_cov_min10, Context = "non-CpG", stringsAsFactors = F))
  
  gg_cpg <- ggplot(qc_summary_cpg, aes(fill=Context, y=Frequency, x=SampleID)) + 
    geom_bar(stat="identity") + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(legend.text=element_text(size=15)) + 
    theme(text = element_text(size=15)) +
    theme_light() + scale_fill_manual(values=qc_colors[3:4])
  #gg_cpg
  ggsave(filename=file.path(config$results_path,"QC_report/",paste0("CpG_nonCpG_plot.",file_type)), plot=gg_cpg, device = file_type)
  
  return(T)
}

