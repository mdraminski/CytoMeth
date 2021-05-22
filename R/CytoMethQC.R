source("./R/main.R")
conf <- readConfig()

# go over all samples and recreate all SAMPLE_QC_summary.yml files
# this step is unnecessary if all calculations succeed
# createAllQCSummary(conf)
  
# read all SAMPLE_QC_summary.yml files as one data.frame and save it to csv file
save <- T
qcsumm <- readAllQCSummary(conf, save = save)
print(qcsumm)

#i=1
for(i in 1:length(conf$meth_tool)){
  meth_tool <- conf$meth_tool[i]
  cat(paste0('##### Creating plots for meth_tool: ', meth_tool, ' #####\n'))
  
  # create all final summary plots
  covSummary <- getCovSummary(conf, meth_tool, min_coverage = c(7,8,9,10,11,12,13), result_format = c('rds'))
  gg <- plotSitesCpG(covSummary, conf, meth_tool, min_coverage = 10, share = F, save = save)
  gg <- plotSitesNonCpG(covSummary, conf, meth_tool, min_coverage = 10, share = F, save = save)
  
  # Context: non-CpG
  methData <- readMethData(conf, meth_tool, context = c("CHG","CHH"), context_label = 'non-CpG', min_coverage = 10, result_format = c('rds'))
  gg <- plotBetaValuesSummary(methData, conf, meth_tool, sample_size = 100000, save = save)
  
  # Context: CpG
  methData <- readMethData(conf, meth_tool, context = c("CG"), context_label = 'CpG', min_coverage = 10)
  #gg <- plotSingleMethStats(methData[[1]])
  gg <- plotMethStats(methData, conf, meth_tool, save = save)
  gg <- plotBetaValuesSummary(methData, conf, meth_tool, sample_size = 100000, save = save)
  gg <- plotMethStatsSummary(methData, conf, meth_tool, log = T, save = save)
  gg <- plotMethStatsSummary(methData, conf, meth_tool, log = F, save = save)
  gg <- plotMethLevels(methData, conf, meth_tool, share = T, save = save)
  gg <- plotCntCommonCpG(methData, conf, meth_tool, save = save)
  gg <- plotCpGAnnotation(methData, conf, meth_tool, hypo_hyper_def = c(0.2,0.8), share = T, save = save)
}
