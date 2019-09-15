source("./R/main.R")
conf <- readConfig()

# go over all samples and recreate all SAMPLE_QC_summary.yml files
# this step is unnecessary if all calculations succeed
# createAllQCSummary(conf)
  
# read all SAMPLE_QC_summary.yml files as one data.frame and save it to csv file
qcsumm <- readAllQCSummary(conf, save = T)
# create all final summary plots
gg <- plotSitesCovBy10(qcsumm, conf, share = F, save = T)
gg <- plotSitesCovBy10CpGnonCpG(qcsumm, conf, share = F, save = T)

# Context: non-CpG
methData <- readMethData(conf, context = 'non-CpG')
gg <- plotBetaValuesSummary(methData, conf, sample_size = 100000, save = T)

# Context: CpG
methData <- readMethData(conf, context = 'CpG')
#gg <- plotSingleMethStats(methData[[1]])
gg <- plotMethStats(methData, conf, save = T)
gg <- plotBetaValuesSummary(methData, conf, sample_size = 100000, save = T)
gg <- plotMethStatsSummary(methData, conf, log = T, save = T)
gg <- plotMethStatsSummary(methData, conf, log = F, save = T)
gg <- plotMethLevels(methData, conf, share = T, save = T)
gg <- plotCpGAnnotation(methData, hypo_hyper_def = c(0.2,0.8), conf, share = T, save = T)
gg <- plotCntCommonCpG(methData, conf, save = T)
