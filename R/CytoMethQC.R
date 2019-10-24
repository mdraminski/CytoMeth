source("./R/main.R")
conf <- readConfig()

# go over all samples and recreate all SAMPLE_QC_summary.yml files
# this step is unnecessary if all calculations succeed
# createAllQCSummary(conf)
  
# read all SAMPLE_QC_summary.yml files as one data.frame and save it to csv file
save <- T
qcsumm <- readAllQCSummary(conf, save = save)
# create all final summary plots
gg <- plotSitesCovBy10(qcsumm, conf, share = F, save = save)
gg <- plotSitesCovBy10CpGnonCpG(qcsumm, conf, share = F, save = save)

# Context: non-CpG
methData <- readMethData(conf, context = c("CHG","CHH"), context_label = 'non-CpG', min_coverage = 10)
gg <- plotBetaValuesSummary(methData, conf, sample_size = 100000, save = save)

# Context: CpG
methData <- readMethData(conf, context = c("CG"), context_label = 'CpG', min_coverage = 10)
#gg <- plotSingleMethStats(methData[[1]])
gg <- plotMethStats(methData, conf, save = save)
gg <- plotBetaValuesSummary(methData, conf, sample_size = 100000, save = save)
gg <- plotMethStatsSummary(methData, conf, log = T, save = save)
gg <- plotMethStatsSummary(methData, conf, log = F, save = save)
gg <- plotMethLevels(methData, conf, share = T, save = save)
gg <- plotCpGAnnotation(methData, hypo_hyper_def = c(0.2,0.8), conf, share = T, save = save)
gg <- plotCntCommonCpG(methData, conf, save = save)

