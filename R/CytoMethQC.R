source("./R/main.R")
conf <- readConfig()

# create all final summary plots
qcsumm <- getQCSummary(conf, save = T)
plotSitesCovBy10(qcsumm, conf, share = F, save = T)
plotSitesCovBy10CpGnonCpG(qcsumm, conf, share = F, save = T)

methData <- readMethData(conf)
plotSingleMethStats(methData[[1]])
plotSingleMethStats(methData[[3]])

plotMethStats(methData, conf, save = T)
plotMethStatsSummary(methData, conf, save = T)
plotMethLevels(methData, conf, share = T, save = T)

plotCpGAnnotation(methData, hypo_hyper_def = c(20,80), conf, share = T, save = T)

plotCntCommonCpG(methData, conf, save = T)
