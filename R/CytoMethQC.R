source("./R/main.R")
conf <- readConfig()

# create all final summary plots
qcsumm <- getQCSummary(conf, save = T)
gg <- plotSitesCovBy10(qcsumm, conf, share = F, save = T)
gg <- plotSitesCovBy10CpGnonCpG(qcsumm, conf, share = F, save = T)

methData <- readMethData(conf)
#gg <- plotSingleMethStats(methData[[1]])
gg <- plotMethStats(methData, conf, save = T)
gg <- plotBetaValuesSummary(methData, conf, save = T)
gg <- plotMethStatsSummary(methData, conf, save = T)
gg <- plotMethLevels(methData, conf, share = T, save = T)
gg <- plotCpGAnnotation(methData, hypo_hyper_def = c(0.2,0.8), conf, share = T, save = T)
gg <- plotCntCommonCpG(methData, conf, save = T)
