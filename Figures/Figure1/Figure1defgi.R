# Plotting the different SigCheck´s experiments by each signature
# ------------------
require(SigCheck)

inputdir <- "data/sigCheck/"

# YACCS
# ===
random <- readRDS(file.path(inputdir, 'sigCheckRandom_ours.rds'))
features <- readRDS(file.path(inputdir, 'sigCheckPermutedFeatures_ours.rds'))
survival <- readRDS(file.path(inputdir, 'sigCheckPermutedSurvival_ours.rds'))
msigdb <- readRDS(file.path(inputdir, 'scKnownMSigDB_ours.rds'))
literature <- readRDS(file.path(inputdir, 'scKnownLiterature_ours.rds'))

sigCheckPlot(random, title = 'Random signatures', nolegend = F)
sigCheckPlot(features, title = 'Permuted features', nolegend = T)
sigCheckPlot(survival, title = 'Permuted survival', nolegend = T)
sigCheckPlot(msigdb, title = 'MSigDB', nolegend = T)
sigCheckPlot(literature, title = 'Literature signatures', nolegend = T)



# ColoPrint
# ===
random <- readRDS(file.path(inputdir,'sigCheckRandom_coloPrint.rds'))
features <- readRDS(file.path(inputdir,'sigCheckPermutedFeatures_coloPrint.rds'))
survival <- readRDS(file.path(inputdir,'sigCheckPermutedSurvival_coloPrint.rds'))

sigCheckPlot(random, title = 'Random signatures', nolegend = F)
sigCheckPlot(features, title = 'Permuted features', nolegend = T)
sigCheckPlot(survival, title = 'Permuted survival', nolegend = T)



# ColoGuideEx
# ===
random = readRDS(file.path(inputdir,'sigCheckRandom_coloGuideEx.rds'))
features = readRDS(file.path(inputdir,'sigCheckPermutedFeatures_coloGuideEx.rds'))
survival = readRDS(file.path(inputdir,'sigCheckPermutedSurvival_coloGuideEx.rds'))

sigCheckPlot(random, title = 'Random signatures', nolegend = F)
sigCheckPlot(features, title = 'Permuted features', nolegend = T)
sigCheckPlot(survival, title = 'Permuted survival', nolegend = T)



# ColoGuidePro
# ===
random = readRDS(file.path(inputdir,'sigCheckRandom_coloGuidePro.rds'))
features = readRDS(file.path(inputdir,'sigCheckPermutedFeatures_coloGuidePro.rds'))
survival = readRDS(file.path(inputdir,'sigCheckPermutedSurvival_coloGuidePro.rds'))

sigCheckPlot(random, title = 'Random signatures', nolegend = F)
sigCheckPlot(features, title = 'Permuted features', nolegend = T)
sigCheckPlot(survival, title = 'Permuted survival', nolegend = T)



# mda114
# ===
random = readRDS(file.path(inputdir,'sigCheckRandom_mda114.rds'))
features = readRDS(file.path(inputdir,'sigCheckPermutedFeatures_mda114.rds'))
survival = readRDS(file.path(inputdir,'sigCheckPermutedSurvival_mda114.rds'))

sigCheckPlot(random, title = 'Random signatures', nolegend = F)
sigCheckPlot(features, title = 'Permuted features', nolegend = T)
sigCheckPlot(survival, title = 'Permuted survival', nolegend = T)



# zhang
# ===
random = readRDS(file.path(inputdir,'sigCheckRandom_zhang.rds'))
features = readRDS(file.path(inputdir,'sigCheckPermutedFeatures_zhang.rds'))
survival = readRDS(file.path(inputdir,'sigCheckPermutedSurvival_zhang.rds'))

sigCheckPlot(random, title = 'Random signatures', nolegend = F)
sigCheckPlot(features, title = 'Permuted features', nolegend = T)
sigCheckPlot(survival, title = 'Permuted survival', nolegend = T)

