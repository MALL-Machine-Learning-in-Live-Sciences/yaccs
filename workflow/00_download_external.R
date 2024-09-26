# Download external cohorts
# ==========
source("requirements.R")
# setwd('~/projects/SurvivalTCGA-COAD/data/Validation/validation_datasets/all-genes/')

# Outputpaths
outputdir <- "extdata/geo/"

# Arguments
geoID <- c("GSE17536", "GSE17537", "GSE29621", "GSE39582")
platform <- "GPL570"
chip <- "hgu133plus2"

# Downloading cohorts
for (i in seq_along(geoID)) {
  gse <- getGEO(geoID[i])[[1]]
  saveRDS(gse, file = paste0(outputdir, geoID[i], ".rds"))
}

# Select sondas
# ==========
annot <- GEOquery::getGEO(platform)@dataTable@table[, 1:12]
probes <- annot$ID
js <- jscores(chip = chip, probeset = probes)

# Add specific annotation genes (our signature and mda)
js[which(rownames(js) == '206400_at'),] <- c(11, '3963', NA, NA, NA, NA, NA, 'LGALS7')
js[which(rownames(js) == '238805_at'),] <- c(11, '91894', NA, NA, NA, NA, NA, 'C11orf52')
js[which(rownames(js) == '223628_at'),] <- c(11, '645426', NA, NA, NA, NA, NA, 'TMEM191C')

js[which(rownames(js) == '204044_at'),] <- c(11, '23475', NA, NA, NA, NA, NA, 'QPRT')
js[which(rownames(js) == '206286_s_at'),] <- c(11, '6997', NA, NA, NA, NA, NA, 'TDGF1') 
js[which(rownames(js) == '209424_s_at'),] <- c(11, '23600', NA, NA, NA, NA, NA, 'AMACR')

# delete NA
js <- js[-which(is.na(js$EntrezID)), ]

# order by gene (entrezID)
js <- js[order(js$EntrezID), ]
js <- js[ !duplicated(js$EntrezID), ]
js$probeID <- rownames(js)

# save annotation jetscore data
saveRDS(js, file = file.path(outputdir, "gpl570-annot-jetscore.rds"))

# Load cohorts
# ==========
# files = list.files(pattern = '^GSE')
# for (i in seq_along(files)) {
#   xx = readRDS(files[i])
#   counts = xx@assayData$exprs
#   counts = counts[js$probeID, ]
  
#   saveRDS(list(genes = js$EntrezID, data = counts), file = paste0('prepro-', files[i]))
# }
