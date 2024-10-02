# Download external cohorts
# ==========
source("requirements.R")

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