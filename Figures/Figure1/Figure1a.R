# Code to reproduce Figure 1a)
# Forest plot of Hazard ratio: multivariate model
# --------------------------
source("requirements.R")

# Inputpath
cox_inputpath <- "data/cox/tcga-train/cox_cvrts_yaccs.rds"
train_inputpath <- "data/cox/tcga-train/train_all.rds"

# Load data
cox <- readRDS(cox_inputpath)
train <- readRDS(train_inputpath)

# Plot
ggforest(
    model = cox,
    data = train,
    cpositions = c(0.01, NA, NA),
    fontsize = 0.8
    )
