# Build Cox model based on TCGA-COAD train dataset
# ----------------------
source("requirements.R")

# Inputpaths
counts_train_inputpath <- "data/data_partitions/counts_train.rds"
counts_test_inputpath <- "data/data_partitions/counts_test.rds"
metadata_train_inputpath <- "data/data_partitions/metadata_train.rds"
metadata_test_inputpath <- "data/data_partitions/metadata_test.rds"
yaccs_inputpath <- "data/yaccs/yaccs_annot.rds"

# Outputpath
outputdir <- "data/cox/tcga-train/"

# Arguments
cvrts <- c(
    "vital_status", "SurvTime",
    "age_at_diagnosis", "gender_male", "Molecular_Subtype_noCIN",
    "Stage_II", "Stage_III", "Stage_IV", "MSI Status_MSI-L/MSS"
    )

# Load data
counts_train <- readRDS(counts_train_inputpath)
counts_test <- readRDS(counts_test_inputpath)
metadata_train <- readRDS(metadata_train_inputpath)
metadata_test <- readRDS(metadata_test_inputpath)
yaccs <- readRDS(yaccs_inputpath)

# Prepare data
train <- counts_train[, yaccs$ENSEMBL]
colnames(train) <- yaccs$SYMBOL
train <- scale(train, scale = FALSE)
metadata_train$age_at_diagnosis <- scale(metadata_train$age_at_diagnosis)

train <- cbind.data.frame(
    metadata_train[, cvrts],
    train
)
rownames(train) <- metadata_train$barcode
train$vital_status <- ifelse(train$vital_status == "Dead", TRUE, FALSE)

# Train model
cox_all <- coxph(Surv(SurvTime, vital_status) ~ ., data = train, x = TRUE)
cox_cvrts <- coxph(Surv(SurvTime, vital_status) ~ ., data = train[,cvrts], x = TRUE)

# Save model
saveRDS(train, file = file.path(outputdir, "train_all.rds"))
saveRDS(cox_all, file = file.path(outputdir, "cox_cvrts_yaccs.rds"))
saveRDS(cox_cvrts, file = file.path(outputdir, "cox_cvrts.rds"))