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
outputdir <- "data/cox/tcga/"

# Arguments
surv_vars <- c("vital_status", "SurvTime")
cvrts <- c(
    "age_at_diagnosis", "gender_male", "Molecular_Subtype_noCIN",
    "Stage_II", "Stage_III", "Stage_IV", "MSI.Status_MSI.L.MSS"
    )

# Load data
counts_train <- readRDS(counts_train_inputpath)
counts_test <- readRDS(counts_test_inputpath)
metadata_train <- readRDS(metadata_train_inputpath)
metadata_test <- readRDS(metadata_test_inputpath)
yaccs <- readRDS(yaccs_inputpath)
yaccs$SYMBOL <- make.names(yaccs$SYMBOL)

# Prepare train data
train <- counts_train[, yaccs$ENSEMBL]
colnames(train) <- yaccs$SYMBOL
train <- scale(train, scale = FALSE)

train <- cbind.data.frame(
    metadata_train[, c(surv_vars, cvrts)],
    train
)
rownames(train) <- metadata_train$barcode
train$vital_status <- ifelse(train$vital_status == "Dead", TRUE, FALSE)

# Train model
cox_all <- coxph(Surv(SurvTime, vital_status) ~ ., data = train, x = TRUE)
cox_cvrts <- coxph(Surv(SurvTime, vital_status) ~ ., data = train[,c(surv_vars, cvrts)], x = TRUE)
cox_yaccs <- coxph(Surv(SurvTime, vital_status) ~ ., data = train %>% dplyr::select(-all_of(cvrts)), x = TRUE)

# Prepare test data
test <- counts_test[, yaccs$ENSEMBL]
colnames(test) <- yaccs$SYMBOL
test <- scale(test, scale = FALSE)

test <- cbind.data.frame(
    metadata_test[, c(surv_vars, cvrts)],
    test
)
rownames(test) <- metadata_test$barcode
test$vital_status <- ifelse(test$vital_status == "Dead", TRUE, FALSE)
names(test) <- make.names(names(test))


# Performances in train ----------------
pred_train_all <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_all, train), train)
pred_train_cvrts <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_cvrts, train[, c(surv_vars, cvrts)]), train)
pred_train_yaccs <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_yaccs, train[, c(surv_vars, yaccs$SYMBOL)]), train)

cindex_train <- data.frame(
    C = c(pred_train_all$concordance, pred_train_cvrts$concordance, pred_train_yaccs$concordance),
    ci_l = c(
        pred_train_all$concordance - 1.96 * pred_train_all$std.err,
        pred_train_cvrts$concordance - 1.96 * pred_train_cvrts$std.err,
        pred_train_yaccs$concordance - 1.96 * pred_train_yaccs$std.err),
    ci_u = c(
        pred_train_all$concordance + 1.96 * pred_train_all$std.err,
        pred_train_cvrts$concordance + 1.96 * pred_train_cvrts$std.err,
        pred_train_yaccs$concordance + 1.96 * pred_train_yaccs$std.err),
    Subset = "Train"
)


# Testing (TCGA-COAD test) ----------------
# Three scenarios:
#   Model with yaccs + cvrts
#   Model with only cvrts
#   Model with only yaccs

pred_test_all <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_all, test), test)
pred_test_cvrts <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_cvrts, test[, c(surv_vars, cvrts)]), test)
pred_test_yaccs <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_yaccs, test[, c(surv_vars, yaccs$SYMBOL)]), test)

cindex_test <- data.frame(
    C = c(pred_test_all$concordance, pred_test_cvrts$concordance, pred_test_yaccs$concordance),
    ci_l = c(
        pred_test_all$concordance - 1.96 * pred_test_all$std.err,
        pred_test_cvrts$concordance - 1.96 * pred_test_cvrts$std.err,
        pred_test_yaccs$concordance - 1.96 * pred_test_yaccs$std.err),
    ci_u = c(
        pred_test_all$concordance + 1.96 * pred_test_all$std.err,
        pred_test_cvrts$concordance + 1.96 * pred_test_cvrts$std.err,
        pred_test_yaccs$concordance + 1.96 * pred_test_yaccs$std.err),
    Subset = "Test"
)

# Results of C-Index in TCGA (train and test)
cindex_tcga <- rbind.data.frame(cindex_train, cindex_test)


# Save model
saveRDS(train, file = file.path(outputdir, "train_all.rds"))
saveRDS(test, file = file.path(outputdir, "test_all.rds"))

saveRDS(cox_all, file = file.path(outputdir, "cox_cvrts_yaccs.rds"))
saveRDS(cox_cvrts, file = file.path(outputdir, "cox_cvrts.rds"))
saveRDS(cox_yaccs, file = file.path(outputdir, "cox_yaccs.rds"))

saveRDS(cindex_tcga, file = file.path(outputdir, "cindex_tcga.rds"))
