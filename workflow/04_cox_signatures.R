# Build Cox model based on TCGA-COAD train dataset
# ----------------------
source("requirements.R")

# Inputpaths
counts_train_inputpath <- "data/data_partitions/counts_train.rds"
counts_test_inputpath <- "data/data_partitions/counts_test.rds"
metadata_train_inputpath <- "data/data_partitions/metadata_train.rds"
metadata_test_inputpath <- "data/data_partitions/metadata_test.rds"
coad_signatures_inputpath <- "extdata/signatures/signatures_coad.rds"

# Outputpath
outputdir <- "data/cox/tcga/signatures/"

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
coad_signatures <- readRDS(coad_signatures_inputpath)

# Convert to gene symbol
signatures_coad_symbol <- list()
for (i in seq_along(coad_signatures)) {
  signatures_coad_symbol[[i]] <- AnnotationDbi::mapIds(
                                                    org.Hs.eg.db,
                                                    keys = coad_signatures[[i]],
                                                    column = "SYMBOL",
                                                    keytype = "ENSEMBL",
                                                    multiVals = "first"
                                                    )                                                  
}
names(signatures_coad_symbol) <- names(coad_signatures)

# Train and test in TCGA-COAD cohort
for (i in seq_along(signatures_coad_symbol)) {

    # Defining signatures
    signature_symbol <- make.names(signatures_coad_symbol[[i]])
    signature_ensembl <- names(signatures_coad_symbol[[i]])
    signature_name <- names(signatures_coad_symbol)[i]
    print(paste0("Runing ", signature_name))

    # Create outputdir
    new_outputdir <- file.path(outputdir, signature_name)
    dir.create(new_outputdir)

    # Define train
    train <- counts_train[, signature_ensembl]
    colnames(train) <- signature_symbol
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
    cox_signature <- coxph(Surv(SurvTime, vital_status) ~ ., data = train %>% dplyr::select(-all_of(cvrts)), x = TRUE)

    # Prepare test data
    test <- counts_test[, signature_ensembl]
    colnames(test) <- signature_symbol
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
    pred_train_signature <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_signature, train[, c(surv_vars, signature_symbol)]), train)

    cindex_train <- data.frame(
        Study = c(paste0(signature_name, " + cvrts"), "cvrts", signature_name),
        C = c(pred_train_all$concordance, pred_train_cvrts$concordance, pred_train_signature$concordance),
        ci_l = c(
            pred_train_all$concordance - 1.96 * pred_train_all$std.err,
            pred_train_cvrts$concordance - 1.96 * pred_train_cvrts$std.err,
            pred_train_signature$concordance - 1.96 * pred_train_signature$std.err),
        ci_u = c(
            pred_train_all$concordance + 1.96 * pred_train_all$std.err,
            pred_train_cvrts$concordance + 1.96 * pred_train_cvrts$std.err,
            pred_train_signature$concordance + 1.96 * pred_train_signature$std.err),
        Subset = "Train"
    )

    # Testing (TCGA-COAD test) ----------------
    # Three scenarios:
    #   Model with yaccs + cvrts
    #   Model with only cvrts
    #   Model with only yaccs

    pred_test_all <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_all, test), test)
    pred_test_cvrts <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_cvrts, test[, c(surv_vars, cvrts)]), test)
    pred_test_signature <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_signature, test[, c(surv_vars, signature_symbol)]), test)

    cindex_test <- data.frame(
        Study = c(paste0(signature_name, " + cvrts"), "cvrts", signature_name),
        C = c(pred_test_all$concordance, pred_test_cvrts$concordance, pred_test_signature$concordance),
        ci_l = c(
            pred_test_all$concordance - 1.96 * pred_test_all$std.err,
            pred_test_cvrts$concordance - 1.96 * pred_test_cvrts$std.err,
            pred_test_signature$concordance - 1.96 * pred_test_signature$std.err),
        ci_u = c(
            pred_test_all$concordance + 1.96 * pred_test_all$std.err,
            pred_test_cvrts$concordance + 1.96 * pred_test_cvrts$std.err,
            pred_test_signature$concordance + 1.96 * pred_test_signature$std.err),
        Subset = "Test"
    )

    # Results of C-Index in TCGA (train and test)
    cindex_tcga <- rbind.data.frame(cindex_train, cindex_test)
    cindex_tcga$signature <- signature_name

    # Save model
    saveRDS(train, file = file.path(new_outputdir, "train_all.rds"))
    saveRDS(test, file = file.path(new_outputdir, "test_all.rds"))

    saveRDS(cox_all, file = file.path(new_outputdir, paste0("cox_cvrts_", signature_name, ".rds")))
    saveRDS(cox_cvrts, file = file.path(new_outputdir, "cox_cvrts.rds"))
    saveRDS(cox_signature, file = file.path(new_outputdir, paste0("cox_", signature_name, ".rds")))

    saveRDS(cindex_tcga, file = file.path(new_outputdir, "cindex_tcga.rds"))
}

