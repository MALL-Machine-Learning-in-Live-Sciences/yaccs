# External validation
# ----------------
source("requirements.R")

# Inputpaths
external_cohorts_inputpath <- "data/pp/geo/"
cox_models_inputpath <- "data/cox/tcga/signatures/"

# Arguments
cvrt <- c(
  "vital_status",
  "SurvTime",
  "gender_male",
  "age_at_diagnosis",
  "Stage_II",
  "Stage_III",
  "Stage_IV",
  "MSI.Status_MSI.L.MSS",
  "Molecular_Subtype_noCIN"
  )

signatures <- c(
  "coloGuideEx",
  "coloGuidePro",
  "coloPrint",
  "mda114",
  "ours",
  "zhang"
  )

cohorts <- c(
  "GSE17536",
  "GSE17537",
  "GSE29621",
  "GSE39582"
  )

# Function to predict in external cohorts
external_prediction <- function(signature, cohort) {

  # Load train
  train <- readRDS(file.path(cox_models_inputpath, signature, "train_all.rds"))

  # Define signature
  genes_signature <- setdiff(colnames(train), cvrt)

  # Load test and format it
  test_clinical <- readRDS(paste0(external_cohorts_inputpath, cohort, "_clinical.rds"))
  names(test_clinical) <- make.names(names(test_clinical))
  test_counts <- readRDS(paste0(external_cohorts_inputpath, cohort, "_counts.rds")) %>% t()
  colnames(test_counts) <- make.names(colnames(test_counts))

  common_genes <- intersect(genes_signature, colnames(test_counts))
  common_cvrts <- intersect(cvrt, names(test_clinical))

  # Define again train set but only with commmon cvrts
  train <- train[, c(common_cvrts, common_genes)]

  # Define test with common genes
  test_counts <- test_counts[, common_genes]
  test <- cbind.data.frame(test_clinical, test_counts)

  # Retrain
  cox_all <- coxph(Surv(SurvTime, vital_status) ~ ., data = train, x = TRUE)
  cox_cvrts <- coxph(Surv(SurvTime, vital_status) ~ ., data = train[,common_cvrts], x = TRUE)
  cox_signature <- coxph(Surv(SurvTime, vital_status) ~ ., data = train[, c("vital_status", "SurvTime", common_genes)], x = TRUE)

  # Batch effect correction
  meta <- data.frame(
    cohort = c(rep("train", nrow(train)), rep("test", nrow(test))),
    patients = c(rownames(train), rownames(test))
  )
  mat <- rbind.data.frame(
    train[, common_genes],
    test[, common_genes]
  )

  modcombat <- model.matrix(~1, data = meta)
  combat <- ComBat(
    dat = t(mat),
    batch = meta$cohort,
    mod = modcombat,
    ref.batch = "train",
    par.prior = TRUE
    ) %>%
    t() %>%
    as.data.frame() %>%
    mutate(
      cohort = meta$cohort
    )

  test[, common_genes] <- combat[which(combat$cohort == "test"), common_genes]
  
  # Format status variable
  test$vital_status <- ifelse(test$vital_status == "Alive", FALSE, TRUE)

  # Prediction
  pred_test_all <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_all, test), test)
  pred_test_cvrts <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_cvrts, test), test)
  pred_test_signature <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox_signature, test), test)

  cindex <- data.frame(
      Study = c(paste0(signature, " + cvrts"), "cvrts", signature),
      C = c(pred_test_all$concordance, pred_test_cvrts$concordance, pred_test_signature$concordance),
      ci_l = c(
          pred_test_all$concordance - 1.96 * pred_test_all$std.err,
          pred_test_cvrts$concordance - 1.96 * pred_test_cvrts$std.err,
          pred_test_signature$concordance - 1.96 * pred_test_signature$std.err),
      ci_u = c(
          pred_test_all$concordance + 1.96 * pred_test_all$std.err,
          pred_test_cvrts$concordance + 1.96 * pred_test_cvrts$std.err,
          pred_test_signature$concordance + 1.96 * pred_test_signature$std.err),
      Subset = cohort,
      signature = signature
  )

  return(cindex)
}

# Iterate across cohorts and signatures
res <- list()
i <- 1
for (s in seq_along(signatures)) {
  for (c in seq_along(cohorts)) {

    print(paste0("Run validaton on: ", signatures[s], " and ", cohorts[c]))

    result <- external_prediction(
      signature = signatures[s],
      cohort = cohorts[c]
    )

    # Store the result in list
    res[[i]] <- result

    i <- i + 1
  }
}

xx <- rbindlist(res)
View(xx)
