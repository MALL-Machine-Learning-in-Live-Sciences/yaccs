# External validation
# ----------------
source("requirements.R")

# Inputpaths
external_cohorts_inputpath <- "data/pp/geo/"
cox_models_inputpath <- "data/cox/tcga/signatures/"

# Outputpath
outputpath <- "data/cox/external/cindex_external_validations.rds"

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

# TO REMOVE
signature <- "coloPrint"
cohort <- "GSE29621"


# Function to predict in external cohorts
external_prediction <- function(signature, cohorts) {

  # Load train
  train <- readRDS(file.path(cox_models_inputpath, signature, "train_all.rds"))

  # Define signature
  genes_signature <- setdiff(colnames(train), cvrt)

  # Load cohorts and format them
  ext_cohorts_counts <- list()
  ext_cohorts_clin <- list()
  cox <- list()

  cohort <- cohorts[1]
  for (cohort in cohorts) {

    test_clinical <- readRDS(paste0(external_cohorts_inputpath, cohort, "_clinical.rds"))
    names(test_clinical) <- make.names(names(test_clinical))
    test_counts <- readRDS(paste0(external_cohorts_inputpath, cohort, "_counts.rds")) %>% t()
    colnames(test_counts) <- make.names(colnames(test_counts))

    common_genes <- intersect(genes_signature, colnames(test_counts))
    common_cvrts <- intersect(cvrt, names(test_clinical))

    # Define again train set but only with commmon cvrts
    train_ <- train[, c(common_cvrts, common_genes)]
    
    # Define test with common genes
    test_counts <- test_counts[, common_genes]
    test <- cbind.data.frame(test_clinical, test_counts)

    # Retrain
    cox_all <- coxph(Surv(SurvTime, vital_status) ~ ., data = train_, x = TRUE)
    cox_cvrts <- coxph(Surv(SurvTime, vital_status) ~ ., data = train_[,common_cvrts], x = TRUE)
    cox_signature <- coxph(Surv(SurvTime, vital_status) ~ ., data = train_[, c("vital_status", "SurvTime", common_genes)], x = TRUE)

    # Store cox models
    cox[[cohort]] <- list(
        all = cox_all,
        cvrts = cox_cvrts,
        signature = cox_signature
    )

    # Store counts of external datasets
    ext_cohorts_counts[[cohort]] <- 
        test_counts %>%
        as.data.frame() %>%
        mutate(
            cohort = cohort
        )

    # Store cvrts of external datasets
    ext_cohorts_clin[[cohort]] <- test_clinical
  }

  # Merge all cohorts to run batch correction
  mat <- rbind.data.frame(
    train[, common_genes] %>% mutate(cohort = "tcga"), 
    ext_cohorts_counts[["GSE17536"]],
    ext_cohorts_counts[["GSE17537"]],
    ext_cohorts_counts[["GSE29621"]],
    ext_cohorts_counts[["GSE39582"]]
    )

  meta <- data.frame(
    cohort = mat$cohort,
    row.names = rownames(mat)
  )

  mat <- mat %>% dplyr::select(-c(cohort))

  modcombat <- model.matrix(~1, data = meta)
  combat <- ComBat(
    dat = t(mat),
    batch = meta$cohort,
    mod = modcombat,
    ref.batch = "tcga",
    par.prior = TRUE
    ) %>%
    t() %>%
    as.data.frame() %>%
    mutate(
      cohort = meta$cohort
    )
    
  # Retrieve counts corrected
  GSE17536 <- combat %>% filter(cohort == "GSE17536") %>% dplyr::select(-c(cohort))
  GSE17537 <- combat %>% filter(cohort == "GSE17537") %>% dplyr::select(-c(cohort))
  GSE29621 <- combat %>% filter(cohort == "GSE29621") %>% dplyr::select(-c(cohort))
  GSE39582 <- combat %>% filter(cohort == "GSE39582") %>% dplyr::select(-c(cohort))

  # Add cvrts to external data
  GSE17536 <- cbind.data.frame(ext_cohorts_clin$GSE17536, GSE17536) %>% mutate(vital_status = ifelse(vital_status == "Alive", FALSE, TRUE))
  GSE17537 <- cbind.data.frame(ext_cohorts_clin$GSE17537, GSE17537) %>% mutate(vital_status = ifelse(vital_status == "Alive", FALSE, TRUE))
  GSE29621 <- cbind.data.frame(ext_cohorts_clin$GSE29621, GSE29621) %>% mutate(vital_status = ifelse(vital_status == "Alive", FALSE, TRUE))
  GSE39582 <- cbind.data.frame(ext_cohorts_clin$GSE39582, GSE39582) %>% mutate(vital_status = ifelse(vital_status == "Alive", FALSE, TRUE))

  external_cohorts <- list(GSE17536, GSE17537, GSE29621, GSE39582)
  names(external_cohorts) <- cohorts
  
  cohort <- cohorts[1]
  cindex <- list()
  # Predict across cohorts  
  for (cohort in cohorts) {

    # Prediction
    pred_test_all <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox[[cohort]]$all, external_cohorts[[cohort]]), external_cohorts[[cohort]])
    pred_test_cvrts <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox[[cohort]]$cvrts, external_cohorts[[cohort]]), external_cohorts[[cohort]])
    pred_test_signature <- survival::survConcordance(Surv(SurvTime, vital_status) ~ predict(cox[[cohort]]$signature, external_cohorts[[cohort]]), external_cohorts[[cohort]])

    cindex[[cohort]] <- data.frame(
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
  }
  
  performances <- data.table::rbindlist(cindex)

  return(performances)
}


res <- lapply(signatures, function(x) external_prediction(x, cohorts))
res <- data.table::rbindlist(res)

saveRDS(res, file = outputpath)

