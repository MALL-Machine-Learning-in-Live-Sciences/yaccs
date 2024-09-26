source("requirements.R")

# Inputpaths
inputpath <- "extdata/geo/GSE17536.rds"
plat_annot_inputpath <- "extdata/geo/gpl570-annot-jetscore.rds"

# Outputpaths
outputdir <- "data/pp/geo/"

# Arguments
geoID <- "GSE17536"
limit <- 365 * 12

# Load data
data <- readRDS(inputpath)
annot <- readRDS(plat_annot_inputpath)

# Format counts
counts <-
    exprs(data) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("probeID") %>%
    left_join(., annot[, c("probeID", "symbol", "overall")], by = "probeID") %>%
    group_by(symbol) %>%
    slice_max(order_by = overall, with_ties = FALSE) %>%
    filter(
        !is.na(symbol)
    ) %>%
    tibble::column_to_rownames("symbol") %>%
    dplyr::select(-c(probeID, overall))

# Format clinical
clinical <-
    pData(data) %>%
    mutate(
        vital_status = case_when(
            `overall_event (death from any cause):ch1` == "death" ~ "Dead",
            `overall_event (death from any cause):ch1` == "no death" ~ "Alive"
        ),
        SurvTime = as.vector(as.numeric(`overall survival follow-up time:ch1`)) * 30,
        age_at_diagnosis = as.vector(as.numeric(`age:ch1`)),
        Stage = case_when(
            `ajcc_stage:ch1` == "1" ~ "I",
            `ajcc_stage:ch1` == "2" ~ "II",
            `ajcc_stage:ch1` == "3" ~ "III",
            `ajcc_stage:ch1` == "4" ~ "IV"
        ),
        gender = `gender:ch1`
    ) %>%
    dplyr::select(c(age_at_diagnosis, gender, Stage, vital_status, SurvTime)) %>%
    filter(
        complete.cases(.)
    ) %>%
    # Relabeling survival times
    mutate(
        SurvTime = if_else(SurvTime > limit, limit, SurvTime),
        vital_status = if_else(SurvTime > limit, "Alive", vital_status)
    )

# Match patients
counts <- counts[, rownames(clinical)]
stopifnot(colnames(counts) == rownames(clinical))

# Predict MSI variable
msi <- predMSI(counts)
clinical$MSI.Status <- msi$msi

# Convert into dummy variables
clinical <- fastDummies::dummy_columns(
    clinical,
    select_columns = c("gender", "Stage", "MSI.Status"),
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE)

# Save data
saveRDS(counts, file = file.path(outputdir, paste0(geoID, "_counts.rds")))
saveRDS(clinical, file = file.path(outputdir, paste0(geoID, "_clinical.rds")))