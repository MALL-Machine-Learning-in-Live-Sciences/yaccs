source("requirements.R")
source("R/MSIpred.R")

# Inputpaths
inputpath <- "extdata/geo/GSE29621.rds"
plat_annot_inputpath <- "extdata/geo/gpl570-annot-jetscore.rds"

# Outputpaths
outputdir <- "data/pp/geo/"

# Arguments
geoID <- "GSE29621"
limit <- 365 * 12

# Load data
data <- readRDS(inputpath)
annot <- readRDS(plat_annot_inputpath)

# Format counts
counts <- exprs(data) %>% as.data.frame()

# Format clinical
clinical <-
    pData(data) %>%
    mutate(
        vital_status = case_when(
            `os event:ch1` == "dead" ~ "Dead",
            `os event:ch1` == "alive" ~ "Alive"
        ),
        SurvTime = as.vector(as.numeric(`overall survival (os):ch1`)) * 30,
        Stage = case_when(
            `ajcc staging:ch1` == "Stage 1" ~ "I",
            `ajcc staging:ch1` == "Stage 2" ~ "II",
            `ajcc staging:ch1` == "Stage 3" ~ "III",
            `ajcc staging:ch1` == "Stage 4" ~ "IV"
        ),
        gender = case_when(
            `gender:ch1` == "MALE" ~ "male",
            `gender:ch1` == "FEMALE" ~ "female"
        )
    ) %>%
    dplyr::select(c(gender, Stage, vital_status, SurvTime)) %>%
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
