source("requirements.R")

# Input and output
inputdir <- "extdata/TCGA-COAD-RNASeq_SumExp.rds"
outdir <- "00_pp/data/"

# Clinical preprocessing
expDat <- readRDS(inputdir)
omics <- assay(expDat)
omics <- as.data.frame(t(omics))

# Get and preprocess clinical variables from TCGAbiolinks
samples <- colData(expDat)
rownames(samples) <- 1:nrow(samples)
ptSamples <- samples[which(samples$definition == 'Primary solid Tumor'), ]
ptSamples <- ptSamples[, c('barcode', 'patient', 'vital_status', 'days_to_last_follow_up', 'days_to_death', 'gender', 'age_at_diagnosis')]
dim(ptSamples)

# Select sample with highest expression mean
SelDupPat <- function(samples, omics){
  
  # Check that sample data have correct variables
  stopifnot('patient' %in% names(samples))
  stopifnot('barcode' %in% names(samples))
  
  # Get patients IDs
  patIDs = samples$patient
  # Return which are duplicated
  dds = unique(patIDs[which(duplicated(patIDs) == T)])
  
  # Iterate throught each duplicate Patient ID
  # and retrieve which one have max mean in omic data
  patDels = list()
  for (i in 1:length(dds)) {
    bcd = samples[which(samples$patient == dds[i]), ]$barcode
    mean = apply(omics[bcd,], 1, function(x) mean(x))
    patDels[[i]] = names(which(mean != max(mean)))
  }
  # Unlist the patients barcode 
  patDels = unlist(patDels)
  res = as.data.frame(samples[! samples$barcode %in% patDels, ])
  
  return(res)
}

clinBiolink = SelDupPat(samples = ptSamples, omics = omics)

# Match and add variables from clinical variables from NIHMS
clinSupp <- read_excel('extdata/NIHMS958067-supplement-2.xlsx', col_names = T, skip = 1)
clinSupp <- clinSupp[which(clinSupp$Organ == 'COAD'),]
clinSupp <- clinSupp[, c('TCGA Participant Barcode' ,'Molecular_Subtype', 'Stage', 'MSI Status')]

# Intersect between two clinical data
int <- intersect(clinSupp$`TCGA Participant Barcode`, clinBiolink$patient)

x <- clinSupp[match(int, clinSupp$`TCGA Participant Barcode`), ]
y <- as.data.frame(clinBiolink[match(int, clinBiolink$patient),])

clinical <- cbind.data.frame(y, x[, c('Molecular_Subtype', 'Stage', 'MSI Status')])

# Defining survival times
follow <- clinical$days_to_last_follow_up
death <- clinical$days_to_death

for (i in 1:length(death)) {
  if(is.na(death[i]) == FALSE){
    follow[i] = death[i]
  }
}
clinical$SurvTime = follow
clinical = subset(clinical, select = -c(days_to_last_follow_up, days_to_death))
clinical = clinical[complete.cases(clinical), ]
clinical = clinical[which(clinical$Stage != 'NA'), ]
clinical = clinical[which(clinical$`MSI Status` != 'NA'), ]
clinical = mutate(clinical, `MSI Status` = ifelse(clinical$`MSI Status` == 'MSI-H', 'MSI-H', 'MSI-L/MSS'))
clinical = mutate(clinical, Molecular_Subtype = ifelse(clinical$Molecular_Subtype == 'CIN', 'CIN', 'noCIN'))
clinical = dummy_cols(clinical, select_columns = c('gender', 'Molecular_Subtype', 'Stage', 'MSI Status'),
                      remove_first_dummy = T, remove_selected_columns = T)
clinical$age_at_diagnosis = round(clinical$age_at_diagnosis / 365, 0)

saveRDS(clinical, file = file.path(outdir, 'clinical.rds'))
