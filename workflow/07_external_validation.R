# External of yaccs in external cohorts
# ---------

source("requirements.R")
source("R/MSIpred.R")

# Inputpaths
inputdir <- "extdata/geo/"


# Outputpaths
# outDir = '~/projects/SurvivalTCGA-COAD/data/Validation/validation_datasets_signatures/'

# Arguments
model = 'all'
geoID = 'GSE17536' #GSE17536 #GSE17537
sig = 'zhang' #ColoGuideEx #ColoGuidePro #MDA114 #Ours #ColoPrint #zhang
model = 'all'
maxY = 12

file <- list.files(pattern = paste0(geoID, '_', sig))

# GSE17536
gse <- readRDS(file)
signature = gse$annotation$ensembl_gene_id

# MSI prediction
msi = predMSI(geoID = geoID, 'GPL570')

# Omic
omic = gse$counts

# Clinical
clinical = gse$clinical
surVars = c("overall_event (death from any cause):ch1", "overall survival follow-up time:ch1")
cvrts = c('age:ch1', 'gender:ch1', 'ajcc_stage:ch1')

cvrts = clinical[, intersect(names(clinical), cvrts)]
surVars = clinical[, intersect(names(clinical), surVars)]

# Change Names of Variables
nms = c('age_at_diagnosis', 'gender', 'Stage')
names(cvrts)[1] = 'age_at_diagnosis'
names(cvrts)[2] = 'Stage'
names(cvrts)[3] = 'gender'

sv = c('vital_status', 'SurvTime')
names(surVars)[1] = 'SurvTime'
names(surVars)[2] = 'vital_status'

# Process cvrts
cvrts$age_at_diagnosis = as.vector(as.numeric(cvrts$age_at_diagnosis))
cvrts$Stage = ifelse(cvrts$Stage == '1', 'I', ifelse(cvrts$Stage == '2', 'II', ifelse(cvrts$Stage == '3', 'III', ifelse(cvrts$Stage == '4', 'IV', NA))))

cvrts = cvrts[, match(nms, names(cvrts))]

# Process survival variables
surVars$SurvTime = as.vector(as.numeric(surVars$SurvTime)) * 30
surVars$vital_status = ifelse(surVars$vital_status == 'death', 'Dead', ifelse(surVars$vital_status == 'no death', 'Alive', NA))

# Complete cases
cc = cbind.data.frame(cvrts, MSI.Status = msi$msi, surVars)
cc = cc[complete.cases(cc), ]

omic = omic[match(rownames(cc), rownames(omic)), ]

omic = omic[, match(gse$annotation$affy_hg_u133_plus_2, names(omic))]
names(omic) = gse$annotation$ensembl_gene_id

val = cbind.data.frame(cc[, intersect(c(nms, 'MSI.Status'), names(cc))], omic, cc[, match(sv, names(cc))])
rnms = rownames(val)
val = dummy_columns(val, select_columns = c('gender', 'Stage', 'MSI.Status'), remove_first_dummy = T, remove_selected_columns = T)
rownames(val) = rnms

names(val) = make.names(names(val))

ll = labelData(years = maxY, clinical = val)
ll$vital_status = ifelse(ll$vital_status == 'Dead', TRUE, FALSE)

res = Scaling(signature, ll)
#res = ll

saveRDS(res, file = paste0(outDir, geoID, '-', model, '-', sig, '.rds'))
