source("requirements.R")

outdir <- "00_pp/data/"

# Data download by barcode (RNASeq and clinical data)
projID <- "TCGA-COAD"

# Expression Data (RNA-Seq)
query <- GDCquery(project = projID, 
                 data.category = 'Transcriptome Profiling',
                 data.type = 'Gene Expression Quantification',
                 experimental.strategy = 'RNA-Seq',
                 workflow.type = 'STAR - Counts')

GDCdownload(query, directory = outdir)
expDat <- GDCprepare(query = query,
                    save = TRUE,
                    directory = outdir,
                    save.filename = "00_preprocessing/data/TCGA-COAD-RNASeq")

saveRDS(expDat, file = paste0(outdir, 'TCGA-COAD-RNASeq_SumExp.rds'))
