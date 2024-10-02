source("requirements.R")

outdir <- "extdata/tcga-coad/"

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
                    save = FALSE,
                    directory = outdir,
                    save.filename = "TCGA-COAD_counts")

saveRDS(expDat, file = paste0(outdir, 'TCGA-COAD-RNASeq_counts.rds'))



