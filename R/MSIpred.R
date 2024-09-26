# install.packages("remotes")
# remotes::install_github("WangX-Lab/PreMSIm")
# see (https://github.com/WangX-Lab/PreMSIm/)

predMSI <- function(counts){
  
  require(PreMSIm)
  require(dplyr)

  # Define signature
  msiSig <- c('DDX27', 'EPM2AIP1', 'HENMT1', 'LYG1', 'MLH1', 'MSH4', 'NHLRC1', 'NOL4L', 'RNLS', 'RPL22L1', 'RTF2', 'SHROOM4', 'SMAP1', 'TTC30A', 'ZSWIM3')
  msiSig[grep('RTF2', msiSig)] <- 'SLC52A3'
  msiSig[grep('NHLRC1', msiSig)] <- 'EPM2A'

  # Subset counts by signature
  counts_ <-
    counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("symbol") %>%
    filter(symbol %in% msiSig) %>%
    tibble::column_to_rownames("symbol")

  # Format gene names
  rownames(counts_)[grep("SLC52A3", rownames(counts_))] <- "RTF2"
  rownames(counts_)[grep("\\bEPM2A\\b", rownames(counts_))] <- "NHLRC1"

  # Save counts to file in tmp dir
  tmppath <- "tmp/temporary_file.txt"
  if (file.exists(tmppath)) {
    file.remove(tmppath)
  }
  write.table(counts_, file =  tmppath, sep = "\t")

  # Calculate MSI
  system.file("extdata", "example.txt", package = "PreMSIm", mustWork = TRUE)
  premsi <- PreMSIm::data_pre(tmppath, type = "Symbol")
  res <-
    PreMSIm::msi_pre(premsi) %>%
    mutate(
      msi = ifelse(MSI_status == 1, "MSI-H", "MSI-L/MSS")
    ) %>%
    tibble::column_to_rownames("Sample")
  
  return(res)
}
