library(TCGAbiolinks)
# Obs: The data in the legacy database has been aligned to hg19
query.met.gbm <- GDCquery(
  project = "TCGA-GBM", 
  legacy = TRUE,
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450", 
  barcode = c("TCGA-76-4926-01B-01D-1481-05", "TCGA-28-5211-01C-11D-1844-05")
)
GDCdownload(query.met.gbm)

met.gbm.450 <- GDCprepare(
  query = query.met.gbm,
  save = TRUE, 
  save.filename = "gbmDNAmet450k.rda",
  summarizedExperiment = TRUE
)

dim(met.gbm.450)



