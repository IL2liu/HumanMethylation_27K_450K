library(TCGAbiolinks)
library(dplyr)
library(DT)
#d
clinical_COAD_query = GDCquery(project = c("TCGA-COAD","TCGA-READ"), data.category = "Clinical")
GDCdownload(clinical_COAD_query)
clinical_coad <- GDCprepare_clinic(clinical_COAD_query, clinical.info = "patient")
write.table(clinical_coad, "/home/smaheshwari/Projects_hdd/TCGA_biolinks_data/data/clinical_COAD_READ.txt", sep="\t")
