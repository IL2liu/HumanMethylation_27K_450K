library(TCGAbiolinks)
library(dplyr)
library(DT)

# colorectal cancer ##########################################################################################

q1 <- GDCquery(project = c("TCGA-COAD","TCGA-READ"), 
               data.category = "DNA Methylation", 
               legacy = FALSE, 
               platform = "Illumina Human Methylation 27")
GDCdownload(q1)
meth_coad <- GDCprepare(q1, summarizedExperiment = FALSE)
write.table(meth_coad, "/home/smaheshwari/Projects/TCGA_biolinks_data/data/meth27_COAD_READ.txt", sep="\t")

# breast cancer ##########################################################################################
q2 <- GDCquery(project = c("TCGA-BRCA"), 
               data.category = "DNA Methylation", 
               legacy = FALSE, platform = "Illumina Human Methylation 27")
GDCdownload(q2)
meth_brca <- GDCprepare(q2, summarizedExperiment = FALSE)
write.table(meth_brca, "/home/smaheshwari/Projects/TCGA_biolinks_data/data/meth27_BRCA_oct2017.txt", sep="\t")

# lung cancer ###############################################################################################
q3 <- GDCquery(project = c("TCGA-LUAD","TCGA-LUSC"), 
               data.category = "DNA Methylation", 
               legacy = FALSE, platform = "Illumina Human Methylation 27")
GDCdownload(q3)
meth_luad <- GDCprepare(q3, summarizedExperiment = FALSE)
write.table(meth_luad, "/home/smaheshwari/Projects/TCGA_biolinks_data/data/meth27_lung_oct2017.txt", sep="\t")

# liver cancer ###############################################################################################
q4 <- GDCquery(project = c("TCGA-LIHC"), 
               data.category = "DNA Methylation", 
               legacy = FALSE, platform = "Illumina Human Methylation 27")
GDCdownload(q4)
meth_lihc <- GDCprepare(q4, summarizedExperiment = FALSE)
write.table(meth_lihc, "/home/smaheshwari/Projects/TCGA_biolinks_data/data/meth27_liver_oct2017.txt", sep="\t")

# kidney cancer  ###############################################################################################
q5 <- GDCquery(project = c("TCGA-KIRC","TCGA-KIRP","TCGA-KICH"), 
               data.category = "DNA Methylation", 
               legacy = FALSE, platform = "Illumina Human Methylation 27")
GDCdownload(q5)
meth_kirc <- GDCprepare(q5, summarizedExperiment = FALSE)
write.table(meth_kirc, "/home/smaheshwari/Projects/TCGA_biolinks_data/data/meth27_kidney_oct2017.txt", sep="\t")

# ovarian cancer  ###############################################################################################
q6 <- GDCquery(project = c("TCGA-OV"), 
               data.category = "DNA Methylation", 
               legacy = FALSE, platform = "Illumina Human Methylation 27")
GDCdownload(q6)
meth_ov <- GDCprepare(q6, summarizedExperiment = FALSE)
write.table(meth_ov, "/home/smaheshwari/Projects/TCGA_biolinks_data/data/meth27_ovarian_oct2017.txt", sep="\t")
