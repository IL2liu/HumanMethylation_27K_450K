library(TCGAbiolinks)
library(dplyr)
library(DT)

q1 <- GDCquery(project = c("TCGA-COAD"), 
               data.category = "Transcriptome Profiling",
 	       data.type = "Gene Expression Quantification", 
               #workflow.type = "HTSeq - FPKM-UQ",
               workflow.type = "HTSeq - Counts",
               legacy = FALSE, 
		)
GDCdownload(q1)
rna_coad <- GDCprepare(q1, summarizedExperiment = FALSE)

write.table(rna_coad, "/home/smaheshwari/Projects/Cancer_Indivumed/TCGA_biolinks/download_rnaseq_rawcounts.txt", sep="\t")
