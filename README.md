# PRAD
source("http://bioconductor.org/bioClite.R")
biocLite()
biocLite("RTCGAToolbox")
biocLite("devtools")
biocLite("PGSEA")
biocLite("GSEABase")
biocLite("pathview")
library(page)
library(PSEA)
library(RTCGAToolbox)
library(BiocInstaller)
library(devtools)
library(DESeq2)
if(require(RTCGAToolbox)){
  si <- devtools::session_info()
  must.install <- FALSE
  if(!grepl("link-ny", si$packages[si$packages[, 1] == "RTCGAToolbox", 5], ignore.case = TRUE)){
    must.install <- TRUE
  }
}else{
  must.install <- TRUE
}
if(must.install){
  biocLite(c("limma", "RCircos", "data.table", "RCurl", "RJSONIO"))
  biocLite("Link-NY/RTCGAToolbox", "vjcitn/MultiAssayExperiment")
}



library(RTCGAToolbox)
rundates <- getFirehoseRunningDates()
analysisdates <- getFirehoseAnalyzeDates()
prad <- getFirehoseData("PRAD", runDate=rundates[1],
                        gistic2_Date=analysisdates[1], RNAseq_Gene=TRUE, 
                        miRNASeq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE, CNA_SNP=TRUE,
                        CNV_SNP=TRUE, CNA_Seq=TRUE, CNA_CGH=TRUE,  Methylation=TRUE,
                        Mutation=TRUE, mRNA_Array=TRUE, miRNA_Array=TRUE, RPPA=TRUE)



choices <- tolower(gsub("_", "", c("RNAseq_Gene", "miRNASeq_Gene",
                                   "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
                                   "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
                                   "miRNA_Array", "RPPA")))

dses <- lapply(choices, function(choice) try(extract(prad, choice, 
                                                     clinic=TRUE),
                                             silent=TRUE))
names(dses) <- choices
dses

set.rnaseq <- extract(prad, "rnaseq2genenorm")
write.csv(exprs(eset.rnaseq), file="prad_rnaseq.csv")
write.csv(pData(eset.rnaseq), file="prad_clinical.csv")
saveRDS(eset.rnaseq, file="prad_eset.rds")

eset.rnaseq <- readRDS("prad_eset.rds")
hist(exprs(eset.rnaseq["KLK3", ]))
hist(log(exprs(eset.rnaseq["KLK3", ])))
