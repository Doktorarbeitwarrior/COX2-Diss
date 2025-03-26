library("glmnet")
library("survival")

#Heatmap and Clustering
library(readr) #import txt
library(readxl) #import excel
library(tidyverse) # data manipulation
library(limma) # DEG
library(edgeR) # DEG
library(AnnotationDbi) # Converting Entrez Gene ID to Symbol
library(org.Hs.eg.db) # Converting Entrez Gene ID to Symbol
library(ComplexHeatmap) #Making the Heatmap
library(circlize) #colours
library(RColorBrewer) # colours
library(dplyr)


TCGA_clinicaldata <- read_excel("~/Desktop/PROMOTION/TCGA-HNSC Master clinical data_v1 Kopie.xlsx")
PTGS_lasso <- TCGA_clinicaldata %>%
  dplyr::select(Patient_ID, OS, OS_months, Gender,HPV16)%>%
  filter(Gender== "Male")%>%
  filter(HPV16== "Negative")%>%
  dplyr::select(Patient_ID, OS, OS_months)




PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS.time))

FPKM <- read_delim("~/Desktop/promotion/TCGA-HNSC/RNAseq/TCGA_HNSC_GDCHarmonized_HTSeqFPKM.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # FPKM data
GeneID_overlap <- read_excel("~/Desktop/PROMOTION/DEG/MALE/DEG ohne HPV/common genes male ohne hpv.xlsx", sheet= 4)

fpkm_data <- FPKM %>% 
  filter(FPKM$Entrez_Gene_ID %in% GeneID_overlap$GeneID_all) %>%
  as.data.frame(.)

colnames(fpkm_data) <- gsub("-01","",colnames(fpkm_data), fixed=TRUE)
rownames(fpkm_data) <- fpkm_data$Entrez_Gene_ID 
fpkm_data$Entrez_Gene_ID <- NULL
fpkm_data <- t(fpkm_data)
fpkm_data <- as.data.frame(fpkm_data)
fpkm_data$Patient_ID <- c(rownames(fpkm_data))


#Lasso data
#FPKM <- read_delim("~/Desktop/promotion/TCGA-HNSC/RNAseq/TCGA_HNSC_GDCHarmonized_HTSeqFPKM.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # FPKM data
#GeneID_overlap <- read_excel("~/Desktop/PROMOTION/DEG/DEG ALL/Heatmap/GeneID_overlap.xlsx")

#GeneID_overlap <- as.data.frame(GeneID_overlap)

#FPKMVIOLIN <- FPKM%>%
  #filter(FPKM$Entrez_Gene_ID %in% GeneID_overlap$GeneID_all) #mit diesem Filter werden die relevanten Gene von GeneID_all aus dem FPKM Datensatz rausgeholt


#colnames(FPKMVIOLIN) <- gsub("-01", "", colnames(FPKMVIOLIN), fixed = T)
#row.names(FPKMVIOLIN) <- FPKMVIOLIN$Entrez_Gene_ID #save genes as row names
#FPKMVIOLIN <- as.data.frame(FPKMVIOLIN)
#FPKMVIOLIN[1] <- NULL
#FPKMVIOLIN <- as.matrix.data.frame(FPKMVIOLIN)

#FPKMVIOLIN <- log(FPKMVIOLIN+1)
#FPKMVIOLIN <- scale(t(FPKMVIOLIN), center = T, scale = T) #subtracting mean and dividing with SD

# Füge die Clinicaldata mit dem Cluster zusammen

PTGS_lasso <- left_join(PTGS_lasso, fpkm_data, by = "Patient_ID") #adding to your clinicaldata



rt <- PTGS_lasso

rt=rt%>%
  filter(!is.na(rt$'6542'))
   
rt$OS <- as.numeric(rt$OS)

x=as.matrix(rt[,c(4:ncol(rt))])
y=data.matrix(Surv(rt$OS.time,rt$OS))

fit <- glmnet(x, y, family = "cox")
plot(fit, xvar = "lambda", label = FALSE, lw=2)

cvfit <- cv.glmnet(x, y, family="cox")
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
plot(cvfit)

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoefmalealleyearshpv.txt",sep="\t",quote=F,row.names=F)

riskScore=predict(cvfit, newx = x, s = "lambda.min",type="response")
outCol=c("OS.time","OS",lassoGene)
outTab=cbind(rt[,outCol],riskScore=as.vector(riskScore),rt$Patient_ID)
write.table(cbind(id=rownames(outTab),outTab),
            file="lassoRiskmaleallyears.txt",
            sep="\t",
            quote=F,
            row.names=F)

###################################################best cutoff
library(maxstat)
library(survival)
stat <- maxstat.test(Surv(outTab$OS.time,outTab$OS)~outTab$riskScore, data = outTab, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
risk=as.vector(ifelse(riskScore>cutoff,"high","low"))
write.table(cbind(id=rownames(outTab),outTab,risk),
            file="MaxstatlassoRiskmaleallyearsHPV.txt",
            sep="\t",
            quote=F,
            row.names=F)


#SPEARMAN CORRELATION

clinicalrisk <- read_excel("Desktop/PROMOTION/DEG/MALE/Lasso cox:spearman/OS:violin:spearman/all/new/RISK MALE OS UNBE.xlsx")


PTGS <- FPKM
#PTGS$Idx <- NULL
colnames(PTGS) <- gsub("-01","",colnames(PTGS), fixed=TRUE)
PTGS <- PTGS %>%
  filter(!is.na(PTGS$Entrez_Gene_ID))
row.names(PTGS) <- PTGS$Entrez_Gene_ID #save genes as row names
PTGS<-PTGS[c("5743"), ]
row.names(PTGS) <- PTGS$Entrez_Gene_ID #save genes as row names
PTGS <- as.data.frame(PTGS)
PTGS[1] <- NULL
PTGS <- as.matrix.data.frame(PTGS)
# converting to a data frame

#Ln+1 Transformation auf die Daten anwenden
PTGS <- log(PTGS+1)
PTGS <- scale(t(PTGS), center = T, scale = T) #subtracting mean and dividing with SD
PTGS <- as.data.frame(PTGS)
PTGS <-rownames_to_column(PTGS,"Patient_ID")

PTGS <- PTGS%>%
  filter(PTGS$Patient_ID %in% clinicalrisk$Patient_ID)
PTGS <- as.matrix.data.frame(PTGS)
PTGS <- as.data.frame(PTGS)


clinicalrisk <- left_join(PTGS, clinicalrisk, by = "Patient_ID") #adding to your clinicaldata

clinicalrisk <- rename(clinicalrisk, c("5743" = "PTGS"))

clinicalrisk $PTGS <- as.numeric(clinicalrisk $PTGS)


res2 <-cor.test(clinicalrisk$riskScore, clinicalrisk$"PTGS",clinicalrisk$"PTGS",  method = "spearman")
res2


library("ggpubr")
ggscatter(clinicalrisk, x = "riskScore", y = "PTGS", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "riskScore HNSC", ylab = "PTGS Expression HNSC")
