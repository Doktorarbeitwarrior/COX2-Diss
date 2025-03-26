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
library(survival)
library(survminer)
library(dplyr)

#import Dataset
counts <- read_delim("~/Desktop/promotion/TCGA-HNSC/RNAseq/TCGA_HNSC_GDCHarmonized_HTSeqCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # Count data
FPKM <- read_delim("~/Desktop/promotion/TCGA-HNSC/RNAseq/TCGA_HNSC_GDCHarmonized_HTSeqFPKM.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # FPKM data

#import Excel table for overlapping genes
GeneID_overlap <- read_excel("~/Desktop/PROMOTION/DEG/MALE/DEG ohne HPV/common genes male ohne hpv.xlsx", sheet= 4)

#create data frame from excel
GeneID_overlap <- as.data.frame(GeneID_overlap)

#Convert Gene IDs to Symbols
#FPKM$Entrez_Gene_ID <- as.character(FPKM$Entrez_Gene_ID)
#FPKM$Entrez_Gene_ID <- mapIds(org.Hs.eg.db, keys = FPKM$Entrez_Gene_ID, column= "SYMBOL", keytype =  "ENTREZID")

#write.table(FPKM, "FPKM.txt", sep="\t", row.names=FALSE, na="")

#

COX2 <- tibble(Patient_ID = str_replace_all(colnames(FPKM),"-01$", ""),
               "5743" = FPKM %>% 
                 filter(Entrez_Gene_ID == 5743) %>%
                 t(.),
               Sample_ID = colnames(FPKM)) 

COX2 <- mutate(COX2, Group = case_when(`5743` >= quantile(`5743`, probs = 0.75, na.rm = T)~ "high", 
                                       `5743` > quantile(`5743`, probs = 0.25, na.rm = T) & `5743` < quantile(`5743`, probs = 0.75, na.rm = T) ~"moderate",
                                       `5743` <= quantile(`5743`, probs = 0.25, na.rm = T)~ "low"))


#Clustering======
heatmap <- FPKM%>%
  filter(FPKM$Entrez_Gene_ID %in% GeneID_overlap$GeneID_all) #mit diesem Filter werden die relevanten Gene von GeneID_all aus dem FPKM Datensatz rausgeholt

colnames(heatmap) <- gsub("-01", "", colnames(heatmap), fixed = T)
row.names(heatmap) <- heatmap$Entrez_Gene_ID #save genes as row names
heatmap <- as.data.frame(heatmap)
heatmap[1] <- NULL
heatmap <- as.matrix.data.frame(heatmap)

# converting to a data frame

#Ln+1 Transformation auf die Daten anwenden
heatmap <- log(heatmap+1)
heatmap <- scale(t(heatmap), center = T, scale = T) #subtracting mean and dividing with SD

#write.table(heatmap, "heatmap_clustvis.txt", sep="\t", row.names=FALSE, na="")

clustering <- hclust(dist(heatmap, method = "euclidean"), method = "ward.D2") #defining Clustering distance and Clustering method
t <- cutree(clustering, 3) # number of Cluster



cluster <- as.data.frame(t) 
cluster$Patient_ID <- as.factor(rownames(cluster))
cluster <- cluster %>%
  dplyr::rename(`Cluster` = `t`)

# Füge die Clinicaldata mit dem Cluster zusammen

TCGA_clinicaldata <- read_excel("~/Desktop/PROMOTION/TCGA-HNSC Master clinical data_v1 Kopie.xlsx")
TCGA_clinicaldata<- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")

TCGA_clinicaldata<- TCGA_clinicaldata %>%
  filter(Gender== "Male")

TCGA_clinicaldata <- left_join(TCGA_clinicaldata, cluster, by = "Patient_ID") #adding to your clinicaldata
write.table(TCGA_clinicaldata, "TCGA_clinicaldata.txt", sep="\t", row.names=FALSE, na="")

#jetzt joine die neuen Data mit COX2, COX2 Tabelle erstmal anpassen
COX2[2:3] <- NULL 
TCGA_clinicaldata <- left_join(COX2, TCGA_clinicaldata, by = "Patient_ID") #adding to your clinicaldata

#Annotations========


clinical <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$Cluster)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select(OS_5y_event, OS_5y_months,Cluster) %>% # selecting your annotations
  as.data.frame(.) #data frame
clinical <- clinical %>%
  filter(!is.na(clinical$OS_5y_months))
clinical <- clinical %>%
  filter(!is.na(clinical$OS_5y_event))
clinical $OS_5y_months <- as.numeric(clinical $OS_5y_months)
clinical $OS_5y_event <- as.numeric(clinical $OS_5y_event)

KM <- survfit( Surv(OS_5y_months, OS_5y_event) ~ Cluster, data = clinical )


#plot(KM, 
     #conf.int=FALSE,
     #xlab = "OS_5y_months",ylab = "Survival Rate",
     #yscale = 100,
     ##las=1,
     #lwd=2,
     #col =c("green","purple","grey"))

#legend("bottomleft",                                      #man kann auch Koordinaten angeben wie (18, 0.95) oder "bottomleft"
       #title= "Clusters",
       #legend=c("A1","B1a","B1b"),
       #col = c("green","purple","grey"),  #man kann auch die Linien zeichnen lassen mit fill/col=c("red1", "blue"...)
      # lty=1,
       #lwd=2,
       #bty = "",
       #cex=0.6)

ggsurvplot(
  fit = survfit(Surv(OS_5y_months, OS_5y_event) ~ Cluster, data = clinical), 
  xlab = "Days", 
  ylab = "Survival Rate",
  risk.table = TRUE,
  pval = TRUE,
  legend.labs =c("B2", "A","B1"),
  palette = c("green","grey","purple"),
  font.legend = c(25, "bold"),
  legend.title="")



survdiff(Surv(OS_5y_months, OS_5y_event) ~ Cluster,data = clinical)

d <- data.frame(time = KM$time, n.risk = KM$n.risk, n.event = KM$n.event,  n.censor = KM$n.censor, surv = KM$surv,upper = KM$upper, lower = KM$lower)
head(d)
