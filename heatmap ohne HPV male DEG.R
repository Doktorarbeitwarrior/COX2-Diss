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

#import Dataset
counts <- read_delim("~/Desktop/promotion/TCGA-HNSC/RNAseq/TCGA_HNSC_GDCHarmonized_HTSeqCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # Count data
FPKM <- read_delim("~/Desktop/promotion/TCGA-HNSC/RNAseq/TCGA_HNSC_GDCHarmonized_HTSeqFPKM.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # FPKM data

#import Excel table for overlapping genes MALE
GeneID_overlap <- read_excel("~/Desktop/PROMOTION/DEG/MALE/DEG ohne HPV/common genes male ohne hpv.xlsx", sheet= 4)

#create data frame from excel
GeneID_overlap <- as.data.frame(GeneID_overlap)


#Cox2 with Filtering to COX2 Gene and mutate a Column for "Groups" as expression
COX2 <- tibble(Patient_ID = str_replace_all(colnames(FPKM),"-01$", ""),
               "5743" = FPKM %>% 
                 filter(Entrez_Gene_ID == 5743) %>%
                 t(.),
               Sample_ID = colnames(FPKM)) 

COX2 <- mutate(COX2, Group = case_when(`5743` >= quantile(`5743`, probs = 0.75, na.rm = T)~ "high", 
                                       `5743` > quantile(`5743`, probs = 0.25, na.rm = T) & `5743` < quantile(`5743`, probs = 0.75, na.rm = T) ~"moderate",
                                       `5743` <= quantile(`5743`, probs = 0.25, na.rm = T)~ "low"))


#create heatmap from FPKM data and filter with the common genes from DEG
heatmap <- FPKM%>%
  filter(FPKM$Entrez_Gene_ID %in% GeneID_overlap$GeneID_all) #mit diesem Filter werden die relevanten Gene von GeneID_all aus dem FPKM Datensatz rausgeholt

#remove the 01 from end of names and save as data matrix

colnames(heatmap) <- gsub("-01", "", colnames(heatmap), fixed = T)
row.names(heatmap) <- heatmap$Entrez_Gene_ID #save genes as row names
heatmap <- as.data.frame(heatmap)
heatmap[1] <- NULL
heatmap <- as.matrix.data.frame(heatmap)

#Ln+1 Transformation auf die Daten anwenden
heatmap <- log(heatmap+1)
heatmap <- scale(t(heatmap), center = T, scale = T) #subtracting mean and dividing with SD

#add data set for clical data
TCGA_clinicaldata <- read_excel("~/Desktop/PROMOTION/TCGA-HNSC Master clinical data_v1 Kopie.xlsx")
TCGA_clinicaldata<- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")
  
TCGA_clinicaldata<- TCGA_clinicaldata %>%
  filter(Gender== "Male")


#jetzt joine die neuen Data mit COX2, COX2 Tabelle erstmal anpassen
COX2[2:3] <- NULL 
TCGA_clinicaldata <- left_join(COX2, TCGA_clinicaldata, by = "Patient_ID") #adding to your clinicaldata
TCGA_clinicaldata <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$Smoking))

#now filter clinical and heatmap

#heatmap2 <- dplyr::select(heatmap, one_of(rownames(TCGA_clinicaldata$Patient_ID)))
#TCGA_clinicaldata2 <- subset(TCGA_clinicaldata, (rownames(TCGA_clinicaldata) %in% colnames(heatmap)))

#create clinical
clinical <- TCGA_clinicaldata%>%
  dplyr::select(Patient_ID, Group, Smoking, Alcohol, HPV, OS_5y_event, Gender) %>% # selecting your annotations
  as.data.frame(.) #data frame

row.names(clinical) <- clinical$Patient_ID  #Tabelle ändern, Pat.ID ganz nach links
clinical$Patient_ID <- NULL #pat_Id als rownames

#filter heatmap to the same amount of pat in clinical, the order of running commands has to be changed!!!

heatmap <- heatmap[rownames(heatmap) %in% rownames(clinical), ]


#clustering after all the data frames/matrices have been adjusted
clustering <- hclust(dist(heatmap, method = "euclidean"), method = "ward.D2") #defining Clustering distance and Clustering method
t <- cutree(clustering, 3) # number of Cluster

cluster <- as.data.frame(t) 
cluster$Patient_ID <- as.factor(rownames(cluster))
cluster <- cluster %>%
  dplyr::rename(`Cluster` = `t`)

#adding cluster to your clinicaldata
TCGA_clinicaldata <- left_join(TCGA_clinicaldata, cluster, by = "Patient_ID") #adding to your clinicaldata
#write.table(TCGA_clinicaldata, "TCGA_clinicaldata.txt", sep="\t", row.names=FALSE, na="")


#Annotations========

clinical <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$Cluster)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select(Patient_ID, Group, Smoking, Alcohol, OS_5y_event, Cluster) %>% # selecting your annotations
  dplyr::rename(PTGS2 = Group)%>% #optional
  as.data.frame(.) #data frame

write.table(clinical, "clinicaldata_clusterchange.txt", sep="\t", row.names=FALSE, na="")

#after changing cluster names from 123 to A, B1,B2

clinical <- read_excel("~/Desktop/PROMOTION/DEG/MALE/DEG ohne HPV/heatmap for DEG/clinicaldata_clusterchange_excel.xls")
clinical <- clinical %>%
  as.data.frame(.) #data frame

#write.table(clinical,file="GSVA MALE TABLE.txt",sep="\t",quote=F,row.names=F)

row.names(clinical) <- clinical$Patient_ID  #Tabelle ändern, Pat.ID ganz nach links
clinical$Patient_ID <- NULL #Hier stimmt Tabelle icht!!!-> patID fehlt.


#clinical <- clinical[rownames(clinical) %in% rownames(heatmap), ]

#Create Heatmap

mypal <- brewer.pal(3, "Set1") # making colours
mypal2 <- brewer.pal(9, "Set1")

annocolours <- list(
  OS_5y_event = c("1"=mypal[1], "0"=mypal[2], "NA"="black"),
  Smoking = c(Yes=mypal[1],No =mypal[2],"NA"="black" ),
  #HPV16 = c(Positive =mypal[1], Negative= mypal[2],"NA"="black"),
  #Gender  = c(Male =mypal[2], Female= mypal[1]),
  Alcohol= c(Yes =mypal[1],No=mypal[2],"NA"="black"),
  Cluster = c("B2"="purple", "A" = "grey", "B1" = "green"),
  PTGS2 = c(high=mypal[1],low=mypal[2], moderate="black"))


#COX2 <- mutate(COX2, Group = case_when(`5743` >= quantile(`5743`, probs = 0.75, na.rm = T)~ "high", 
#`5743` > quantile(`5743`, probs = 0.25, na.rm = T) & `5743` < quantile(`5743`, probs = 0.75, na.rm = T) ~"moderate",
#`5743` <= quantile(`5743`, probs = 0.25, na.rm = T)~ "low"))

Gene_expression <- GeneID_overlap
Gene_expression <- mutate(Gene_expression, Regulation= case_when(GeneID_overlap$GeneID_all %in% GeneID_overlap$GeneID_low~"low",
                                                                 GeneID_overlap$GeneID_all %in% GeneID_overlap$GeneID_high~"high"))

#jetzt die High und Low columns entfernen

Gene_expression[2:3] <- NULL

#GeneID_all jetzt nach links verschieben und den Column löschen
gene <- Gene_expression
#write.table(gene,file="GSVA MALE TABLE.txt",sep="\t",quote=F,row.names=F)
row.names(gene) <- Gene_expression$GeneID_all
gene$GeneID_all <- NULL 
gene <- as.data.frame(gene)


#Creating annotations
clinical_annotations = HeatmapAnnotation(df = clinical, which = "column", col = annocolours, annotation_name_side = "right")
gene_annotation = HeatmapAnnotation(df = gene , which = "row", col = list(Regulation = c(high=mypal[1], low=mypal[2])))

#Match data
genomic_idx <- match(rownames(gene), colnames(heatmap))
heatmap  <- heatmap[,genomic_idx]

#col_rna = colorRamp2(seq(min(heatmap), max(heatmap), length = 3), c("blue", "#EEEEEE", "red"), space = "LAB", transparency = 0)
#col= colorRamp2(c(min(heatmap, na.rm=T), median(heatmap, na.rm=T), max(heatmap, na.rm=T)), c("blue", "white", "red"))
col_rna = colorRamp2(c(-4, 0, 4), c("blue", "#EEEEEE", "red"), space = "LAB", transparency = 0)




#Heatmap========


HM_rna <- ComplexHeatmap::Heatmap(t(heatmap),
                                  col = col_rna,
                                  name = "Expression",
                                  clustering_distance_rows = "euclidean",
                                  clustering_method_rows = "ward.D2",
                                  cluster_rows = TRUE,
                                  cluster_row_slices = TRUE,
                                  cluster_columns = clustering,
                                  #clustering_distance_columns = "euclidean",
                                  #clustering_method_columns = "ward.D2",
                                  column_split = 3, 
                                  column_title = c("A", "B1", "B2"),
                                  column_title_side = "bottom",
                                  row_split = 2,
                                  row_title = c("", ""),
                                  row_title_rot = 0,
                                  show_column_names = FALSE,
                                  show_row_names = FALSE, 
                                  use_raster = TRUE, 
                                  raster_device = "png",
                                  top_annotation = clinical_annotations,
                                  left_annotation = gene_annotation)

draw(HM_rna)


#### VIOLIN

FPKMVIOLIN <- FPKM%>%
  filter(FPKM$Entrez_Gene_ID %in% GeneID_overlap$GeneID_all) #mit diesem Filter werdendie relevanten Gene von GeneID_all aus dem FPKM Datensatz rausgeholt

FPKMVIOLIN <- tibble(Patient_ID = str_replace_all(colnames(FPKM),"-01$", ""),
                     "5743" = FPKM %>% 
                       filter(Entrez_Gene_ID == 5743) %>%
                       t(.))

FPKMVIOLIN <- FPKMVIOLIN[-c(1), ] 

FPKMVIOLIN <- as.data.frame(FPKMVIOLIN)

colnames(FPKMVIOLIN)[colnames(FPKMVIOLIN) == "5743"] <- "Gene"

#row.names(FPKMVIOLIN) <- FPKMVIOLIN$Patient_ID#save genes as row names
#FPKMVIOLIN <- as.data.frame(FPKMVIOLIN)
#FPKMVIOLIN[1] <- NULL
#FPKMVIOLIN <- as.matrix.data.frame(FPKMVIOLIN)

#create high low table for gene
#COX2 <- tibble(Patient_ID = str_replace_all(colnames(FPKM),"-01$", ""),
#"5743" = FPKM %>% 
#filter(Entrez_Gene_ID == 5743) %>%
# t(.),
#Sample_ID = colnames(FPKM)) 
#COX2 <- mutate(COX2, Group = case_when(`5743` >= quantile(`5743`, probs = 0.75, na.rm = T)~ "high", 
#`5743` > quantile(`5743`, probs = 0.25, na.rm = T) & `5743` < quantile(`5743`, probs = 0.75, na.rm = T) ~"moderate",
#`5743` <= quantile(`5743`, probs = 0.25, na.rm = T)~ "low"))
#COX2[2:3] <- NULL
#COX2 <- COX2[-c(1), ] 

#create cluster table

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




#join 2 tables

Violinplot <- left_join(cluster, FPKMVIOLIN, by = "Patient_ID")

#Filter for high low

Violinplot <-Violinplot %>% 
  dplyr::select("Cluster","Gene")

#save file

write.table(Violinplot, "Violinplotchangedfinal.txt", sep="\t", row.names=FALSE, na="")

#import Excel, which has been adjusted

Violinplotchanged <- read_excel("~/Desktop/PROMOTION/DEG/MALE/DEG ohne HPV/heatmap for DEG/Violon final.xlsx")

#VIOLIN
ggplot(data = Violinplotchanged,aes(x=Cluster  ,fill=Cluster ,y=Gene))+geom_violin()+ylim(0,35)+geom_boxplot(width=0.4)+ theme_classic()+scale_fill_manual(values = c("grey","green","purple"))



x=Violinplotchanged$Cluster

y=Violinplotchanged$Gene

#t.test(y~x, var.equal = TRUE, alternative = "two.sided")


anovatest <- aov(Gene ~ Cluster, data = Violinplotchanged)
# Summary of the analysis
summary(anovatest)
TukeyHSD(anovatest)