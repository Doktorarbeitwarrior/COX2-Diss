library(readr)
library(edgeR)
library(limma)
library(dplyr)
library(Biobase)
library(readr) #import txt
library(readxl) #import excel
library(tidyverse) # data manipulation
library(AnnotationDbi) # Converting Entrez Gene ID to Symbol
library(org.Hs.eg.db) # Converting Entrez Gene ID to Symbol
library(DESeq2)

TCGA_HNSC_GDCHarmonized_HTSeqFPKM <- read_delim("~/Desktop/promotion/TCGA-HNSC/RNAseq/TCGA_HNSC_GDCHarmonized_HTSeqFPKM.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # FPKM data

TCGA_HNSC_GDCHarmonized_HTSeqCounts <- read_delim("~/Desktop/promotion/TCGA-HNSC/RNAseq/TCGA_HNSC_GDCHarmonized_HTSeqCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # Count data

TCGA_clinicaldata <- read_excel("~/Desktop/PROMOTION/TCGA-HNSC Master clinical data_v1.xlsx")

COX2 <- tibble(Patient_ID = str_replace_all(colnames(TCGA_HNSC_GDCHarmonized_HTSeqFPKM),"-01$", ""),
               "5743" = TCGA_HNSC_GDCHarmonized_HTSeqFPKM %>% 
                 filter(Entrez_Gene_ID == 5743) %>%
                 t(.),
               Sample_ID = colnames(TCGA_HNSC_GDCHarmonized_HTSeqFPKM)) 

COX2 <- mutate(COX2, Group = case_when(`5743` >= quantile(`5743`, probs = 0.75, na.rm = T)~ "high", 
                                       `5743` > quantile(`5743`, probs = 0.25, na.rm = T) & `5743` < quantile(`5743`, probs = 0.75, na.rm = T) ~"moderate",
                                       `5743` <= quantile(`5743`, probs = 0.25, na.rm = T)~ "low"))
COX2 <- COX2 %>%
 dplyr::rename("PTGS" = "5743")
#Join Clinical data to add gender
COX2 <- left_join(TCGA_clinicaldata, COX2, by = "Patient_ID")

COX2 <- COX2 %>%
  dplyr::select(Patient_ID, Sample_ID, PTGS, Gender, Group) %>% # selecting your annotations
  as.data.frame(.) #data frame

#Filter everything that is NA
COX2 <- COX2 %>%
  filter(!is.na(COX2$Sample_ID))

#now Filter for Male

COX2 <- COX2 %>%
  filter(Gender == "Male")

#epression set 
exprSet <- TCGA_HNSC_GDCHarmonized_HTSeqCounts %>% # For DEG use counts 
  as.data.frame(.)
rownames(exprSet) <- exprSet$Entrez_Gene_ID
exprSet$Entrez_Gene_ID <- NULL

#now remove female from Expressionset
#exprSet <- exprSet%>%
  #filter(exprSet$Entrez_Gene_ID %in% GeneID_overlap$GeneID_all) #mit diesem Filter werden die relevanten Gene von GeneID_all aus dem FPKM Datensatz rausgeholt

temp <- COX2 %>%
  dplyr::select(Sample_ID, Group)
temp <- as.data.frame(colnames(exprSet)) %>%
  left_join(.,temp, by = c("colnames(exprSet)" = "Sample_ID")) %>%
  na.omit(.)

exprSet <- exprSet %>%
  dplyr::select(temp$`colnames(exprSet)`)

#Add Column at the beginning for reference for Pat ID and change Column name for understanding

(colData <- data.frame(row.names=colnames(exprSet), temp = temp ))
colnames(colData) <- c("PatID", "Group")



# Create DESeq2 datasets
x <- DESeqDataSetFromMatrix(countData = exprSet, colData = colData, 
                              design = ~ Group )

keep <- rowSums(counts(x)) >= 10
x <- x[keep,]

x <- DESeq(x)


# Plot dispersion estimates
plotDispEsts(x)

#define treshold

p.threshold <- 0.05

## DESeq2 ##

rescoxdq2 <- results(x, contrast=c("Group", "high","low"))
rescoxdq2$threshold <- as.logical(rescoxdq2$padj < p.threshold)
genes.deseq <- row.names(rescoxdq2)[which(rescoxdq2$threshold)]

resOrdered <- rescoxdq2[order(rescoxdq2$pvalue),]
resOrdered

DEG = as.data.frame(resOrdered)
COX_DESeq2 = na.omit(DEG)

#write table

write.table(COX_DESeq2, file="COX2_DESeq2_DEG.txt", sep="\t", col.names = NA, quote = F)

#add column ID

COX_DESeq2 <- cbind(ID = rownames(COX_DESeq2), COX_DESeq2)
COX_DESeq2
rownames(COX_DESeq2) <- NULL
COX_DESeq2

#Convert Gene IDs to Symbols

COX_DESeq2$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = COX_DESeq2$ID, column= "SYMBOL", keytype =  "ENTREZID")
write.table(COX_DESeq2, "COX_DESeq2_DEG_FINAL.txt", sep="\t", row.names=FALSE, na="")

#Volcano

volcanoplot <- COX_DESeq2 %>%
  mutate(class = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Up",
                           log2FoldChange <= -1 & padj < 0.05 ~ "Down",
                           (log2FoldChange < 1 & log2FoldChange > -1) | padj >= 0.05 ~ "Not Sig"))

DEG_Volcanoplot <- ggplot(data = volcanoplot, mapping = aes(x=log2FoldChange, y=-log(padj), colour = class)) +
  geom_point() +
  scale_colour_manual(values = c(Up = "green", Down = "red", `Not Sig` = "grey")) +
  geom_vline(xintercept = 1, color ="gray") +
  geom_vline(xintercept = -1, color ="gray") +
  geom_hline(yintercept = -log(0.05), color = "gray") + 
  coord_cartesian(ylim = c(0,60),xlim = c(-7, 7)) + 
  xlab("log2 Fold Change") + 
  ylab("-log10 adj.P") +  
  theme(axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16)) +
  theme(legend.position = "none") 

DEG_vp <- DEG_Volcanoplot + theme_bw() + theme(panel.grid = element_blank())

DEG_vp
ggsave("DEG_COX_DESeq2_MALE.tiff", plot = DEG_Volcanoplot, width = 14, height= 14, units = "cm", dpi = 75)