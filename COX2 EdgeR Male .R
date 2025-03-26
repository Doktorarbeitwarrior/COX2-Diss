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
  dplyr::select(Patient_ID, Sample_ID, PTGS, Gender, Group,HPV16) %>% # selecting your annotations
  as.data.frame(.) #data frame

#Filter everything that is NA
COX2 <- COX2 %>%
  filter(!is.na(COX2$Sample_ID))

#now Filter for Male

COX2 <- COX2 %>%
  filter(Gender == "Male")

COX2 <- COX2 %>%
  filter(HPV16 == "Negative")

exprSet <- TCGA_HNSC_GDCHarmonized_HTSeqCounts %>% # For DEG use counts 
  as.data.frame(.)
rownames(exprSet) <- exprSet$Entrez_Gene_ID
exprSet$Entrez_Gene_ID <- NULL

temp <- COX2 %>%
  dplyr::select(Sample_ID, Group)
temp <- as.data.frame(colnames(exprSet)) %>%
  left_join(.,temp, by = c("colnames(exprSet)" = "Sample_ID")) %>%
  na.omit(.)

exprSet <- exprSet %>%
  dplyr::select(temp$`colnames(exprSet)`)

x <- DGEList(counts= exprSet, lib.size = NULL,
             norm.factors =NULL, samples = NULL,
             group = temp$Group, genes = NULL, remove.zeros = FALSE)



design <- model.matrix(~0+x$samples$group)
colnames(design) <- levels(x$samples$group)
design

#filter low expressed genes before normalization step = common practice

keep <- rowSums(cpm(x)>=1) >= 2 #Entweder mussen 2x Werte über 1CPM sein oder man macht´s mit "isexpr <- filterByExpr(d, group=Group)"
x <- x[keep, , keep.lib.sizes=FALSE]
#Apply TMM normalization to the count data
x <- calcNormFactors(x, method = "TMM")


#EdgeR
# Estimate dispersion parameter for GLM

x <- estimateGLMCommonDisp(x, design)

x <- estimateGLMTrendedDisp(x, design, method="power")

x <- estimateGLMTagwiseDisp(x,design)

#plot mean variance

plotBCV(x)

# Design matrix

design <- model.matrix(~ 0 + x$samples$group)
colnames(design) <- c("high", "low", "moderate")


# Model fitting

fit.edgeR <- glmFit(x, design)

# Differential expression

contrasts.edgeR <- makeContrasts(high- low, levels=design)

lrt.edgeR <- glmLRT(fit.edgeR, contrast=contrasts.edgeR)

# Access results tables

edgeR_results <- lrt.edgeR$table

sig.edgeR <- topTags(lrt.edgeR, n = 25300, adjust.method = "BH", sort.by = "PValue")

edgeR_results <- sig.edgeR$table


write.table(edgeR_results, file="edgeR_results.txt", sep="\t", col.names = NA, quote = F)

#add column ID

edgeR_results <- cbind(ID = rownames(edgeR_results), edgeR_results)
edgeR_results
rownames(edgeR_results) <- NULL
edgeR_results
#Convert Gene IDs to Symbols

edgeR_results$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = edgeR_results$ID, column= "SYMBOL", keytype =  "ENTREZID")
write.table(edgeR_results, "EdgeR_COX2_MALE.txt", sep="\t", row.names=FALSE, na="")


# Volcano Plot 

volcanoplot <- edgeR_results %>%
  mutate(class = case_when(logFC >= 1 & FDR < 0.05 ~ "Up",
                           logFC <= -1 & FDR < 0.05 ~ "Down",
                           (logFC < 1 & logFC > -1) | FDR >= 0.05 ~ "Not Sig"))

DEG_Volcanoplot <- ggplot(data = volcanoplot, mapping = aes(x=logFC, y=-log(FDR), colour = class)) +
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
ggsave("EdgeR_COX2_MAlE_VP.tiff", plot = DEG_vp)

