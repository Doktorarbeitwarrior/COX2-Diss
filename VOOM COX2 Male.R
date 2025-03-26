library(readr) #import txt
library(readxl) #import excel
library(tidyverse) # data manipulation
library(limma) # DEG
library(edgeR) # DEG
library(AnnotationDbi) # Converting Entrez Gene ID to Symbol
library(org.Hs.eg.db) # Converting Entrez Gene ID to Symbol


TCGA_HNSC_GDCHarmonized_HTSeqFPKM <- read_delim("~/Desktop/promotion/TCGA-HNSC/RNAseq/TCGA_HNSC_GDCHarmonized_HTSeqFPKM.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # FPKM data

TCGA_HNSC_GDCHarmonized_HTSeqCounts <- read_delim("~/Desktop/promotion/TCGA-HNSC/RNAseq/TCGA_HNSC_GDCHarmonized_HTSeqCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE) # Count data

TCGA_clinicaldata <- read_excel("~/Desktop/PROMOTION/TCGA-HNSC Master clinical data_v1 Kopie.xlsx")

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
colnames(design) <- gsub("x$samples$group", "", colnames(design), fixed = T)
design

#filter low expressed genes before normalization step = common practice
keep <- filterByExpr(x, design = design, group = NULL, lib.size = NULL,
                     min.count = 10, min.total.count = 15)
x <- x[keep,,keep.lib.sizes=FALSE]

#Apply TMM normalization to the count data
x <- calcNormFactors(x, method = "TMM")

#making a contrast matrix, which defines the groups you want to compare, in our case high - low. (But of course you can also compare for example moderate vs low)
contr.matrix <- makeContrasts(
  HvsL = high-low, 
  levels = colnames(design))
contr.matrix

#voom transformation is necessary for RNA-Seq data
v <- voom(x, design, plot=TRUE)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend") ## Anzeige Plots

#Getting a data frame of the results
DEG_COX2_highvslow <- topTable(efit, n=Inf, genelist = rownames(efit), sort.by = "logFC")

#Convert Gene IDs to Symbols
DEG_COX2_highvslow$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = DEG_COX2_highvslow$ID, column= "SYMBOL", keytype =  "ENTREZID")

write.table(DEG_COX2_highvslow, "DEG_COX2_MALE.txt", sep="\t", row.names=FALSE, na="")

#volcanoplot
volcanoplot <- DEG_COX2_highvslow %>%
  mutate(class = case_when(logFC >= 1 & adj.P.Val < 0.05 ~ "Up",
                           logFC <= -1 & adj.P.Val < 0.05 ~ "Down",
                           (logFC < 1 & logFC > -1) | adj.P.Val >= 0.05 ~ "Not Sig"))




DEG_Volcanoplot <- ggplot(data = volcanoplot, mapping = aes(x=logFC, y=-log(adj.P.Val), colour = class)) +
  geom_point() +
  scale_colour_manual(values = c(Up = "green", Down = "red", `Not Sig` = "grey")) +
  geom_vline(xintercept = 1, color ="gray") +
  geom_vline(xintercept = -1, color ="gray") +
  geom_hline(yintercept = -log(0.05), color = "gray") + 
  coord_cartesian(ylim = c(0,50),xlim = c(-6, 6)) + 
  xlab("log2 Fold Change") + 
  ylab("-log10 adj.P") +  
  theme(axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16)) +
  theme(legend.position = "none") 

DEG_vp <- DEG_Volcanoplot + theme_bw() + theme(panel.grid = element_blank())

DEG_vp
ggsave("COX2 VOOM Male.tiff", plot = DEG_Volcanoplot, width = 14, height= 14, units = "cm", dpi = 75)


