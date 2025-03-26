library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(magrittr)
library(DESeq2)
library(airway)
library(readxl)
library(org.Hs.eg.db, lib.loc = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library")
library(dplyr)
library(tidyverse)

if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}

if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}

if (!("org.Mm.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}

BiocManager::install("pathview")
BiocManager::install("enrichplot")


# phad-ordner erstellen
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots" # Can replace with path to desired output plots directory

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results" # Can replace with path to desired output results directory

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}


# Define the url to your differential expression results file
#dge_url <- "https://refinebio-examples.s3.us-east-2.amazonaws.com/03-rnaseq/results/SRP123625/SRP123625_differential_expression_results.tsv"

dge_results_file <- file.path("~/Desktop/PROMOTION/DEG/DEG ALL/GSEA/DSEQ New.txt")

#download.file(
# dge_url,
# The file will be saved to this location and with this name
#destfile = dge_results_file
#)

# Read in the contents of the differential expression results file
dge_df <- read_excel("~/Desktop/PROMOTION/DEG/MALE/GSEA/DSEQ GSEA MALE.xlsx")

dge_df

mm_hallmark_sets <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "C7")

head(mm_hallmark_sets)
head(dge_df)

#GSEA

keytypes(org.Hs.eg.db)

#wo doppelten?  
dup_gene_Names <- dge_df %>%
  dplyr::filter(duplicated(gene_Names)) %>%
  dplyr::pull(gene_Names)

dge_df %>%
  dplyr::filter(gene_Names %in% dup_gene_Names) %>%
  dplyr::arrange(gene_Names)

#as.numeric(gsub(",", ".", dge_df$log2FoldChange))


#manuell über excel logfold von hoch nach niederig ordnen

filtered_dge_df <- dge_df %>%
  # Sort so that the highest absolute values of the log2 fold change are at the
  # top
  arrange(desc(abs(log2FoldChange)))  %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(gene_Names, .keep_all = TRUE)

any(duplicated(filtered_dge_df$gene_Names))

#ranken
lfc_vector <- filtered_dge_df$log2FoldChange
names(lfc_vector) <- filtered_dge_df$gene_Names

#absteigende reihenfolge
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# Look at first entries of the ranked log2 fold change vector
head(lfc_vector)

# Set the seed so our results are reproducible:
set.seed(5480)

gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 1, # Minimum gene set size
  maxGSSize = 1000, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_hallmark_sets,
    gs_name,
    gene_symbol
  )
)

# We can access the results from our `gsea_results` object using `@result`
head(gsea_results@result)

gsea_result_df <- data.frame(gsea_results@result)
write.table(gsea_result_df, file="gsea_result_male_c7.txt", sep="\t", col.names = NA, quote = F)

#visualisieren:

gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 10)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "FLETCHER_PBMC_BCG_10W_INFANT_PPD_STIMULATED_VS_UNSTIMULATED_10W_UP",
  title = "FLETCHER_PBMC_BCG_10W_INFANT_PPD_STIMULATED_VS_UNSTIMULATED_10W_UP",
  color.line = "#0d76ff")
most_positive_nes_plot

ggsave("ptgs_gsea_enrich_positive_plot.png", width = 20, height = 20, units= "cm")

##neg:
gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(NES, n = 3)

most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "ZAK_PBMC_MRKAD5_HIV_1_GAG_POL_NEF_AGE_20_50YO_1DY_UP",
  title = "ZAK_PBMC_MRKAD5_HIV_1_GAG_POL_NEF_AGE_20_50YO_1DY_UP",
  color.line = "#0d76ff")

most_negative_nes_plot
ggsave("most_negative_nes_plot_ptgs.png", width = 20, height = 20, units= "cm")


readr::write_tsv(
  gsea_result_df,
  file.path(
    results_dir,
    "ptgs_gsea_results_c7.tsv"))

ggplot(gsea_result_df, aes(reorder(Description, NES), NES)) +
  geom_col(aes(fill=p.adjust<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
ggplot2::ggsave("ptgs_gsea_GSEA_barplot_c7.png")  #BARPLOT


dotplot(gsea_results, showCategory= 50)       #DOTPLOT
