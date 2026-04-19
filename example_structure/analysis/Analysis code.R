#Only for the first time using the code run these 4 code lines:
#install.packages("BiocManager")
#BiocManager::install(c("tximport", "DESeq2", "GenomicFeatures","txdbmaker"))
#BiocManager::install("txdbmaker")
#install.packages("ggplot2")

library(GenomicFeatures)
library(txdbmaker)

#Write the annotation file name correctly, the folder path should be already correct
gtf <- "../annotation/AnnotationFileName.gtf.gz"

txdb <- txdbmaker::makeTxDbFromGFF(gtf, format = "gtf")

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

head(tx2gene)

#Change depending on the Sample Names (Salmon folder files names without the .tabular)
samples <- c(
  "WT1","WT2","WT3","WT4","WT5","WT6",
  "5XFAD1","5XFAD2","5XFAD3","5XFAD4","5XFAD5","5XFAD6"
)

files <- file.path("../Salmon", paste0(samples, ".tabular"))
names(files) <- samples

files

library(tximport)

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

#Change depending on the study as the condition compared and the number of samples per condition
condition <- factor(c(
  rep("WT", 6),
  rep("5XFAD", 6)
))

coldata <- data.frame(
  row.names = samples,
  condition = condition
)

library(DESeq2)

#Create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

#Relevel before running DESeq
dds$condition <- relevel(dds$condition, ref = "WT")

#Execute DESeq
dds <- DESeq(dds)

#Results in the example case (5xFAD vs WT)
res <- results(dds)

#Order based on the p-adj value
resOrdered <- res[order(res$padj), ]

#Save as CSV
write.csv(as.data.frame(resOrdered), "DESeq2_results.csv")

#Run PCA to determine if the data is good
vsd <- vst(dds)
plotPCA(vsd, intgroup = "condition")

#Volcano Plot
library(ggplot2)

df <- data.frame(
  log2FC = res$log2FoldChange,
  negLogP = -log10(res$pvalue),
  padj = res$padj
)

#Move the padj and log2FC values in order to change the values considered significative
df$threshold <- "Not sig"
df$threshold[df$padj < 0.05 & df$log2FC > 0] <- "Up"
df$threshold[df$padj < 0.05 & df$log2FC < 0] <- "Down"

#Now generate the volcano plot
ggplot(df, aes(x = log2FC, y = negLogP, color = threshold)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Downregulated in 5XFAD" = "blue", "Not sig" = "grey", "Upregulated in 5XFAD" = "red")) +
  theme_minimal() +
  labs(
    x = "log2(Fold Change)",
    y = "-log10(p-value)",
    color = "Regulation"
  ) +
  theme(
    text = element_text(size = 14),
    legend.position = "right"
  )

#Heatmap run the 3 lines below only the first time to install packages
#install.packages("pheatmap")
#install.packages("org.Mm.eg.db")
#install.packages("dplyr")

library(dplyr)
library(pheatmap)

#Choose the Ensembl codes for the desired genes to visualize
genes_final <- c(
  "ENSMUSG00000023992"
)

genes_present <- genes_final[genes_final %in% rownames(vsd)]

mat <- assay(vsd)[genes_present, ]

mat_scaled <- t(scale(t(mat)))

pheatmap(
  mat_scaled,
  cluster_rows = FALSE,     
  cluster_cols = FALSE,      
  show_rownames = TRUE,
  fontsize_row = 8,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Expression levels - Hipocampus"
)


#To analyze isoform expression (possible with Salmon)
library(tximport)
library(dplyr)
library(tidyr)      
library(ggplot2)
library(tibble)

txi_tx <- tximport(files,
                   type = "salmon",  
                   txOut = TRUE)

tx_counts <- txi_tx$counts

rownames(tx_counts) <- sub("\\..*", "", rownames(tx_counts))
tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)

#Choose the desired gene to analyze
my_gene <- "ENSMUSG00000023992"

tx_my_gene <- tx2gene %>%
  filter(GENEID == my_gene) %>%
  pull(TXNAME)

mat_gene <- tx_counts[tx_my_gene, , drop = FALSE]

df_gene <- as.data.frame(mat_gene) %>%
  rownames_to_column("transcript") %>%
  pivot_longer(cols = -transcript,
               names_to = "sample",
               values_to = "counts") %>%
  left_join(coldata %>% as.data.frame() %>% rownames_to_column("sample"),
            by = "sample")

df_gene_prop <- df_gene %>%
  group_by(sample) %>%
  mutate(gene_total = sum(counts),
         prop = counts / gene_total) %>%
  ungroup()

ggplot(df_gene_prop, aes(x = condition, y = prop, fill = transcript)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(title = paste("Uso de transcritos para", my_gene),
       y = "Proporción dentro del gen",
       x = "") +
  theme_bw()

ggplot(df_gene_prop, aes(x = condition, y = prop, fill = transcript)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = transcript),
              width = 0.15, size = 2, alpha = 0.8) +
  labs(title = paste("Uso de transcritos para", my_gene),
       y = "Proporción dentro del gen",
       x = "") +
  theme_bw()
