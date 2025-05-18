# miRNA Expression Analysis Pipeline:

This project analyzes miRNA expression in breast cancer by comparing 3 tumor and 3 control samples. It performs adapter trimming, alignment to known miRNA references, read quantification, differential expression analysis using DESeq2, and biological annotation of significant miRNAs.

The goal is to identify key miRNAs that are differentially expressed in tumors and may serve as biomarkers or have functional roles in cancer biology.

This repository provides a complete pipeline for analyzing miRNA expression from raw FASTQ files through differential expression analysis and visualization, using Bash, Bowtie, R (DESeq2, EnhancedVolcano), and Python (pandas).

---

## 🔬 Steps Overview

### 1. Trimming Adapters

```bash
samples=("Tumor1" "Tumor2" "Tumor3" "Control1" "Control2" "Control3")

for sample in "${samples[@]}"; do
    cutadapt -a TGGAATTCTCGGGTGCCAAGG \
             -m 18 -M 30 \
             -o trimmed_${sample}.fastq \
             ${sample}.fastq > ${sample}_cutadapt.log
done
```

### 2. Reference Indexing

```bash
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
gunzip mature.fa.gz
# Rename to match expected name
mv mature.fa hsa_mature.fa
bowtie-build hsa_mature.fa hsa_mature_index
```

### 3. Alignment

```bash
mkdir -p aligned

for sample in tumor1 tumor2 tumor3 control1 control2 control3; do
    bowtie -v 1 -S hsa_mature_index trimmed_${sample}.fastq > aligned/${sample}.sam
done
```

### 4. Extract Mapped Reads & Convert to BAM

```bash
mkdir -p bam

for sample in tumor1 tumor2 tumor3 control1 control2 control3; do
    samtools view -h -F 4 aligned/${sample}.sam > bam/${sample}_mapped.sam
    samtools view -S -b bam/${sample}_mapped.sam > bam/${sample}_mapped.bam
    samtools sort bam/${sample}_mapped.bam -o bam/${sample}_sorted.bam
    samtools index bam/${sample}_sorted.bam
done
```

### 5. Quantification with featureCounts

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz
awk '$0 ~ /gene_type "miRNA"/' gencode.v44.annotation.gtf > hsa_miRNA.gtf

featureCounts \
  -a references/hsa_miRNA.gtf \
  -o results/miRNA_counts.txt \
  -t exon \
  -g gene_id \
  results/bam/*.bam
```

### 6. Filtering with Python (Pandas)

```python
import pandas as pd

df = pd.read_csv("results/miRNA_counts.txt", sep="\t", skiprows=1, header=None)
count_data = df.iloc[:, 6:].apply(pd.to_numeric, errors='coerce')
count_data.index = df.iloc[:, 0]
filtered = count_data[count_data.max(axis=1) >= 10]
print(filtered.shape)
print(filtered.head())
```

### 7. Differential Expression with DESeq2 (R)

```R
library(DESeq2)

counts <- read.table("counts.txt", header=TRUE, row.names=1)
metadata <- read.table("metadata.txt", header=TRUE, row.names=1)
all(rownames(metadata) == colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

resSig <- res[which(res$padj < 0.05), ]
```

### 8. Volcano Plot with ggplot2 or EnhancedVolcano

```R
# ggplot2 volcano
library(ggplot2)
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05, "Significant", "Non-Significant")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Non-Significant" = "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Significant vs Non-Significant miRNAs",
       x = "log2 Fold Change", y = "-log10 Adjusted P-Value")

# EnhancedVolcano
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj')
```

### 9. Annotation with biomaRt (R)

```R
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- gsub("\\..*", "", rownames(resSig))
annotations <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                     filters = "ensembl_gene_id",
                     values = ensembl_ids,
                     mart = ensembl)
resSig$ensembl_gene_id <- gsub("\\..*", "", rownames(resSig))
res_annotated <- merge(resSig, annotations, by = "ensembl_gene_id")
```

### 🔥 Heatmap of Top Significant miRNAs

```R
if (!require(pheatmap)) install.packages("pheatmap")
library(pheatmap)

significant_mirnas_expression <- counts_data[rownames(counts_data) %in% res_annotated$ensembl_gene_id, ]
significant_mirnas_expression_log <- log2(significant_mirnas_expression + 1)

pheatmap(significant_mirnas_expression_log,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = coldata,
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "row",
         main = "Heatmap of Significant miRNAs")
```

---

## 📁 Folder Structure

```
tumor_mirna/
├── references/                 # Reference genome and annotation
├── trimmed/                   # Trimmed fastq files
├── aligned/                   # SAM files from Bowtie
├── bam/                       # BAM files (sorted, indexed)
├── results/                   # Quantification and plots
├── scripts/                   # Custom analysis scripts (Python/R)
└── README.md
```

## 📦 Dependencies

* fastq
* cutadapt
* bowtie /
* samtools
* featureCounts (subread)
* pandas (Python)
* DESeq2, ggplot2, EnhancedVolcano, biomaRt, pheatmap (R)

---

Feel free to contribute, fork, or open an issue. 🎓
