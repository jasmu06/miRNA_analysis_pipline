samples=("Tumor1" "Tumor2" "Tumor3" "Control1" "Control2" "Control3")

for sample in "${samples[@]}"; do
    cutadapt -a TGGAATTCTCGGGTGCCAAGG \
             -m 18 -M 30 \
             -o trimmed_${sample}.fastq \
             ${sample}.fastq > ${sample}_cutadapt.log
done

# Reference preparation

wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
gunzip mature.fa.gz
The file was unzipped as hsa_mature.fa and process for indexing. Mature miRNA sequences from miRbase and tool Bowtie tool was carried out.

bowtie-build mature_human_miRNA.fa mature_miRNA_index

#Alignmemnt

Alignment loop for all samples:
for sample in tumor1 tumor2 tumor3 control1 control2 control3
do
    bowtie -v 1 -S hsa_mature_index trimmed/${sample}.fastq > aligned/${sample}.sam
done

# extracted maapped alignment

for sample in Tumor1 Tumor2 Tumor3 Control1 Control2 Control3; do
    samtools view -h -F 4 "${sample}_trimmed.sam" > "${sample,}_mapped.sam"
done

#sam to bam and sorted

for sample in tumor1 tumor2 tumor3 controll1 control2 cotrol3; do
samtools view -S -b ${sample}_mapped.sam > ${sample}_mapped.bam
samtools sort ${sample}_mapped.bam -o ${sample}_sorted.bam
samtools index ${sample}_sorted.bam
done

# For extraction of only mapped reads:
for sample in tumor1 tumor2 tumor3 control1 control2 control3
do
    echo "Processing $sample..."
    samtools view -h -F 4 ${sample}_trimmed.sam > ${sample}_mapped.sam
samtools view -Sb ${sample}_mapped.sam | samtools sort -o ${sample}_sorted.bam
done

# quantification 

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz

awk '$0 ~ /gene_type "miRNA"/' gencode.v44.annotation.gtf > hsa_miRNA.gtf

featureCounts \
  -a /Users/jasmu/tumor_mirna/references/hsa_miRNA.gtf \
  -o /Users/jasmu/tumor_mirna/results/miRNA_counts.txt \
  -t exon \
  -g gene_id \
  /Users/jasmu/tumor_mirna/results/aligned/*.bam

# filtering data by used panda in python

import pandas as pd

# Load the file (adjust path and skiprows as needed)
df = pd.read_csv("results/miRNA_counts.txt", sep="\t", skiprows=1, header=None)

# Keep only the count columns (assuming gene_id is column 0 and counts start from column 6)
count_data = df.iloc[:, 6:]

# Convert all columns to numeric (coerce errors just in case)
count_data = count_data.apply(pd.to_numeric, errors='coerce')

# Set gene_id as row index
count_data.index = df.iloc[:, 0]

# Now filter genes with at least 10 reads in any sample
filtered = count_data[count_data.max(axis=1) >= 10]

# Optional: Check the result
print(filtered.shape)
print(filtered.head())

# Differential gene expression uses DEseq in R

library(DESeq2)

# Read in count matrix and sample metadata
counts <- read.table("counts.txt", header=TRUE, row.names=1)
metadata <- read.table("metadata.txt", header=TRUE, row.names=1)

# Ensure that the columns in counts match the rows in metadata
all(rownames(metadata) == colnames(counts))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              metaData = metadata,
                              design = ~ condition)

# Run the differential expression analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# data was filtered to eliminate non-numeric characters and low-expressed 
# Filter to keep miRNAs with at least 10 reads in at least one sample
filtered_data <- res[rowSums(res[, -1] >= 10) > 0, ]
# Convert the count data column to numeric
res$log2FoldChange <- as.numeric(res$log2FoldChange)

# Differentially expressed genes were identified by applying filtered and adjusted p-value (padj) threshold of < 0.05 using DESeq2 
res <- results(dds)
resSig <- res[which(res$padj < 0.05), ].

# Visualization of Differential Expression and Creation of Volcano and MA Plots

# MA plot
 plotMA(res, main="MA Plot - Tumor vs Control", ylim=c(-2,2))
  pdf("MA_plot.pdf", width=800, height=600)
  plotMA(res, main="MA Plot - Tumor vs Control", ylim=c(-2,2))
  dev.off()

#Generate volcano plot
library(ggplot2)
ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5, color = "blue") +
  theme_minimal() +
  ggtitle("Volcano Plot") +
  xlab("Log2 Fold Change") + 
  ylab("-log10 Adjusted P-value")

# Volcano plot using EnhancedVolcano package
# install.packages("BiocManager")
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
# List of top significant miRNAs to label in the volcano plot
top_miRNAs <- c("ENSG00000284190.1", "ENSG00000251856.1", "ENSG00000198983.1", "ENSG00000199158.1", "ENSG00000207691.1")

#  Visualization of Differentially Expressed miRNAs
# Convert DESeq2 results to a data frame
res_df <- as.data.frame(res)

# Add a column to indicate significance based on padj < 0.05
res_df$significant <- ifelse(res_df$padj < 0.05, "Significant", "Non-Significant")

# Load ggplot2 library
library(ggplot2)

# Create the volcano plot with color coding
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +  # Adjust transparency for better visibility
  scale_color_manual(values = c("Significant" = "red", "Non-Significant" = "blue")) +  # Set color for significant vs non-significant
  theme_minimal() +  # Minimal theme for better visualization
  labs(title = "Volcano Plot: Significant vs Non-Significant miRNAs", 
       x = "log2 Fold Change", 
       y = "-log10 Adjusted P-Value")

# Ensembl Gene IDs for miRNA Gene Annotation:

if (!requireNamespace("biomaRt", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)

#then 
# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Example: get gene names for your Ensembl IDs
ensembl_ids <- gsub("\\..*", "", rownames(res_ordered)[1:20])  # strip version like ".1"

annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Merge annotations with DESeq2 results: 

significant_mirnas$ensembl_gene_id <- gsub("\\..*", "", rownames(significant_mirnas))
res_annotated <- merge(significant_mirnas, annotations, by = "ensembl_gene_id")

# To Visualize the top 20 miRNAs across the sample, a  heatmap was generated with the  following command in R,

# Install and load pheatmap 
if (!require(pheatmap)) {
  install.packages("pheatmap")
}
library(pheatmap)

# 1. Filter the significant miRNAs
log2fc_threshold <- 1  # Define the log2 fold change threshold
padj_threshold <- 0.05  # Define the adjusted p-value threshold

# Filter to get the significant miRNAs
significant_mirnas <- res_annotated[res_annotated$padj < padj_threshold & abs(res_annotated$log2FoldChange) > log2fc_threshold, ]

# 2. Create an expression matrix for the significant miRNAs
# Extract only the expression values of the significant miRNAs
significant_mirnas_expression <- counts_data[rownames(counts_data) %in% significant_mirnas$ensembl_gene_id, ]

# Optional: log-transform the counts (common practice for RNA-Seq data)
significant_mirnas_expression_log <- log2(significant_mirnas_expression + 1)  # Adding 1 to avoid log(0)

# 3. Plot the heatmap
pheatmap(significant_mirnas_expression_log, 
         cluster_rows = TRUE,  # Cluster miRNAs (rows)
         cluster_cols = TRUE,  # Cluster samples (columns)
         annotation_col = coldata,  # Add sample annotations (conditions)
         show_rownames = TRUE,  # Show miRNA names
         show_colnames = TRUE,  # Show sample names
         scale = "row",  # Normalize the data by row (miRNA)
         main = "Heatmap of Significant miRNAs")



