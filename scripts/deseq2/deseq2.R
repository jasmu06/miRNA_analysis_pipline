
cat > scripts/deseq2.R << 'EOF'
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]

library(DESeq2)
counts <- read.table(counts_file, header=TRUE, row.names=1)
metadata <- read.csv(metadata_file, row.names=1)

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Tumor","Control"))

dir.create(output_dir, recursive=TRUE)
write.csv(as.data.frame(res), file.path(output_dir, "DEG_results.csv"))

pdf(file.path(output_dir, "pca_plot.pdf"))
plotPCA(vst(dds), intgroup="condition")
dev.off()
EOF

chmod +x scripts/deseq2.R

q

