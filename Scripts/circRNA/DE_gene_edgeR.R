rm(list=ls())
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(optparse))

setwd("D:/Fan Lab/Research-data/AV-RNA-seq/CircRNA/DE circRNA/mm10/quant_2_3/av_wt")
# Load data
lib_mtx <- read.csv(file = "library_info.csv", row.names = 1)
gene_mtx <- read.csv("gene_count_matrix.csv", row.names = 1, check.names=FALSE)
gene_mtx <- gene_mtx[,rownames(lib_mtx)]
bsj_mtx <- read.csv("circRNA_bsj.csv", row.names = 1, check.names=FALSE)

gene_DGE <- DGEList(counts = gene_mtx, group = lib_mtx$Group)
gene_idx <- filterByExpr(gene_DGE)
gene_DGE <- gene_DGE[gene_idx, keep.lib.sizes=FALSE]
gene_DGE <- calcNormFactors(gene_DGE)


if ("Subject" %in% colnames(lib_mtx)) {
  subject <- factor(lib_mtx$Subject)
  treat <- factor(lib_mtx$Group, levels=c("C", "T"))
  
  design <- model.matrix(~subject + treat)
} else {
  treat <- factor(lib_mtx$Group, levels=c("C", "T"))
  design <- model.matrix(~treat)
}

design <- model.matrix(~factor(lib_mtx$Group))
gene_DGE <- estimateDisp(gene_DGE, design, robust = TRUE)
gene_fit <- glmFit(gene_DGE, design)
gene_lrt <- glmLRT(gene_fit)
#
gene_df <- gene_lrt$table
gene_order <- order(gene_lrt$table$PValue)
gene_df$DE <- decideTestsDGE(gene_lrt)
gene_df <- gene_df[gene_order, ]
gene_df$FDR <- p.adjust(gene_df$PValue, method="fdr")

write.table(gene_df, "D:/Fan Lab/Research-data/AV-RNA-seq/CircRNA/DE circRNA/mm10/quant_2_3/av_wt/AV_WT_gene.tsv", 
            sep = "\t", quote = F, row.names = T)

circ_DGE <- DGEList(counts = bsj_mtx,
                    group = lib_mtx$Group,
                    lib.size = gene_DGE$samples[, "lib.size"],
                    norm.factors = gene_DGE$samples[, "norm.factors"])

# circ_idx <- filterByExpr(circ_DGE, min.count=)
# circ_DGE <- circ_DGE[circ_idx, , keep.lib.sizes=TRUE]
# head(circ_df[rowSums(bsj_mtx >= 2) >= nrow(lib_mtx) / 2,])

circ_DGE <- estimateDisp(circ_DGE, design, robust = TRUE)
circ_fit <- glmFit(circ_DGE, design)
circ_lrt <- glmLRT(circ_fit)

circ_df <- circ_lrt$table
circ_order <- order(circ_lrt$table$PValue)
circ_df$DE <- as.vector(decideTestsDGE(circ_lrt))
circ_df <- circ_df[circ_order, ]
circ_df$FDR <- p.adjust(circ_df$PValue, method="fdr")
write.csv(circ_df, file=opt$out, quote = FALSE)