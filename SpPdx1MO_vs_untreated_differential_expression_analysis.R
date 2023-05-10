library("DESeq2")
library("IHW")
library("tximport")
library("tidyverse")
library("readr")
library("ggplot2")
library("scales")
library("ggpie")
setwd("~/Work/Collaborations/downstream_pdx1/RNA-seq")
set.seed(255)
dirs_sp48 <- list.files(file.path("./sp48lox_wt_1Lox/"))
quant_files_sp48 <- paste0("./sp48lox_wt_1Lox/", dirs_sp48, "/quant.sf")
renamefiles_sp48 <- c("sp48lox1", "sp48lox2", "sp48lox3", "sp48wt1", "sp48wt2", "sp48wt3")
names(quant_files_sp48) <- renamefiles_sp48
file.exists(quant_files_sp48) #true
transcript2whl <- read_csv(file.path("~/Work/Collaborations/downstream_pdx1/RNA-seq", "whl_tx2gene_oneLox.csv"))
tx_sp48 <- tximport(quant_files_sp48, type = "salmon", countsFromAbundance = "no", txOut = FALSE, tx2gene = transcript2whl)

sampleTable_sp48 <- data.frame(condition = factor(rep(c("Sp-Pdx1 MO", "Untreated"), each = 3)))
rownames(sampleTable_sp48) <- colnames(tx_sp48$counts)


dds <- DESeqDataSetFromTximport(tx_sp48, sampleTable_sp48, ~condition) 

dds$condition <- relevel(dds$condition, ref = "Untreated")
dds <- estimateSizeFactors(dds)
deseq_sp48 <- DESeq(dds, fitType = "local")
#plotDispEsts(deseq_sp48)
rlog_sp48 <- rlog(deseq_sp48)



names <- read.delim("./WHL_names.tsv", header=FALSE, stringsAsFactors=FALSE)
names <- names %>% rename("gene" = "V1", "name"="V2")
sp48_hindgut_cells_whls <- scan("./sp48_hindgut_cells_whls.txt", character(), quote = "")
norm_counts <- as.data.frame(counts(dds, normalized= TRUE))
norm_counts$gene <- row.names(norm_counts)

res <- results(deseq_sp48, contrast=c("condition", "Sp-Pdx1 MO", "Untreated"), filterFun = ihw)
write.table(subset(res, padj < 0.05), file = "Sp48_lox_wt_gc_seq_bias_oneLox.tsv", quote = F, sep = "\t")
res_genes <- as.data.frame(subset(res,padj < 0.05))
res_genes$gene <- row.names(res_genes)
res_genes <- inner_join(res_genes, names, by ="gene")

hindgut_res_genes <- res_genes[res_genes$gene %in% sp48_hindgut_cells_whls,]
write.table(hindgut_res_genes, file = "hindgut_Sp48_lox_wt_gc_seq_bias_oneLox.tsv", quote = F, sep = "\t")


write.table(res, file = "Sp48_lox_wt_full.tsv", quote = F, sep = "\t")




