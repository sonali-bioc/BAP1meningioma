## Read in gene lists 

gene_list = read.delim("GENE_LISTS_bap1_prc_target_genes.txt", header=T, stringsAsFactors = FALSE)
bap1_targets = toupper(gene_list[,1])
bap1_targets = bap1_targets[-c(which(bap1_targets==""))]
SUZ12.targets = toupper(gene_list[,2])

library(readxl)
res = as.data.frame(read_xlsx(file.path(resdir, "DESeq2_analysis_BAP1_vs_NF2.xlsx"), sheet = 1))

colnames(res)[2] = "logFc"
res$log10padj = -log10(res$padj)

res$delabel = NA_character_
res$delabel[which( res$padj<.05 & abs(res$logFc)>lfc) ] = "sig"
up1 = res[ which(res$logFc > lfc &  res$padj < 0.05), ]
down1 = res[ which(res$logFc < -lfc & res$padj < 0.05), ]

# direct bap1 targets among differentially expressed genes in BAP1-altered meningioma
pdf("Volcano_plot_bap1.pdf", width = 10)
chosen_goi = intersect( bap1_targets, res$gene)
with(res, plot(logFc, log10padj, pch=20, col ="grey80",main=paste0("Volcano plot : BAP1 targets") ))
with(subset(res, res$padj<.05 & logFc>lfc_cutoff), points(logFc, log10padj, pch=20, col="lightsalmon"))
with(subset(res, res$padj<.05 & logFc < -lfc_cutoff), points(logFc, log10padj, pch=20, col="lightgreen"))
abline(v = lfc_cutoff , col = "black", lty = 2)
abline(v = -lfc_cutoff, col = "black", lty = 2)
abline(h = -log10(0.05), col = "black", lty = 2)
with(res[which(res[,1] %in% chosen_goi), ], points(logFc, log10padj, pch=20,  col="black"))
dev.off()

# direct SUZ12 targets among differentially expressed genes in BAP1-altered meningioma
pdf("Volcano_plot_suz12.pdf", width = 10)
chosen_goi = intersect( SUZ12.targets, res$gene)
with(res, plot(logFc, log10padj, pch=20, col ="grey80",main=paste0("Volcano plot : SUZ12.targets") ))
with(subset(res, res$padj<.05 & logFc>lfc_cutoff), points(logFc, log10padj, pch=20, col="lightsalmon"))
with(subset(res, res$padj<.05 & logFc < -lfc_cutoff), points(logFc, log10padj, pch=20, col="lightgreen"))
abline(v = lfc_cutoff , col = "black", lty = 2)
abline(v = -lfc_cutoff, col = "black", lty = 2)
abline(h = -log10(0.05), col = "black", lty = 2)
with(res[which(res[,1] %in% chosen_goi), ], points(logFc, log10padj, pch=20,  col="black"))
dev.off()


