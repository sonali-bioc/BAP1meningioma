prot_mat = readRDS("BAP1_hg38_raw_protein_coding_genes.rds")
sampleinfo = read.delim("updated_sampleinfo.txt", header=T, stringsAsFactors=F)

want = sapply(sampleinfo[,3], function(x) grep(x, colnames(prot_mat)) )
raw_mat = prot_mat[, want]

resdir = "anova_analysis_v2_10_3_2024"

#filter genes with low expression in each group
mat1 = raw_mat[, which(sampleinfo$Group=="BAP1")]
mat2 = raw_mat[, which(sampleinfo$Group=="H3K27me3 loss")]
mat3= raw_mat[, which(sampleinfo$Group=="NF2")]

ch1 = apply(mat1, 1 , function(x) length(which(x!=0)))
ch2= apply(mat2, 1 , function(x) length(which(x!=0)))
ch3 = apply(mat3, 1 , function(x) length(which(x!=0)))

c1 = which(ch1 >=3)
c2 = which(ch2 >=3)
c3 = which(ch3 >=3)

keep = intersect(intersect(c1, c2) , c3)
length(keep)
raw_mat= raw_mat[keep, ]

# deseq2 analysis

library(DESeq2)
coldata = data.frame(sampleGroup = sampleinfo$Group, 
                     sampleName = colnames(raw_mat))
rownames(coldata) = colnames(raw_mat)
dds <- DESeqDataSetFromMatrix(countData =round (raw_mat), 
                              colData = sampleinfo, design = ~ Group)
dds = DESeq(dds)
lfc = log2(1.5)

# bap1 vs nf2            
res1 <- results(dds, alpha = 0.05,  contrast = c("Group", "BAP1", "NF2"))

resdf = as.data.frame(res1)
raw_counts = raw_mat[rownames(resdf), ]
library(edgeR)
cpm_counts = cpm(raw_counts)
colnames(raw_counts) = paste0("raw_", colnames(raw_counts))
colnames(cpm_counts) = paste0("cpm_", colnames(cpm_counts))

avg_mat = cbind( rowMeans(cpm_counts[, which(sampleinfo$Group=="BAP1")]), 
                 rowMeans(cpm_counts[, which(sampleinfo$Group=="NF2")]) ) 
colnames(avg_mat) = c("CPM_Avg_BAP1", "CPM_Avg_NF2")
resdf = cbind(gene = rownames(resdf), resdf[, c("log2FoldChange", "pvalue", "padj")] , 
              avg_mat)

up_genes = resdf[ which(resdf$log2FoldChange > lfc &  resdf$padj < 0.05), ]
down_genes = resdf[which(resdf$log2FoldChange < -lfc &  resdf$padj < 0.05), ]  

up_genes = up_genes[order(up_genes$log2FoldChange, decreasing=T), ]
down_genes = down_genes[order(down_genes$log2FoldChange), ]
library(writexl)
lst = list(all_genes = resdf, up_genes = up_genes, down_genes = down_genes)
fname =paste0(c("DESeq2_analysis_BAP1_vs_NF2"), collapse ="")
write_xlsx(lst, file.path( paste0(fname  , ".xlsx")))

up_genes1 = up_genes
down_genes1 = down_genes
            
# bap1 vs h3k27me3 loss

res1 <- results(dds, alpha = 0.05,  contrast = c("Group", "BAP1", "H3K27me3 loss"))
resdf = as.data.frame(res1)
raw_counts = raw_mat[rownames(resdf), ]
library(edgeR)
cpm_counts = cpm(raw_counts)
colnames(raw_counts) = paste0("raw_", colnames(raw_counts))
colnames(cpm_counts) = paste0("cpm_", colnames(cpm_counts))

avg_mat = cbind( rowMeans(cpm_counts[, which(sampleinfo$Group=="BAP1")]), 
                 rowMeans(cpm_counts[, which(sampleinfo$Group=="H3K27me3 loss")]) ) 
colnames(avg_mat) = c("CPM_Avg_BAP1", "CPM_Avg_H3K27me3 loss")
resdf = cbind(gene = rownames(resdf), resdf[, c("log2FoldChange", "pvalue", "padj")] , 
              avg_mat)

up_genes = resdf[ which(resdf$log2FoldChange > lfc &  resdf$padj < 0.05), ]
down_genes = resdf[which(resdf$log2FoldChange < -lfc &  resdf$padj < 0.05), ]  

up_genes = up_genes[order(up_genes$log2FoldChange, decreasing=T), ]
down_genes = down_genes[order(down_genes$log2FoldChange), ]
library(writexl)
lst = list(all_genes = resdf, up_genes = up_genes, down_genes = down_genes)
fname =paste0(c("DESeq2_analysis_BAP1_vs_H3K27me3 loss"), collapse ="")
write_xlsx(lst, file.path( paste0(fname  , ".xlsx")))

up_genes2 = up_genes
down_genes2 = down_genes

# h3k27me3loss vs nf2
res1 <- results(dds, alpha = 0.05,  contrast = c("Group",  "H3K27me3 loss", "NF2"))
resdf = as.data.frame(res1)
raw_counts = raw_mat[rownames(resdf), ]
library(edgeR)
cpm_counts = cpm(raw_counts)
colnames(raw_counts) = paste0("raw_", colnames(raw_counts))
colnames(cpm_counts) = paste0("cpm_", colnames(cpm_counts))

avg_mat = cbind( rowMeans(cpm_counts[, which(sampleinfo$Group=="NF2")]), 
                 rowMeans(cpm_counts[, which(sampleinfo$Group=="H3K27me3 loss")]) ) 
colnames(avg_mat) = c("CPM_Avg_nf2", "CPM_Avg_H3K27me3 loss")
resdf = cbind(gene = rownames(resdf), resdf[, c("log2FoldChange", "pvalue", "padj")] , 
              avg_mat)

up_genes = resdf[ which(resdf$log2FoldChange > lfc &  resdf$padj < 0.05), ]
down_genes = resdf[which(resdf$log2FoldChange < -lfc &  resdf$padj < 0.05), ]  

up_genes = up_genes[order(up_genes$log2FoldChange, decreasing=T), ]
down_genes = down_genes[order(down_genes$log2FoldChange), ]
library(writexl)
lst = list(all_genes = resdf, up_genes = up_genes, down_genes = down_genes)
fname =paste0(c("DESeq2_analysis_H3K27me3loss_vs_NF2"), collapse ="")
write_xlsx(lst, file.path(paste0(fname  , ".xlsx")))

up_genes3= up_genes
down_genes3 = down_genes

## Enrichment analysis ussing up and down-regulated genes


library(enrichR)
my_enrichment_function = function(goi, title, res_folder){
    dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", 
             "GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022")
    enriched = enrichr( goi, dbs)
    enriched= lapply(enriched, function(z){
        rm_cols = match(c("Old.P.value", "Old.Adjusted.P.value"), colnames(z))
        if(length(rm_cols)!=0){
            z = z[, -c(rm_cols)]
        }
        z = z[which(z$Adjusted.P.value < 0.05), ]
        
    })
    
    write_xlsx(x = enriched, path = file.path(res_folder, paste0(title)) )
    enriched
}

a1  = my_enrichment_function(up_genes1$gene, 
                             title = paste0("top_",x,"_genesEnrichr_up_reg_BAP1_vs_NF2.xlsx"), resdir)
b1  = my_enrichment_function(down_genes1$gene, 
                             title = paste0("top_",x,"_genesEnrichr_down_reg_BAP1_vs_NF2.xlsx"), resdir)

a2  = my_enrichment_function(up_genes2$gene, 
                             title =  paste0("top_",x,"_genesEnrichr_up_reg_BAP1_vs_H3K27me3loss.xlsx"), resdir)
b2  = my_enrichment_function(down_genes2$gene, 
                             title =  paste0("top_",x,"_genesEnrichr_down_reg_BAP1_vs_H3K27me3loss.xlsx"), resdir)

a3  = my_enrichment_function(up_genes3$gene, 
                             title =  paste0("top_",x,"_genesEnrichr_up_reg_H3K27me3loss_vs_NF2.xlsx"), resdir)  
b3  = my_enrichment_function(down_genes3$gene, 
                             title =  paste0("Enrichr_down_reg_H3K27me3loss_vs_NF2.xlsx"), resdir)


            
