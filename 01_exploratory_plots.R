
setwd("~/HollandLabShared/Sonali/Frank/BAP1_signaling")
prot_mat = readRDS("BAP1_hg38_raw_protein_coding_genes.rds")

sampleinfo = read.delim("updated_sampleinfo.txt", header=T, stringsAsFactors=F)

library(ggplot2)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = raw_mat, 
                              colData = sampleinfo, design = ~ Group)
norm_data <-vst(data.matrix(raw_mat))

sampleDists <- dist( t( norm_data ) )
sampleDistMatrix <- as.matrix( sampleDists )
hc = hclust(sampleDists)

sampleName = colnames(raw_mat)

mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(sampleinfo))
m1 <- ggplot(mds, aes(X1,X2,color=Group)) +
    geom_point(size=4)+ theme_bw() +
    ggtitle(paste0("MDS plot : ")) +
    theme(plot.title = element_text(lineheight=.8, face="bold")) +
    geom_text(aes(label=sampleName),hjust="inward", vjust=2, size=4)

pc= prcomp(t(norm_data))
pc_data1 = data.frame(PC1=pc$x[,1], PC2=pc$x[,2],
                      Group=(sampleinfo[,"Group"]) )
percentVar <- (pc$sdev^2 / sum( pc$sdev^2 ) )*100
percentVar= round(percentVar[1:2], 2)
p1 = ggplot(pc_data1, aes(PC1, PC2, color=Group)) +
    geom_point(size=4) +theme_bw() +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))+
    ggtitle(paste0("PCA plot : ")) +
    geom_text(aes(label=sampleName),hjust="inward", vjust=2, size=4)+
    theme(plot.title = element_text(lineheight=.8, face="bold"))


pdf("exploratory_plots.pdf", width =10)
plot(hc)
print(m1)
print(p1)
dev.off()
