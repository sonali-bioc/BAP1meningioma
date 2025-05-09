
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


# boxplots 

lst = list(c("ERBB2", "ERBB3", "NRG3"),
            c("FGFR2", "FGF1", "FGF14", "FGF17", "FGF18"),  
            c("NTRK1", "NGF", "BDNF", "NTF3"), 
            c("WNT3A", "WNT10A", "WNT10B", "WNT11", "WNT16"),# "WNT8B", 
           c("NOTCH1", "NOTCH2"))


library(reshape2)
pdf("boxplot_gene_families_10_24_2024.pdf")
l1 = lapply(lst, function(mygoi){
    mydf = norm_data[mygoi, ]
    ggdf2 = melt(mydf)
    ggdf2$group = unlist(lapply(sampleinfo$Group, function(x) rep(x, length(mygoi))) )
    ggdf2$Var1  = factor(ggdf2$Var1, levels = mygoi)
    
    p1 = ggplot(ggdf2, aes(x=Var1, y=value,  fill=group)) +
        geom_boxplot() +
        ylab("VST counts per gene") + xlab("") +
        ggtitle("") + theme_bw() + 
        theme(
            legend.text=element_text(size=12),
            legend.position ='bottom',
            legend.justification = 'left',
            legend.spacing.y = unit(0.5, 'cm'),
            axis.title=element_text(size=12),
            axis.text=element_text(size=12), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(p1)
})
dev.off()


goi_lst = list(c("ERBB2", "ERBB3", "NRG3"),
           c("FGFR2", "FGF1", "FGF14", "FGF17", "FGF18"),  #"FGF7", 
           c("NTRK1", "NGF", "BDNF", "NTF3"), 
           c("PDGFA", "PDGFD"), 
           c("TGFBR3", "TGFA"), #"TGFB3", 
           c("RET", "GFRA2", "GFRA4", "GDNF"), 
           c("VEGFA", "VEGFC"), 
           c("WNT3A", "WNT10A", "WNT10B", "WNT11", "WNT16"),# "WNT8B", 
           c("NOTCH1", "NOTCH2"))
                         
