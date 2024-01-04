module load gcc/8.2.0-fasrc01 R/3.6.1-fasrc02
export R_LIBS_USER=$HOME/app/R_3.6:$R_LIBS_USER

BiocManager::install("DESeq2")
BiocManager::install("reshape")
BiocManager::install("ggplot2")
BiocManager::install("ggrepel")
BiocManager::install("DEGreport")
BiocManager::install("RColorBrewer")
BiocManager::install("pheatmap")


library("DESeq2")
library("reshape")
library("ggplot2")
library("ggrepel")
library("DEGreport")
library("RColorBrewer")
library("pheatmap")

#black gonads
filePath <- "/n/holyscratch01/edwards_lab/nmejia/"
sampleNames <- c("SAMPLE-PREP0157_FTermREA12618A_C01v1_N003_S3_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_E01v1_N006_S5_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_G01v1_N008_S7_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_H01v1_N009_S8_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_C02v1_N012_S11_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_D02v1_N013_S12_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_A03v1_N018_S17_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_E03v1_N022_S21_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_F03v1_N023_S22_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_D04v1_N029_S28_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_B01v1_N002_S2_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_F01v1_N007_S6_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_F02v1_N015_S14_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_C03v1_N020_S19_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_G03v1_N024_S23_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_H03v1_N025_S24_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_C04v1_N028_S27_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_E02v1_N014_S13_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_G02v1_N016_S15_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_B03v1_N019_S18_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_D03v1_N021_S20_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_A04v1_N026_S25_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_B04v1_N027_S26_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_D01v1_N005_S4_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_A02v1_N010_S9_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_B02v1_N011_S10_L004Aligned.toTranscriptome.out.bam",
                 "SAMPLE-PREP0157_FTermREA12618A_H02v1_N017_S16_L004Aligned.toTranscriptome.out.bam")
                 
)
#
countData.list <- sapply(sampleNames, function(x) read.table(file=paste0(filePath, x, ".isoforms.results"), header=T, sep="\t"), simplify=F)
countData.df <- do.call("cbind", countData.list)
colsToKeep <- c(1,grep("expected_count", names(countData.df)))
ct <- countData.df[,colsToKeep]
names(ct) <- c("transcript_id", sampleNames)
ct[,2:28] <- round(ct[,2:28])
sampleMetaData <- data.frame(sex=rep(c(rep(c("females"), 17),rep(c("males"), 10)),1))
#
rownames(sampleMetaData) <- sampleNames
rsem.in <- DESeqDataSetFromMatrix(countData = ct, colData = sampleMetaData, design = ~ sex, tidy = T)
rsem.de <- DESeq(rsem.in)
rsem.res <- results(rsem.de)
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
threshold <- rsem.res$padj < padj.cutoff & abs(rsem.res$log2FoldChange) > lfc.cutoff
length(which(threshold))
rsem.res$threshold <- threshold                
sigOE <- data.frame(subset(rsem.res, threshold==TRUE))
sigOE_ordered <- sigOE[order(sigOE$padj), ]
### Set a color palette
heat.colors <- brewer.pal(6, "YlOrRd")
### Run pheatmap
normalized_counts <- counts(rsem.de, normalized=T)
norm_OEsig <- normalized_counts[rownames(sigOE),]
heat.colors <- brewer.pal(6, "YlOrRd")
pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=F,border_color=NA, fontsize = 10, scale="row",
      fontsize_row = 10, height=20)
#
hp <- pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=F,border_color=NA, fontsize = 10, scale="row",
      fontsize_row = 10, height=20)
setEPS()
postscript("WTHS-SEX-heatplot.eps")
hp
dev.off()
######################
rsem.res <-results(rsem.de, alpha = 0.05, lfcThreshold = 0.58)
summary(rsem.res)
outlier <- subset(rsem.res, padj < 0.05)
write.table(file="WTSH-outliers-sex.txt", norm_OEsig, quote=FALSE, row=TRUE)
#volcano plot#########
resOE_df_ordered <- norm_OEsig[order(norm_OEsig$padj), ] 

resOE_df_ordered$genelabels <- rownames(resOE_df_ordered) %in% rownames(resOE_df_ordered[1:10,])

View(resOE_df_ordered)
resOE_df <- data.frame(rsem.res)

View(resOE_df)
# Volcano plot
pp <- ggplot(resOE_df) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) + theme_bw() +
  ggtitle("Differential geen expression between sex") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
#
setEPS()
postscript("Volcanoplot-Differential geen expression between sex.eps")
pp
dev.off()
