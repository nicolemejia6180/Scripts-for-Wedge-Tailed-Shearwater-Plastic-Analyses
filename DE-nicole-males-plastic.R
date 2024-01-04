module load gcc/8.2.0-fasrc01 R/3.6.1-fasrc02
export R_LIBS_USER=$HOME/app/R_3.6:$R_LIBS_USER

library(DESeq2)
library(reshape)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
filePath <- "/n/holyscratch01/edwards_lab/nmejia/"
sampleNames <- c("SAMPLE-PREP0157_FTermREA12618A_E02v1_N014_S13_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_G02v1_N016_S15_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_B03v1_N019_S18_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_D03v1_N021_S20_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_A04v1_N026_S25_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_B04v1_N027_S26_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_D01v1_N005_S4_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_H02v1_N017_S16_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_A02v1_N010_S9_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_B02v1_N011_S10_L004Aligned.toTranscriptome.out.bam")
#save all in a pdf
pdf("WTSH-males-plastic.pdf")
#
countData.list <- sapply(sampleNames, function(x) read.table(file=paste0(filePath, x, ".isoforms.results"), header=T, sep="\t"), simplify=F)
countData.df <- do.call("cbind", countData.list)
colsToKeep <- c(1,grep("expected_count", names(countData.df)))
ct <- countData.df[,colsToKeep]
names(ct) <- c("transcript_id", sampleNames)
ct[,2:11] <- round(ct[,2:11])
sampleMetaData <- data.frame(plastic=rep(c(rep(c("negative"), 6),rep(c("positive"), 4)),1))
#
rownames(sampleMetaData) <- sampleNames
rsem.in <- DESeqDataSetFromMatrix(countData = ct, colData = sampleMetaData, design = ~ plastic, tidy = T)
#exploring and filtering steps remove low low counts close to 0
nrow(rsem.in)
keep <- rowSums(counts(rsem.in)) > 1
#more restricted filtering
keep <- rowSums(counts(rsem.in) >= 10) >= 3
rsem.in <- rsem.in[keep,]
nrow(rsem.in)
#The variance stabilizing transformation and the rlog
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)
#for logaritmic transfromation
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)
#what to chose method? n< 30 rlog, but n<30 vst. but explore both
vsd <- vst(rsem.in, blind = FALSE)
head(assay(vsd), 3)
rld <- rlog(rsem.in, blind = FALSE)
head(assay(rld), 3)
#estimate size facor and plot to compare log2 with normalized count, vsd and rlog methods
rsem.in <- estimateSizeFactors(rsem.in)
#
df <- bind_rows(as_data_frame(log2(counts(rsem.in, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
# 
colnames(df)[1:2] <- c("x", "y")  
#
lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)
#
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  
#calculate distance between samples without experimental desing. Which samples are similar to each other, which are different? Does this fit to the expectation from the experiment’s design?
sampleDists <- dist(t(assay(vsd)))
sampleDists
#Heatmap of sample-to-sample distances using the variance stabilizing transformed values
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$plastic, vsd$sampleNames,sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
#pca to explore groups by experimental desing
plotPCA(vsd, intgroup = c("plastic"))
##########################################finally, Running the differential expression pipeline###############3
rsem.de <- DESeq(rsem.in)
rsem.res <- results(rsem.de)
#specifying the 2 gorups of comparition. HOwever specifiyng or tnot geave the same results
res <- results(rsem.de, contrast=c("plastic","positive","negative"))
#see what have been caluclated on those results
mcols(res, use.names = TRUE)
mcols(rsem.res, use.names = TRUE)
#to prof null hypothesis, give the threshold pval
res.05 <- results(rsem.de, alpha = 0.05)
table(res.05$padj < 0.05)
#raise the log2 fold change threshold
resLFC1 <- results(rsem.de, lfcThreshold=1)
table(resLFC1$padj < 0.1)
#In high-throughput biology, we are careful to not use the p values directly as evidence against the null, but to correct for multiple testing
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
#DESeq2 uses the Benjamini-Hochberg (BH) adjustment (Benjamini and Hochberg 1995) with padj function
sum(res$padj < 0.1, na.rm=TRUE)
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
#with the strongest up-regulation
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
#visualize the counts for a particular gene
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(rsem.de, gene = topGene, intgroup=c("plastic"))
#MA-plot (Dudoit et al. 2002) provides a useful overview for the distribution of the estimated coefficients in the model, e.g. the comparisons of interest, across all genes
#Before making the MA-plot, we use the lfcShrink function to shrink the log2 fold changes for the comparison of dex treated vs untreated samples. There are three types of shrinkage estimators in DESeq2, which are covered in the DESeq2 vignette. Here we specify the apeglm method for shrinking coefficients, which is good for shrinking the noisy LFC estimates while giving low bias LFC estimates for true large differences (Zhu, Ibrahim, and Love 2018)
library("apeglm")
#show the names of the treatmen
resultsNames(rsem.de)
res <- lfcShrink(rsem.de, coef="plastic_positive_vs_negative", type="apeglm")
plotMA(res, ylim = c(-5, 5))
#DESeq2 package uses a Bayesian procedure to moderate (or “shrink”) log2 fold changes from genes with very low counts and highly variable counts.
#in general means to keep ploting the noise removed before
res.noshr <- results(rsem.de, name="plastic_positive_vs_negative")
plotMA(res.noshr, ylim = c(-5, 5))
#label the dots with the "topgene" transcript name
plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
#Another useful diagnostic plot is the histogram of the p values
#Histogram of p values for genes with mean normalized count larger than 1.
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
#gene clustering: 20 genes with the highest variance across samples
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
#anno <- as.data.frame(colData(vsd)[, c("plastic")])
cdata <- colData(vsd)
pheatmap(mat,cluster_rows = T,
    show_rownames = T,
    cluster_cols = T,annotation_col = as.data.frame(cdata),border_color=NA, fontsize = 5,
      fontsize_row = 5, height=100)
#
cdata <- colData(vsd)
pheatmap(mat,
    cluster_rows = T,
    show_rownames = T,
    cluster_cols = F,
    annotation_col = as.data.frame(cdata),border_color=NA, fontsize = 5,
      fontsize_row = 5, height=100)
#
dev.off()
#
setEPS()
postscript("WTHS-MALES-plastic-heatplot2.eps")
hp
dev.off()

