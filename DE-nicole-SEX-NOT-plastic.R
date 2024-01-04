library(DESeq2)
library(reshape)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(pheatmap)
#
#save all in a pdf
pdf("WTSH-SEX-NOT-plastic.pdf")
#
filePath <- "/n/holyscratch01/edwards_lab/nmejia/"
sampleNames <- c("SAMPLE-PREP0157_FTermREA12618A_A01v1_N001_S1_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_C01v1_N003_S3_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_C02v1_N012_S11_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_D02v1_N013_S12_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_H01v1_N009_S8_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_F03v1_N023_S22_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_B01v1_N002_S2_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_E01v1_N006_S5_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_G01v1_N008_S7_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_A03v1_N018_S17_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_E03v1_N022_S21_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_D04v1_N029_S28_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_E02v1_N014_S13_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_G02v1_N016_S15_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_B03v1_N019_S18_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_D03v1_N021_S20_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_A04v1_N026_S25_L004Aligned.toTranscriptome.out.bam","SAMPLE-PREP0157_FTermREA12618A_B04v1_N027_S26_L004Aligned.toTranscriptome.out.bam")
##prepare matrix
countData.list <- sapply(sampleNames, function(x) read.table(file=paste0(filePath, x, ".isoforms.results"), header=T, sep="\t"), simplify=F)
countData.df <- do.call("cbind", countData.list)
colsToKeep <- c(1,grep("expected_count", names(countData.df)))
ct <- countData.df[,colsToKeep]
names(ct) <- c("transcript_id", sampleNames)
ct[,2:19] <- round(ct[,2:19])
sampleMetaData <- data.frame(sex=rep(c(rep(c("females"), 12),rep(c("males"), 6)),1))
#
rownames(sampleMetaData) <- sampleNames
rsem.in <- DESeqDataSetFromMatrix(countData = ct, colData = sampleMetaData, design = ~ sex, tidy = T)
#exploring and filtering steps remove low low counts close to 0
nrow(rsem.in)
keep <- rowSums(counts(rsem.in)) > 1
#more restricted filtering
#keep <- rowSums(counts(rsem.in) >= 10) >= 3
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
rownames(sampleDistMatrix) <- paste( vsd$sex, vsd$sampleNames,sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
#pca to explore groups by experimental desing
plotPCA(vsd, intgroup = c("sex"))
##########################################finally, Running the differential expression pipeline###############3
rsem.de <- DESeq(rsem.in)
rsem.res <- results(rsem.de)
#specifying the 2 gorups of comparition. HOwever specifiyng or tnot geave the same results
res <- results(rsem.de, contrast=c("sex","females","males"))
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
plotCounts(rsem.de, gene = topGene, intgroup=c("sex"))
#MA-plot (Dudoit et al. 2002) provides a useful overview for the distribution of the estimated coefficients in the model.
#Before making the MA-plot, we use the lfcShrink function.
#
library("apeglm")
#
#show the names of the treatmen
resultsNames(rsem.de)
res <- lfcShrink(rsem.de, coef="sex_males_vs_females", type="apeglm")
plotMA(res, ylim = c(-5, 5))
#DESeq2 package uses a Bayesian procedure to moderate (or “shrink”) log2 fold changes from genes with very low counts and highly variable counts.
#in general means to keep ploting the noise removed before
res.noshr <- results(rsem.de, name="sex_males_vs_females")
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
#anno <- as.data.frame(colData(vsd)[, c("sex")])
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
