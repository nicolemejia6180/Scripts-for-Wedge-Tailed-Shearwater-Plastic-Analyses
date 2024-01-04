#prepare gene list
for i in *txt; do awk '{print $1}' $i >> $i-transcripts.txt;done
GTF=/n/holyscratch01/edwards_lab/nmejia/ref/GCA_013401115.1_ASM1340111v1_genomic.gff
for i in *transcripts.txt; do cat $GTF | grep -w -f $i >> $i-GENES-transcript-outliers.txt;done
for i in *-GENES-transcript-outliers.txt; do awk '{print $9}' $i >> $i-GENES-list-outliers.txt;done
#after this parsing,need to remove some characters from the file before allenricher
#to get just the gene ID
grep -Eo "gene=[^;]*" females-plastic-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt
for i in *-GENES-list-outliers.txt; do grep -Eo "gene=[^;]*" $i >> $i-JUST-GENES-list-outliers.txt;done
#after this before runing allenricher need to remove the "gene="
for i in *-JUST-GENES-list-outliers.txt; do sed -i 's/gene=//g' $i;done
#to get the uniq genes
for i in *-JUST-GENES-list-outliers.txt; do uniq -c $i > $i-uniq.txt;done
############divide up and dowregulated genes###########################3
cat males-plastic-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -f males-plastic-upreg.txt >> males-plastic-upreg-GENE-list.txt
cat males-plastic-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -v -f males-plastic-upreg.txt >> males-plastic-downreg-GENE-list.txt
#
cat sex-All-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -f sex-all-downreg.txt >> sex-all-downreg-GENE-list.txt
cat sex-All-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -v -f sex-all-downreg.txt >> sex-all-upreg-GENE-list.txt
#
cat SEX-noplastic-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -f sex-noPlastic-upreg.txt >> sex-noPlastic-upreg-GENE-list.txt
cat SEX-noplastic-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -v -f sex-noPlastic-upreg.txt >> sex-noPlastic-Downreg-GENE-list.txt
#
cat SEX-Plastic-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -f sex-plastic-upreg.txt >> sex-plastic-upreg-GENE-list.txt
cat SEX-Plastic-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -v -f sex-plastic-upreg.txt >> sex-plastic-downreg-GENE-list.txt
#
cat weight-3factor-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -f weight-all-upreg.txt >> weight-all-upreg-GENE-list.txt
cat weight-3factor-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -v -f weight-all-upreg.txt >> weight-all-downreg-GENE-list.txt
#
cat weight-3factor-females-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -f weight-females-upreg.txt >> weight-females-upreg-GENE-list.txt
cat weight-3factor-females-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt | grep -w -v -f weight-females-upreg.txt >> weight-females-downreg-GENE-list.txt
#
cat weight-2factor-males-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt |grep -w -f weight-males-upreg.txt >> weight-males-upreg-GENE-list.txt
cat weight-2factor-males-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt |grep -w -v -f weight-males-upreg.txt >> weight-males-downreg-GENE-list.txt
#
GTF=/n/holyscratch01/edwards_lab/nmejia/ref/GCA_013401115.1_ASM1340111v1_genomic.gff
cat $GTF | grep -w -f general-free-upreg.txt >> general-free-upreg-GENE-list.txt
#get just the gene list
for i in *GENE-list.txt; do grep -Eo "gene=[^;]*" $i >> $i-JUST-GENES.txt;done
#after this before runing allenricher need to remove the "gene="
for i in *-JUST-GENES.txt; do sed -i 's/gene=//g' $i;done
#
cat $GTF | grep  -w -f weight-females-upreg.txt >> weight-females-upreg-GENE-list.txt
grep -Eo "gene=[^;]*" weight-females-upreg-GENE-list.txt >> weight-females-upreg-GENE-list.txt-JUST-GENES.txt
sed -i 's/gene=//g' weight-females-upreg-GENE-list.txt-JUST-GENES.txt
################Enrichment analysis################################
#need to create an envirionment to run this program AllEnricher
conda activate AllEnricher
module load gcc/8.2.0-fasrc01 R/3.6.1-fasrc02
export R_LIBS_USER=$HOME/app/R_3.6:$R_LIBS_USER
#
#
AE=/n/holylfs04/LABS/edwards_lab/Lab/ftermignoni/BWA-chr/bamOUTsex/BAYESCAN/40CHR/AllEnricher-master/AllEnricher
#run AllEnrichment to the gene list you create before with bedtools
perl $AE -l SEX-noplastic-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt -s gga -i GO -o ./allenricher2 -r $(command -v Rscript)
###
for i in *-JUST-GENES-list-outliers.txt; do perl $AE -l $i -s gga -i GO+KEGG -c 0.1 -n 1 -o ./allenricher -r $(command -v Rscript);done
#########PLOT RESULTS###
library(DESeq2)
library(airway)
library(gprofiler2)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
#
geneslist <- read.table("SEX-noplastic-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt")
geneslist <- read.table("males-plastic-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt")
geneslist <- read.table("females-plastic-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt")
geneslist <- read.table("weight-3factor-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt")
geneslist <- read.table("sex-All-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt")
#
geneslist <- read.table("weight-3factor-females-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt")
geneslist <- read.table("weight-2factor-males-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt")
#gene list of experiments
weight-3factor-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
weight-3factor-females-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
weight-3factor-females-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt
weight-2factor-males-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
weight-2factor-males-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
SEX-Plastic-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
SEX-Plastic-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
SEX-noplastic-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
SEX-noplastic-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
sex-All-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
sex-All-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
males-plastic-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
males-plastic-outliers-0.1.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
females-plastic-outliers.txt-transcripts.txt-GENES-transcript-outliers.txt-GENES-list-outliers.txt-JUST-GENES-list-outliers.txt
#
geneslist2 <- unlist(geneslist)
gp_up = gost(geneslist2, organism = "hsapiens")
gp_up = gost(geneslist2, organism = "mmusculus")
gp_up = gost(geneslist2, organism = "tguttata")
gp_up = gost(geneslist2, organism = "ggallus")
#for "tguttata" nothing the same for ggallus. just if we ask not to only take significants
#
gp_up = gost(geneslist2, organism = "ggallus", significant=F)
#make a plot and save it
pdf("Terms-source-weight-3factor-all-outliers-0.1_significant")
gostplot(gp_up, interactive =F)
dev.off()
head(gp_up$result)
#save table results
gp_up = gost(geneslist2, organism = "ggallus", significant=F)
R <- gp_up$result
write.table(file="Terms-source-weight-3factor-females-outliers-0.1_non-significant-GALLUS.txt", R, quote=FALSE, row=TRUE)
#########################################
mouse_genes = gorth(geneslist2, source_organism = "ggallus",target_organism = "mmusculus")
#
# modify the g:Profiler data 
gp_mod = gp_up$result[,c("query", "source", "term_id","term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size")]
gp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
#
gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
names(gp_mod) = c("term_name", "term_id","p_value", "query_size", "term_size", "effective_domain_size", "GeneRatio", "BgRatio")
gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
row.names(gp_mod) = gp_mod$ID
# define as compareClusterResult object
gp_mod_cluster = new("compareClusterResult", compareClusterResult = gp_mod)
# defi ne as enrichResult object
gp_mod_enrich  = new("enrichResult", result = gp_mod)
enrichplot::dotplot(gp_mod_cluster)
barplot(gp_mod_enrich, showCategory =40, font.size =16) +  ggplot2::facet_grid(~Cluster) +  ggplot2::ylab("Intersection size")
#
###################separate down & up regulated gene lists#####################
general-free-upreg-GENE-list.txt-JUST-GENES.txt     sex-plastic-upreg-GENE-list.txt-JUST-GENES.txt
males-plastic-downreg-GENE-list.txt-JUST-GENES.txt  weight-all-downreg-GENE-list.txt-JUST-GENES.txt
males-plastic-upreg-GENE-list.txt-JUST-GENES.txt    weight-all-upreg-GENE-list.txt-JUST-GENES.txt
sex-all-downreg-GENE-list.txt-JUST-GENES.txt	    weight-females-downreg-GENE-list.txt-JUST-GENES.txt
sex-all-upreg-GENE-list.txt-JUST-GENES.txt	    weight-females-upreg-GENE-list.txt-JUST-GENES.txt
sex-noPlastic-Downreg-GENE-list.txt-JUST-GENES.txt  weight-males-downreg-GENE-list.txt-JUST-GENES.txt
sex-noPlastic-upreg-GENE-list.txt-JUST-GENES.txt    weight-males-upreg-GENE-list.txt-JUST-GENES.txt
sex-plastic-downreg-GENE-list.txt-JUST-GENES.txt
#
geneslist <- read.table("general-free-upreg-GENE-list.txt-JUST-GENES.txt")
#
geneslist2 <- unlist(geneslist)
#
gp_up = gost(geneslist2, organism = "ggallus")
gp_up = gost(geneslist2, organism = "tguttata")
gp_up = gost(geneslist2, organism = "mmusculus")
gp_up = gost(geneslist2, organism = "hsapiens")
#
gp_up = gost(geneslist2, organism = "ggallus", significant=F)
#make a plot of terms and save it
pdf("Terms-general-upreg-GALLOS-NONSIG.pdf")
pdf("Terms-general-upreg-MUSCULUS.pdf")
gostplot(gp_up, interactive =F)
dev.off()
head(gp_up$result)
#save table results for plot later
df_Place2 = data.frame(lapply(gp_up$result, as.character), stringsAsFactors=FALSE)
write.table(file="Terms-general-upreg-MUSCULUS.txt", df_Place2, quote=FALSE, row=TRUE)
#
df_Place2 = data.frame(lapply(gp_up$result, as.character), stringsAsFactors=FALSE)
write.table(file="Terms-general-upreg-GALLOS-NONSIG.txt", df_Place2, quote=FALSE, row=TRUE)
#########################################################################################
geneslist <- read.table("males-plastic-downreg-GENE-list.txt-JUST-GENES.txt")
#
geneslist2 <- unlist(geneslist)
#
gp_up = gost(geneslist2, organism = "ggallus")
gp_up = gost(geneslist2, organism = "tguttata")
gp_up = gost(geneslist2, organism = "mmusculus")
gp_up = gost(geneslist2, organism = "hsapiens")
#
gp_up = gost(geneslist2, organism = "ggallus", significant=F)
#make a plot of terms and save it
pdf("Terms-general-downreg-GALLOS.pdf")
gostplot(gp_up, interactive =F)
dev.off()
#
pdf("Terms-general-downreg-MUSCULUS.pdf")
gostplot(gp_up, interactive =F)
dev.off()
head(gp_up$result)
#save table results for plot later
df_Place2 = data.frame(lapply(gp_up$result, as.character), stringsAsFactors=FALSE)
write.table(file="Terms-general-downreg-MUSCULUS-down.txt", df_Place2, quote=FALSE, row=TRUE)
#
df_Place2 = data.frame(lapply(gp_up$result, as.character), stringsAsFactors=FALSE)
write.table(file="Terms-general-downreg-GALLOS-down.txt", df_Place2, quote=FALSE, row=TRUE)
################################################
weight-females-downreg-GENE-list.txt-JUST-GENES.txt
geneslist <- read.table("weight-females-upreg-GENE-list.txt-JUST-GENES.txt")
geneslist2 <- unlist(geneslist)
#
gp_up = gost(geneslist2, organism = "ggallus")
gp_up = gost(geneslist2, organism = "tguttata")
gp_up = gost(geneslist2, organism = "mmusculus")
gp_up = gost(geneslist2, organism = "hsapiens")
#
gp_up = gost(geneslist2, organism = "ggallus", significant=F)
#make a plot of terms and save it
head(gp_up$result)
#save table results for plot later
gp_up = gost(geneslist2, organism = "ggallus")
df_Place2 = data.frame(lapply(gp_up$result, as.character), stringsAsFactors=FALSE)
write.table(file="Terms-weight-females-upreg-GALLOS.txt", df_Place2, quote=FALSE, row=TRUE)
pdf("Terms-weight-females-upreg-GALLOS.pdf")
gostplot(gp_up, interactive =F)
dev.off()
#
gp_up = gost(geneslist2, organism = "tguttata")
df_Place2 = data.frame(lapply(gp_up$result, as.character), stringsAsFactors=FALSE)
write.table(file="Terms-weight-females-upreg-TGUTTATA.txt", df_Place2, quote=FALSE, row=TRUE)
pdf("Terms-weight-females-upreg-TGUTTATA.pdf")
gostplot(gp_up, interactive =F)
dev.off()
#
gp_up = gost(geneslist2, organism = "mmusculus")
df_Place2 = data.frame(lapply(gp_up$result, as.character), stringsAsFactors=FALSE)
write.table(file="Terms-weight-females-upreg-MUSCULUS.txt", df_Place2, quote=FALSE, row=TRUE)
pdf("Terms-weight-females-upreg-MUSCULUS.pdf")
gostplot(gp_up, interactive =F)
dev.off()
