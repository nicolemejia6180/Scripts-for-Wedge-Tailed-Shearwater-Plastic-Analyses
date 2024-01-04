grep -H '% of reads mapped to multiple loci' *final.out > WTSH-multiple-map.txt

grep -H 'Uniquely mapped reads %' *final.out > WTSH-uniq-map.txt

#check headers for R and format as needed
#
library(ggplot2)
library(scales)
library(tidyverse)

CN <- read.table("Format-Cyncas-uniq-map.txt.txt", header=F,stringsAsFactors = F)
#
BM$V1 <-as.numeric(BM$V1)
sample<-c(1:(nrow(BM)))
BM2 <- Bm[order(BM$V1),]
ordered = BM[order(BM$species,BM$V1),]

p<- ggplot(ordered) +
  geom_bar( aes(y=V1, x=sample, fill=species), stat="identity") +
  labs(x="", y = "Uniquely mapped reads %")+ theme_bw() +
  theme(axis.text.x=element_blank(), axis.line = element_line(colour = "black"),
	axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
	axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
	axis.title=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+ scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_fill_manual("Species", values = c("Brown jay" = "chocolate4", "Green jay" = "green4", "Yucatan jay" = "deepskyblue3"))
ggsave(plot=p,"Uniq-mapped-reads-jays.eps",device=cairo_ps)
#
CN <- read.table("Format-Cyncas-multiple-map.txt.txt", header=F,stringsAsFactors = F)
#
BM$V1 <-as.numeric(BM$V1)
sample<-c(1:(nrow(BM)))
BM2 <- Bm[order(BM$V1),]
ordered = BM[order(BM$species,BM$V1),]

p<- ggplot(ordered) +
  geom_bar( aes(y=V1, x=sample, fill=species), stat="identity") +
  labs(x="", y = "Multiple mapped reads %")+ theme_bw() +
  theme(axis.text.x=element_blank(), axis.line = element_line(colour = "black"),
	axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
	axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
	axis.title=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+ scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + scale_fill_manual("Species", values = c("Brown jay" = "chocolate4", "Green jay" = "green4", "Yucatan jay" = "deepskyblue3"))
ggsave(plot=p,"multiple-mapped-reads-jays.eps",device=cairo_ps)

