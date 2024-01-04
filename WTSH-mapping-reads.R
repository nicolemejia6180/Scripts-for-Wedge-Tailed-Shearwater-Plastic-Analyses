library(ggplot2)
library(scales)
library(tidyverse)
setwd <-("/home/flavia/Desktop/")
CN <- read.table("WTSH-uniq-map.txt", header=F,stringsAsFactors = F)
#
CN$V2 <-as.numeric(CN$V2)
CN$V1 <-as.factor(CN$V1)
#
sample<-c(1:(nrow(CN)))
ordered <- CN[order(CN$V2),]

p<- ggplot(ordered,aes(y=V2, x=reorder(V1, -V2))) +
  geom_bar(stat="identity") +
  labs(x="Samples", y = "Uniquely mapped reads %")+ theme_bw() +
  theme(axis.text.x=element_text(size = 5, angle = 090, hjust = 1, vjust = 0, face = "plain"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
        axis.title.y = element_text(size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
        axis.title=element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+ scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

#
#ggsave(plot=p,"Uniq-mapped-reads-jays.eps",device=cairo_ps)
#
CN <- read.table("WTSH-multiple-map.txt", header=F,stringsAsFactors = F)
#
CN$V2 <-as.numeric(CN$V2)
CN$V1 <-as.factor(CN$V1)
#
sample<-c(1:(nrow(CN)))
ordered <- CN[order(CN$V2),]
p<- ggplot(ordered,aes(y=V2, x=reorder(V1, -V2))) +
  geom_bar(stat="identity") +
  labs(x="Samples", y = "Multiple mapped reads %")+ theme_bw() +
  theme(axis.text.x=element_text(size = 5, angle = 090, hjust = 1, vjust = 0, face = "plain"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
        axis.title.y = element_text(size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
        axis.title=element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+ scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
