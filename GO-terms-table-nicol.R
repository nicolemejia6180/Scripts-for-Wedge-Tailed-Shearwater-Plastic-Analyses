library(ggplot2)
library(scales)
library(tidyverse)
setwd <-("/home/flavia/Desktop/")
CN <- read.table("Terms-big-table.txt", header=T,stringsAsFactors = F)
#Weight-all,Weight-females
x <- subset(CN, Experimental == "Weight-all")
x <- x[!duplicated(x$term_name), ]
UP <- subset(x, Up.down == "up_reg")
DO <- subset(x, Up.down == "down_reg")


#x$Up.down <- as.factor(x$Up.down)
#
#is.factor(x$Up.down)
#x$Up.down2[x$Up.down == "up_reg"] <- "1"
#x$Up.down2[x$Up.down == "down_reg"] <- "2"
#
#x$Up.down2 <- as.factor(x$Up.down2)
#is.factor(x$Up.down2)
#x$Up.down <- factor(x$Up.down,levels = c("up_reg", "down_reg"))
#x$Up.down2 <- factor(x$Up.down2,levels = c("1", "2"))
#
#x2 <- x[order(x$Up.down2),]
#
ggplot(UP, aes(x=p_value, y=term_name))+
  geom_bar(stat='identity')+
   theme_bw() +
  theme(axis.text.x=element_text(size = 14, angle = 090, hjust = 1, vjust = 0, face = "plain"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
        axis.title.y = element_text(size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
        axis.title=element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
#
p <- ggplot(DO, aes(x=p_value, y=reorder(term_name, -p_value)))+
  geom_bar(stat='identity')+
  theme_bw() +
  theme(axis.text.x=element_text(size = 12, angle = 090, hjust = 1, vjust = 0, face = "plain"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
        axis.title.y = element_text(size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
        axis.title=element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
#
p + coord_cartesian(xlim=c(0.0000000000001,0.000001))
