setwd("~/Documents/Thesis_Docs")
library(ggplot2)
#read files
MOR <- read.table("WTSH_MORPHO.txt", header=TRUE)
chem <- read.table("WTSH_CHEM.txt", header=TRUE)
#merge files into one table
juntosOR <- merge(MOR,CHEM, by.x=6, by.y=1, sort=F)
write.table(juntosOR, file="ALL-dataset-WTSH-01-07-2022.txt",quote=F)
#prepare variables
#check info of object to declare it as you wish
str(juntosOR)
juntosOR$Bird_ID <- as.factor(juntosOR$Bird_ID)
juntosOR$Plastic.y <- as.factor(juntosOR$Plastic.y)
#keep complate cases
juntosOR2<-juntosOR[complete.cases(juntosOR),]
#order
juntosOR$Plastic.y <- factor(juntosOR$Plastic.y, levels = sort(unique(juntosOR$Plastic.y)), ordered=TRUE)
juntosOR$Plastic.y
#order with dplyr
library(dplyr)
juntosOR2 <- juntosOR %>% group_by(Plastic.y) 
#make plot with all samples per category
p<-ggplot(juntosOR2)
pp <-p + geom_bar(aes(y=BUNmg.dL,x=Plastic.y,fill=Plastic.y),stat="identity") +
theme_bw()+ labs(y = "BUNmg.dLe",x = "Bird_ID")+
  theme(axis.text.x=element_text(size = 12, angle = 90, hjust = 1, vjust = 0.3, face = "plain"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
        axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
#
pp <-p + geom_bar(aes(y=BUNmg.dL,x=Bird_ID, fill=Plastic.y),stat="identity") +
  theme_bw()+ labs(y = "BUNmg.dLe",x = "Bird_ID")+
  theme(axis.text.x=element_text(size = 12, angle = 90, hjust = 1, vjust = 0.3, face = "plain"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
        axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

#facet wrap type
pp <-p + geom_bar(aes(y=BUNmg.dL,x=Bird_ID, fill=Plastic.y),stat="identity") + facet_wrap(~Plastic.y, ncol=1) 
#
ggsave(plot=pp,"WTSH-plot-plstic-CHEM.eps",device=cairo_ps, width=13, height=5, limitsize = FALSE)
#############MORPHOMETRICS##########################
#make PCA##########################################
#apply PCA - scale. = TRUE is highly
#advisable, but default is FALSE
#log transform
##log.ir <- log(juntosOR[, 8:13]) 
##ir.ID <- juntosOR[, 18]
##ir.pca <- prcomp(log.ir,
                 #center = TRUE,
                 #scale. = TRUE) 
log.ir <- log(juntosOR[,c(9:13,16)])
ir.ID <- juntosOR[, 18]
ir.pca <- prcomp(log.ir,
                 center = TRUE,
                 scale. = TRUE)
#print method
print(ir.pca)
table2 <- (ir.pca$x)
write.table(table2, file="PCA-Morfometric-SCORES-WTSH.txt")
print(log.ir)
# esta figura nos permitira decidir cuantos PCs retendremos para subsecuentes analisis
#plot metho
plot(ir.pca, type= "l")
#summary method
summary(ir.pca)
#predict PCs
predict(ir.pca,
        newdata=tail(juntosOR, 28))

table <- predict(ir.pca,
                 newdata=tail(juntosOR, 28))
write.table(table, file="PCs-Morfometric-DATA-sinNA.txt")
####plot con ggbiplot
install.packages("devtools")
library(devtools)
##install_github("vqv/ggbiplot")
##install_github("ggbiplot", "vqv")
library(ggbiplot)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1,
              groups = ir.ID,
              ellipse = TRUE, 
              circle = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
g1 <- ggbiplot(ir.pca, choices=3:4, obs.scale = 1, var.scale = 1, 
               groups = ir.ID, 
               ellipse = TRUE, 
               circle = FALSE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(g1)
g2 <- ggbiplot(ir.pca, choices=2:3, obs.scale = 1, var.scale = 1, 
               groups = ir.ID, 
               ellipse = TRUE, 
               circle = FALSE)
g2 <- g2 + scale_color_discrete(name = '')
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(g2)
##############pca WITH CHEMICALS###########
#make PCA##################################
#first remove samples with no data
juntosOR2<-juntosOR[complete.cases(juntosOR),]
#apply PCA - scale. = TRUE is highly advisable, but default is FALSE
#log transform
log.ir <- log(juntosOR2[, 19:28])
ir.ID <- juntosOR2[, 30]
ir.pca <- prcomp(log.ir,
                 center = TRUE,
                 scale. = TRUE) 
#print method
print(ir.pca)
table2 <- (ir.pca$x)
write.table(table2, file="PCA-Chemicals-SCORES-WTSH.txt")
print(log.ir)
# this figure shows how many PCs we should retain for the next analyses
#plot metho
plot(ir.pca, type= "l")
#summary method
summary(ir.pca)
#predict PCs
predict(ir.pca,
        newdata=tail(juntosOR, 28))

table <- predict(ir.pca,
                 newdata=tail(juntosOR, 28))
#write.table(table, file="PCs-chemicals-DATA-sinNA.txt")
####plot con ggbiplot
library(devtools)
#install_github("ggbiplot", "vqv")
library(ggbiplot)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1,
              groups = ir.ID,
              ellipse = TRUE, 
              circle = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
g1 <- ggbiplot(ir.pca, choices=3:4, obs.scale = 1, var.scale = 1, 
               groups = ir.ID, 
               ellipse = TRUE, 
               circle = FALSE)
g1 <- g1 + scale_color_discrete(name = '')
g1 <- g1 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(g1)
g2 <- ggbiplot(ir.pca, choices=2:3, obs.scale = 1, var.scale = 1, 
               groups = ir.ID, 
               ellipse = TRUE, 
               circle = FALSE)
g2 <- g2 + scale_color_discrete(name = '')
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(g2)
#################GLM analyses to explore variable with associations###########
#all chemicals
fit <- glm(Weight.24 ~ BUNmg.dL+ Creamg.dL+ Hct.perctPCU+Hb.g.dL +AnGapmmol.L+Nammol.L+ K.mmol.L+ Clmmol.L+ iCammol.L+ TCO2+ Glumg.dL , data=juntosOR, family="poisson")
summary(fit)
#
fit <- glm(Plastic.y ~ Weight.24+BUNmg.dL+ Creamg.dL+ Hct.perctPCU+Hb.g.dL +AnGapmmol.L+Nammol.L+ K.mmol.L+ Clmmol.L+ iCammol.L+ TCO2+ Glumg.dL , data=juntosOR, family = binomial(link = "logit"))
summary(fit)
#simplify the model remove colineal variables
fit <- glm(Plastic.y ~ Sex.y, data=juntosOR, family = binomial(link = "logit"))
summary(fit)
#
fit <- glm(Plastic.y ~ Weight.24+BUNmg.dL+ Creamg.dL, data=juntosOR, family = binomial(link = "logit"))
summary(fit)
#
fit <- glm(Sex.y ~ Weight.24, data=juntosOR, family = "poisson")
summary(fit)

