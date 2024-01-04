setwd("/home/flavia/Desktop/Nicole-Thesis/Data-R/")
library(ggplot2)
pdf("metadata-nicol.pdf")
#read files
MOR <- read.table("WTSH_MORPHO.txt", header=TRUE)
CHEM <- read.table("WTSH_CHEM.txt", header=TRUE)
ID <- read.table("Sample_Individual-IDs.txt", header=TRUE)
#merge files into one table
juntosOR <- merge(MOR,CHEM, by.x=6, by.y=1, sort=F)
juntosOR
write.table(juntosOR, file="ALL-dataset-WTSH-01-07-2022.txt",quote=F)
#include sample ID and transcriptomic sample ID
juntosOR3 <- merge(juntosOR,ID, by.x=1, by.y=3, sort=F)
write.table(juntosOR3, file="sampleID-ALL-dataset-WTSH-01-19-2022.txt",quote=F)
juntosOR3
#prepare variables
#check info of object to declare it as you wish
str(juntosOR)
juntosOR$Bird_ID <- as.factor(juntosOR$Bird_ID)
juntosOR$Plastic.y <- as.factor(juntosOR$Plastic.y)
#keep complate cases
juntosOR2<-juntosOR[complete.cases(juntosOR),]
juntosOR2

#
#############MORPHOMETRICS##########################
#make PCA##########################################
#apply PCA - scale. = TRUE is highly
#advisable, but default is FALSE
#log transform morpho variables
##log.ir <- log(juntosOR[,c(9:11,16)])
##ir.ID <- juntosOR2[, 16]
##ir.pca <- prcomp(log.ir,
#                center = TRUE,
#                  scale. = TRUE) 
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
library(devtools)
##install_github("ggbiplot", "vqv")
install.packages("ggbiplot")
library(ggbiplot)
g <- ggbiplot(ir.pca, choices=1:2,obs.scale = 1, var.scale = 1,
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
# pca for sex (morpho)
log.irsm <- log(juntosOR2[,c(9:13,16)])
log.irsm
ir.IDsm <- juntosOR2[, 14]
juntosOR2
ir.pcasm <- prcomp(log.irsm,
                  center = TRUE,
                  scale. = TRUE)


print(ir.pcasm)
table2 <- (ir.pcasm$x)
write.table(table2, file="PCA-Chemicals-SCORES-WTSH.txt")
print(log.irsm)
plot(ir.pcasm, type= "l")
summary(ir.pcasm)
predict(ir.pcasm,
        newdata=tail(juntosOR, 28))
table <- predict(ir.pcasm,
                 newdata=tail(juntosOR, 28))
gsm <- ggbiplot(ir.pcasm, obs.scale = 1, var.scale = 1,
               groups = ir.IDsm,
               ellipse = TRUE, 
               circle = FALSE)
gsm <- gsm + scale_color_discrete(name = '')
gsm <- gsm + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(gsm)
gsm2 <- ggbiplot(ir.pcasm, choices=3:4, obs.scale = 1, var.scale = 1,
                groups = ir.IDsm,
                ellipse = TRUE, 
                circle = FALSE)
gsm2 <- gsm2 + scale_color_discrete(name = '')
gsm2 <- gsm2 + theme(legend.direction = 'horizontal', 
                   legend.position = 'top')
print(gsm2)

##############pca WITH CHEMICALS###########
#make PCA##################################
#first remove samples with no data
juntosOR2<-juntosOR[complete.cases(juntosOR),]
#apply PCA - scale. = TRUE is highly advisable, but default is FALSE
#log transform
juntosOR2
log.ir2 <- log1p(juntosOR2[, 19:29])
ir.ID2 <- juntosOR2[, 18]
ir.pca2 <- prcomp(log.ir,
                 center = TRUE,
                 scale. = TRUE)
#print method
print(ir.pca2)
table2 <- (ir.pca2$x)
write.table(table2, file="PCA-Chemicals-SCORES-WTSH.txt")
print(log.ir2)
# this figure shows how many PCs we should retain for the next analyses
#plot metho
plot(ir.pca2, type= "l")
#summary method
summary(ir.pca2)
#predict PCs
predict(ir.pca2,
        newdata=tail(juntosOR, 28))

table <- predict(ir.pca2,
                 newdata=tail(juntosOR, 28))
#write.table(table, file="PCs-chemicals-DATA-sinNA.txt")
####plot con ggbiplot
library(devtools)
#install_github("ggbiplot", "vqv")
library(ggbiplot)
g2 <- ggbiplot(ir.pca2, obs.scale = 1, var.scale = 1,
              groups = ir.ID2,
              ellipse = TRUE, 
              circle = FALSE)
g2 <- g2 + scale_color_discrete(name = '')
g2 <- g2 + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g2)
g3 <- ggbiplot(ir.pca2, choices=3:4, obs.scale = 1, var.scale = 1, 
               groups = ir.ID2, 
               ellipse = TRUE, 
               circle = FALSE)
g3 <- g3 + scale_color_discrete(name = '')
g3 <- g3 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(g3)
g4 <- ggbiplot(ir.pca2, choices=2:3, obs.scale = 1, var.scale = 1, 
               groups = ir.ID2, 
               ellipse = TRUE, 
               circle = FALSE)
g4 <- g4 + scale_color_discrete(name = '')
g4 <- g4 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(g4)
# pca for sex (chem)
log.irs <- log1p(juntosOR2[, 19:29])
ir.IDs <- juntosOR2[, 30]
ir.pcas <- prcomp(log.ir,
                  center = TRUE,
                  scale. = TRUE)
print(ir.pcas)
table2 <- (ir.pcas$x)
write.table(table2, file="PCA-Chemicals-SCORES-WTSH.txt")
print(log.irs)
plot(ir.pcas, type= "l")
summary(ir.pcas)
predict(ir.pcas,
        newdata=tail(juntosOR, 28))
table <- predict(ir.pcas,
                 newdata=tail(juntosOR, 28))
gs <- ggbiplot(ir.pcas, obs.scale = 1, var.scale = 1,
               groups = ir.IDs,
               ellipse = TRUE, 
               circle = FALSE)
gs <- gs + scale_color_discrete(name = '')
gs <- gs + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(gs)
gs2 <- ggbiplot(ir.pcas, choices=3:4, obs.scale = 1, var.scale = 1,
               groups = ir.IDs,
               ellipse = TRUE, 
               circle = FALSE)
gs2 <- gs2 + scale_color_discrete(name = '')
gs2 <- gs2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(gs2)
gs3 <- ggbiplot(ir.pcas, choices=2:3, obs.scale = 1, var.scale = 1,
                groups = ir.IDs,
                ellipse = TRUE, 
                circle = FALSE)
gs3 <- gs3 + scale_color_discrete(name = '')
gs3 <- gs3 + theme(legend.direction = 'horizontal', 
                   legend.position = 'top')
print(gs3)
#################GLM analyses to explore variable with associations###########
#
fit <- glm(Sex.y ~ Weight.24 + Bill_W+ Bill_D+Wing_chrd+ Tarsus+ Bill_L,data=juntosOR, family = binomial(link = "logit"))
summary(fit)
# plastic and weight
fit1 <- glm(Plastic.y ~ Weight.24 + Bill_W+ Bill_D+Wing_chrd+ Tarsus+ Bill_L,data=juntosOR, family = binomial(link = "logit"))
summary(fit1)
#
#Weight with all chemicals
fit2 <- glm(Weight.24 ~ BUNmg.dL+ Creamg.dL+ Hct.perctPCU+Hb.g.dL +AnGapmmol.L+Nammol.L+ K.mmol.L+ Clmmol.L+ iCammol.L+ TCO2+ Glumg.dL , data=juntosOR, family="poisson")
summary(fit2)
#Plastic with all chemicals
fit3 <- glm(Plastic.y ~ Weight.24+BUNmg.dL+ Creamg.dL+ Hct.perctPCU+Hb.g.dL +AnGapmmol.L+Nammol.L+ K.mmol.L+ Clmmol.L+ iCammol.L+ TCO2+ Glumg.dL , data=juntosOR, family = binomial(link = "logit"))
summary(fit3)
#
juntosOR$Sex.y <- as.factor(juntosOR$Sex.y)
fit4 <- glm(Sex.y ~ Weight.24+BUNmg.dL+ Creamg.dL+ Hct.perctPCU+Hb.g.dL +AnGapmmol.L+Nammol.L+ K.mmol.L+ Clmmol.L+ iCammol.L+ TCO2+ Glumg.dL , data=juntosOR, family = binomial(link = "logit"))
summary(fit4)
######################
#regression plots
ggplot(juntosOR, aes(x=BUNmg.dL, y=Weight.24)) + geom_point()+geom_smooth(method=lm) + 
labs(y= "Weight",x="BUNmg") + theme_bw() + theme(axis.line = element_line(colour = "black")) #
#
ggplot(juntosOR, aes(x=Hct.perctPCU, y=Weight.24)) + geom_point()+geom_smooth(method=lm) + 
  labs(y= "Weight",x="Hct.perctPCU") + theme_bw() + theme(axis.line = element_line(colour = "black"))
#
ggplot(juntosOR, aes(x=K.mmol.L, y=Weight.24)) + geom_point()+geom_smooth(method=lm) + 
  labs(y= "Weight",x="K.mmol.L") + theme_bw() + theme(axis.line = element_line(colour = "black"), axis.text.x = element_blank(),
                                                          panel.grid.major = element_blank(),
                                                          panel.grid.minor = element_blank(),
                                                          panel.background = element_blank())

#
#ALL WITH PLASTIC POSITIVE AND NEGATIVE: JUST SHOWING THE SIGNIFICANT CURVES
fit5 <- glm(Weight.24 ~ Plastic.y , data=juntosOR, family="poisson")
summary(fit5)
#
lm5 <-lm(Weight.24 ~ Plastic.y, data=juntosOR)
summary(lm5)
#weight and plastic 
ggplot(juntosOR, aes(Weight.24, as.numeric(Plastic.y) -1)) + geom_point() + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se=F) + 
  labs(y="Plastic",x="Weight")
#plastic and tco2
ggplot(juntosOR, aes(TCO2, as.numeric(Plastic.y) -1)) + geom_point() + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se=F) + 
  labs(y="Plastic",x="TCO2")
#
lm2 <-lm(TCO2 ~ Plastic.y, data=juntosOR)
summary(lm2)
#
############Weight####

juntosOR$Sex.y <- as.factor(juntosOR$Sex.y)
fit2 <- glm(Sex.y ~ Weight.24 +BUNmg.dL+ Creamg.dL+ Hct.perctPCU+Hb.g.dL +AnGapmmol.L+Nammol.L+ K.mmol.L+ Clmmol.L+ iCammol.L+ TCO2+ Glumg.dL , data=juntosOR, family = binomial(link = "logit"))
summary(fit2)
#
lm3 <-lm(Weight.24 ~ Sex.y, data=juntosOR)
summary(lm3)
#
#make categories for weight
wfact <- cut(juntosOR2$Weight.24,3,labels=c('Low','Medium','High'))
table(wfact)
juntosOR2$wfact <- wfact
write.table(juntosOR2, file="WTSH-with-weight-as-factors-3.txt")
#
wfact <- cut(juntosOR2$Weight.24,2,labels=c('Low','High'))
table(wfact)
juntosOR2$wfact <- wfact
write.table(juntosOR2, file="WTSH-with-weight-as-factors-2.txt")
library("devtools")
dev.off()
help(dev.off)
