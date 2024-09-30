## This script contains code for morphometric and chemical data analyses ##
library(ggplot2)
MOR <- read.table("WTSH_MORPHO.txt", header=TRUE)
chem <- read.table("WTSH_CHEM.txt", header=TRUE)
juntosOR <- merge(MOR,chem, by.x=6, by.y=1, sort=F)
write.table(juntosOR, file="ALL-dataset-WTSH-01-08-2023.txt",quote=F)
str(juntosOR)
juntosOR$Plastic.y <- factor(juntosOR$Plastic.y, levels = sort(unique(juntosOR$Plastic.y)), ordered=TRUE)
juntosOR$Plastic.y
library(dplyr)
juntosOR2<-juntosOR[complete.cases(juntosOR),]
juntosOR$Plastic.y <- factor(juntosOR$Plastic.y, levels = sort(unique(juntosOR$Plastic.y)), ordered=TRUE)
juntosOR$Plastic.y
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
pp <-p + geom_bar(aes(y=BUNmg.dL,x=Bird_ID, fill=Plastic.y),stat="identity") + facet_wrap(~Plastic.y, ncol=1)
ggsave(plot=pp,"WTSH-plot-plstic-CHEM.eps",device=cairo_ps, width=13, height=5, limitsize = FALSE)
ggsave(plot=pp,"WTSH-plot-plstic-CHEM.pdf",device=cairo_ps, width=13, height=5, limitsize = FALSE)
#Linear Model
fit2 <- glm(Weight.24 ~ BUNmg.dL+ Creamg.dL+ Hct.perctPCU+Hb.g.dL +AnGapmmol.L+Nammol.L+ K.mmol.L+ Clmmol.L+ iCammol.L+ TCO2+ Glumg.dL , data=juntosOR, family="poisson")
summary(fit2)
##BUN final
ggplot(juntosOR, aes(x=BUNmg.dL, y=Weight.24)) + geom_point()+geom_smooth(method=lm) + 
  labs(y= "Body weight (g)",x="BUN (mg)") + theme_bw() + theme(axis.line = element_line(colour = "black"))
BUN3 <- ggplot(juntosOR, aes(x=BUNmg.dL, y=Weight.24)) + geom_point()+geom_smooth(method=lm) + 
  labs(y= "Body weight (g)",x="BUN (mg)") + theme_bw() + theme(axis.line = element_line(colour = "black"))
BUN_final <- BUN_final + theme(panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank())


##HCT
ggplot(juntosOR, aes(x=Hct.perctPCU, y=Weight.24)) + geom_point()+geom_smooth(method=lm) + 
  labs(y= "Body weight (g)",x="Hematocrit (%)") + theme_bw() + theme(axis.line = element_line(colour = "black"))
HCT <- ggplot(juntosOR, aes(x=Hct.perctPCU, y=Weight.24)) + geom_point()+geom_smooth(method=lm) + 
  labs(y= "Body weight (g)",x="Hematocrit (%)") + theme_bw() + theme(axis.line = element_line(colour = "black"))
HCT_final <- HCT_final + theme(panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank())

##Potassium
ggplot(juntosOR, aes(x=K.mmol.L, y=Weight.24)) + geom_point()+geom_smooth(method=lm) + 
  labs(y= "Body weight (g)",x="Potassium (mmol/L)") + theme_bw() + theme(axis.line = element_line(colour = "black"))
K <- ggplot(juntosOR, aes(x=Hct.perctPCU, y=Weight.24)) + geom_point()+geom_smooth(method=lm) + 
  labs(y= "Body weight (g)",x="Potassium (mmol/L)") + theme_bw() + theme(axis.line = element_line(colour = "black"))
K_final <- Kwin_final + theme(panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank())

##together chem 
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,3)))
print(BUN_final, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(HCT_final, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(K_final, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
##lm
weight <- ggplot(juntosOR, aes(Weight.24, as.numeric(Plastic.y) -1)) + geom_point() + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se=F) + 
  labs(y="Presence of Plastic",x="Body weight (g)") 
weight_final <- weight + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))

weight_final
#tco2
tco2 <- ggplot(juntosOR, aes(TCO2, as.numeric(Plastic.y) -1)) + geom_point() + 
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se=F) + 
  labs(y="Presence of Plastic",x="TCO2 (mmol)")
tco2_final <- tco2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))
tco2_final
#together
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))
print(weight_final, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(tco2_final, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

##fit)
##PCAmorpho##
log.ir <- log(juntosOR[,c(9:13,16)])
ir.ID <- juntosOR[, 18]
ir.pca <- prcomp(log.ir,
                 center = TRUE,
                 scale. = TRUE)
print(ir.pca)
table2 <- (ir.pca$x)
write.table(table2, file="PCA-Morphometric-SCORES-WTSH.txt")
print(log.ir)
plot(ir.pca, type= "l")
#summary method
summary(ir.pca)
#predict PCs
predict(ir.pca,
        newdata=tail(juntosOR, 28))
write.table(table, file="PCs-Morphometric-DATA-sinNA.txt")
####plot con ggbiplot
library(devtools)
##install_github("ggbiplot", "vqv")
install.packages("ggbiplot")
library(ggbiplot)
g <- ggbiplot(ir.pca, choices=1:2,obs.scale = 1, var.scale = 1, varname.size = 2, varname.adjust = 1.1,
              groups = ir.ID,
              ellipse = TRUE, 
              circle = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

g1 <- ggbiplot(ir.pca, choices=3:4, obs.scale = 1, var.scale = 1, varname.size = 2, varname.adjust = 1.1,  
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
gsm <- ggbiplot(ir.pcasm, obs.scale = 1, var.scale = 1, varname.size = 2, varname.adjust = 1.1,
                groups = ir.IDsm,
                ellipse = TRUE, 
                circle = FALSE)
gsm <- gsm + scale_color_discrete(name = '')
gsm <- gsm + theme(legend.direction = 'horizontal', 
                   legend.position = 'top')
print(gsm)
gsm2 <- ggbiplot(ir.pcasm, choices=3:4, obs.scale = 1, var.scale = 1, varname.size = 2, varname.adjust = 1.1,
                 groups = ir.IDsm,
                 ellipse = TRUE, 
                 circle = FALSE)
gsm2 <- gsm2 + scale_color_discrete(name = '')
gsm2 <- gsm2 + theme(legend.direction = 'horizontal', 
                     legend.position = 'top')
print(gsm2)
#PCA together (g, g1, gsm, gsm2)
grid.newpage()
morphotog <- pushViewport(viewport(layout = grid.layout(2,2)))
print(g, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(g1, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(gsm, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(gsm2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

## chemicalsPCA##
juntosOR2<-juntosOR[complete.cases(juntosOR),]
#apply PCA - scale. = TRUE is highly advisable, but default is FALSE
#log transform
juntosOR2
log.ir2 <- log1p(juntosOR2[, 19:29])
#juntosOR2

ir.ID2 <- juntosOR2[, 18]
ir.pca2 <- prcomp(log.ir2,
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
        newdata=tail(juntosOR2, 28))


table3 <- predict(ir.pca2,
                 newdata=tail(juntosOR2, 28))


#write.table(table, file="PCs-chemicals-DATA-sinNA.txt")
####plot con ggbiplot
library(devtools)
#install_github("ggbiplot", "vqv")
library(ggbiplot)
g2 <- ggbiplot(ir.pca2, choices=1:2, obs.scale = 1, var.scale = 1, varname.size = 2, varname.adjust = 1.1, 
               groups = ir.ID2,
               ellipse = TRUE, 
               circle = FALSE)
g2 <- g2 + scale_color_discrete(name = '')
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(g2)
g3 <- ggbiplot(ir.pca2, choices=3:4, obs.scale = 1, var.scale = 1, varname.size = 2, varname.adjust = 1.1,
               groups = ir.ID2, 
               ellipse = TRUE, 
               circle = FALSE)
g3 <- g3 + scale_color_discrete(name = '')
g3 <- g3 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(g3)
# pca for sex and chemistry panel)
log.irs <- log1p(juntosOR2[, 19:29])
ir.IDs <- juntosOR2[, 30]
ir.pcas <- prcomp(log.irs,
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
gs <- ggbiplot(ir.pcas, obs.scale = 1, var.scale = 1, varname.size = 2, varname.adjust = 1.1,
               groups = ir.IDs,
               ellipse = TRUE, 
               circle = FALSE)
gs <- gs + scale_color_discrete(name = '')
gs <- gs + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(gs)
gs2 <- ggbiplot(ir.pcas, choices=3:4, obs.scale = 1, var.scale = 1, varname.size = 2, varname.adjust = 1.1,
                groups = ir.IDs,
                ellipse = TRUE, 
                circle = FALSE)
gs2 <- gs2 + scale_color_discrete(name = '')
gs2 <- gs2 + theme(legend.direction = 'horizontal', 
                   legend.position = 'top')
print(gs2)

##chem and plastic
log.ir2 <- log(juntosOR2[, 19:28])
ir.ID2 <- juntosOR2[, 30]
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

table <- predict(ir.pca,
                 newdata=tail(juntosOR, 28))
g2 <- ggbiplot(ir.pca2, obs.scale = 1, var.scale = 1, varname.size = 2, varname.adjust = 1.1,
              groups = ir.ID2,
              ellipse = TRUE, 
              circle = FALSE)
g2 <- g2 + scale_color_discrete(name = '')
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(g2)
##with pca3 and 4
g3 <- ggbiplot(ir.pca2, choices=3:4, obs.scale = 1, var.scale = 1, varname.size = 2, varname.adjust = 1.1,
               groups = ir.ID2, 
               ellipse = TRUE, 
               circle = FALSE)
g3 <- g3 + scale_color_discrete(name = '')
g3 <- g3 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
print(g3)
##PCA together for chemical analyses gs, gs2, g2, g3
grid.newpage()
morphotog <- pushViewport(viewport(layout = grid.layout(2,2)))
print(g2, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(g3, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(gs, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(gs2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

############## bar plots######
#t-test weight male and females 
females <-which(MOR$Sex == "F")
males <- which(MOR$Sex == "M")
weight_f<-(MOR$Weight.24[females])
weight_m<-(MOR$Weight.24[males])
t.test(weight_f,weight_m)
#mean value for different columns, ignoring "NA" cells
mean(chem$iCammol.L,na.rm = T)
#how many males and females in data set?
length(which(chem$Sex == "F"))
length(which(chem$Sex == "M"))
#how many have plastic in them?
length(which(chem$Plastic == "pos"))
#how many males and females have plastic?
length(which(chem$Plastic == "pos" & chem$Sex == "M"))
length(which(chem$Plastic == "pos" & chem$Sex == "F"))
#mean value of chemicals in birds with and without plastic
#first find the rows of birds with plastic:
plasticbirds<-which(chem$Plastic == "pos")
plasticbirds
mean(chem$Nammol.L[plasticbirds],na.rm = T)
mean(chem$Nammol.L[-plasticbirds],na.rm = T)
sd(chem$Nammol.L[plasticbirds],na.rm = T)
sd(chem$Nammol.L[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
Na_pos<-(chem$Nammol.L[plasticbirds])
Na_neg<-(chem$Nammol.L[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(Na_pos,Na_neg)
## barplot(y1, col = rainbow(c(2,3)))
## help(barplot)
meanplast<-mean(chem$Nammol.L[plasticbirds],na.rm = T)
meannoplast<-mean(chem$Nammol.L[-plasticbirds],na.rm = T)
y<-c(meanplast,meannoplast)
sdplast<-sd(chem$Nammol.L[plasticbirds],na.rm = T)
sdnoplast<-sd(chem$Nammol.L[-plasticbirds],na.rm = T)
y.sd<-c(sdplast,sdnoplast)
mid<-barplot(y,col = c("dark grey", "white"), plot = F)
mid
barplot(y,col = c("dark grey", "white"),ylim=range(0,180),yaxp=c(0, 180, 18),names.arg = c("+","-"),ylab = "Sodium (mmol)",las=1) + arrows(x0=mid, y0=y-y.sd, x1=mid, y1=y+y.sd, code=3, angle=90, length=0.1)

## Bar Plot for Kmmol#
#now find the mean and s.d. of Kmmol/L for birds in these rows and in rows outside this set:
mean(chem$K.mmol.L[plasticbirds],na.rm = T)
mean(chem$K.mmol.L[-plasticbirds],na.rm = T)
sd(chem$K.mmol.L[plasticbirds],na.rm = T)
sd(chem$K.mmol.L[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
K_pos<-(chem$K.mmol.L[plasticbirds])
K_neg<-(chem$K.mmol.L[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(K_pos,K_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastk<-mean(chem$K.mmol.L[plasticbirds],na.rm = T)
meannoplastk<-mean(chem$K.mmol.L[-plasticbirds],na.rm = T)
k<-c(meanplastk,meannoplastk)
sdplastk<-sd(chem$K.mmol.L[plasticbirds],na.rm = T)
sdnoplastk<-sd(chem$K.mmol.L[-plasticbirds],na.rm = T)
y.sdk<-c(sdplastk,sdnoplastk)

midk<-barplot(k,col = c("dark grey", "white"), plot = F)
midk
kmmol <- barplot(k,col = c("dark grey", "white"),ylim=range(0,7),yaxp=c(0, 7, 2),names.arg = c("+","-"),ylab = "Potassium (mmol)",las=1) + arrows(x0=midk, y0=k-y.sdk, x1=midk, y1=k+y.sdk, code=3, angle=90, length=0.1)
print(kmmol)

## Bar plot for Clmmol 
#now find the mean and s.d. of Clmmol/L for birds in these rows and in rows outside this set:
mean(chem$Clmmol.L[plasticbirds],na.rm = T)
mean(chem$Clmmol.L[-plasticbirds],na.rm = T)
sd(chem$Clmmol.L[plasticbirds],na.rm = T)
sd(chem$Clmmol.L[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
Cl_pos<-(chem$Clmmol.L[plasticbirds])
Cl_neg<-(chem$Clmmol.L[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(Cl_pos,Cl_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastcl<-mean(chem$Clmmol.L[plasticbirds],na.rm = T)
meannoplastcl<-mean(chem$Clmmol.L[-plasticbirds],na.rm = T)
cl<-c(meanplastcl,meannoplastcl)
sdplastcl<-sd(chem$Clmmol.L[plasticbirds],na.rm = T)
sdnoplastcl<-sd(chem$Clmmol.L[-plasticbirds],na.rm = T)
y.sdcl<-c(sdplastcl,sdnoplastcl)
y.sdcl
midcl<-barplot(cl,col = c("dark grey", "white"), plot = F)
midcl
barplot(cl,col = c("dark grey", "white"),ylim=range(0,140),yaxp=c(0, 180, 18),names.arg = c("+","-"),ylab = "Chloride (mmol)",las=1) + arrows(x0=midcl, y0=cl-y.sdcl, x1=midcl, y1=cl+y.sdcl, code=3, angle=90, length=0.1)
## Bar Plot for iCammol
#now find the mean and s.d. of iCammol/L for birds in these rows and in rows outside this set:
mean(chem$iCammol.L[plasticbirds],na.rm = T)
mean(chem$iCammol.L[-plasticbirds],na.rm = T)
sd(chem$iCammol.L[plasticbirds],na.rm = T)
sd(chem$iCammol.L[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
iCa_pos<-(chem$iCammol.L[plasticbirds])
iCa_neg<-(chem$iCammol.L[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(iCa_pos,iCa_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastIC<-mean(chem$iCammol.L[plasticbirds],na.rm = T)
meannoplastIC<-mean(chem$iCammol.L[-plasticbirds],na.rm = T)
ic<-c(meanplastIC,meannoplastIC)
sdplastIC<-sd(chem$iCammol.L[plasticbirds],na.rm = T)
sdnoplastIC<-sd(chem$iCammol.L[-plasticbirds],na.rm = T)
y.sdic<-c(sdplastIC,sdnoplastIC)

midIC<-barplot(ic,col = c("dark grey", "white"), plot = F)
midIC
barplot(ic,col = c("dark grey", "white"),ylim=range(0,2),yaxp=c(0, 2, 1),names.arg = c("+","-"),ylab = "Ionized Calcium (mmol)",las=1) + arrows(x0=midIC, y0=ic-y.sdic, x1=midIC, y1=ic+y.sdic, code=3, angle=90, length=0.1)

## Bar Plot for TCO2
#now find the mean and s.d. of TCO2 for birds in these rows and in rows outside this set:
mean(chem$TCO2[plasticbirds],na.rm = T)
mean(chem$TCO2[-plasticbirds],na.rm = T)
sd(chem$TCO2[plasticbirds],na.rm = T)
sd(chem$TCO2[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
CO2_pos<-(chem$TCO2[plasticbirds])
CO2_neg<-(chem$TCO2[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(CO2_pos,CO2_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastCO2<-mean(chem$TCO2[plasticbirds],na.rm = T)
meannoplastCO2<-mean(chem$TCO2[-plasticbirds],na.rm = T)
CO2<-c(meanplastCO2,meannoplastCO2)
sdplastCO2<-sd(chem$TCO2[plasticbirds],na.rm = T)
sdnoplastCO2<-sd(chem$TCO2[-plasticbirds],na.rm = T)
y.sdCO2<-c(sdplastCO2,sdnoplastCO2)

midCO2<-barplot(CO2,col = c("dark grey", "white"), plot = F)
midCO2

barplot(CO2,col = c("dark grey", "white"),ylim=range(0,17),yaxp=c(0, 180, 18),names.arg = c("+","-"),ylab = "TCO2 (mmol)",las=1) + arrows(x0=midCO2, y0=CO2-y.sdCO2, x1=midCO2, y1=CO2+y.sdCO2, code=3, angle=90, length=0.1)

## Bar Plot for Glumg-dL
#now find the mean and s.d. of Glu for birds in these rows and in rows outside this set:
mean(chem$Glumg.dL[plasticbirds],na.rm = T)
mean(chem$Glumg.dL[-plasticbirds],na.rm = T)
sd(chem$Glumg.dL[plasticbirds],na.rm = T)
sd(chem$Glumg.dL[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
Gl_pos<-(chem$Glumg.dL[plasticbirds])
Gl_neg<-(chem$Glumg.dL[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(Gl_pos,Gl_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastGl<-mean(chem$Glumg.dL[plasticbirds],na.rm = T)
meannoplastGl<-mean(chem$Glumg.dL[-plasticbirds],na.rm = T)
Gl<-c(meanplastGl,meannoplastGl)
sdplastGl<-sd(chem$Glumg.dL[plasticbirds],na.rm = T)
sdnoplastGl<-sd(chem$Glumg.dL[-plasticbirds],na.rm = T)
y.sdGl<-c(sdplastGl,sdnoplastGl)

midGl<-barplot(Gl,c("dark grey", "white"), plot = F)
midGl
barplot(Gl,col = c("dark grey", "white"),ylim=range(0,300), yaxp=c(0, 300, 10), names.arg = c("+","-"),ylab = "Glucose (mg/dL)",las=1) + arrows(x0=midGl, y0=Gl-y.sdGl, x1=midGl, y1=Gl+y.sdGl, code=3, angle=90, length=0.1)
## Bar Plot for BUN
#now find the mean and s.d. of BUN for birds in these rows and in rows outside this set:
mean(chem$BUNmg.dL[plasticbirds],na.rm = T)
mean(chem$BUNmg.dL[-plasticbirds],na.rm = T)
sd(chem$BUNmg.dL[plasticbirds],na.rm = T)
sd(chem$BUNmg.dL[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
BUN_pos<-(chem$BUNmg.dL[plasticbirds])
BUN_neg<-(chem$BUNmg.dL[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(BUN_pos,BUN_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastBUN<-mean(chem$BUNmg.dL[plasticbirds],na.rm = T) 
meannoplastBUN<-mean(chem$BUNmg.dL[-plasticbirds],na.rm = T) 
BUN<-c(meanplastBUN,meannoplastBUN)
sdplastBUN<-sd(chem$BUNmg.dL[plasticbirds],na.rm = T)
sdnoplastBUN<-sd(chem$BUNmg.dL[-plasticbirds],na.rm = T)
y.sdBUN<-c(sdplastBUN,sdnoplastBUN)

midBUN<-barplot(BUN,col = c("dark grey", "white"), plot = F)
midBUN
barplot(BUN,col = c("dark grey", "white"),ylim=range(0,20),yaxp=c(0, 180, 18),names.arg = c("+","-"),ylab = "BUN (mg/dl)",las=1) + arrows(x0=midBUN, y0=BUN-y.sdBUN, x1=midBUN, y1=BUN+y.sdBUN, code=3, angle=90, length=0.1)

## Bar Plot for Hct
#now find the mean and s.d. of Hct for birds in these rows and in rows outside this set:
mean(chem$Hct.perctPCU[plasticbirds],na.rm = T)
mean(chem$Hct.perctPCU[-plasticbirds],na.rm = T)
sd(chem$Hct.perctPCU[plasticbirds],na.rm = T)
sd(chem$Hct.perctPCU[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
Hct_pos<-(chem$Hct.perctPC[plasticbirds])
Hct_neg<-(chem$Hct.perctPC[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(Hct_pos,Hct_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastHct<-mean(chem$Hct.perctPC[plasticbirds],na.rm = T)
meannoplastHct<-mean(chem$Hct.perctPC[-plasticbirds],na.rm = T)
Hct<-c(meanplastHct,meannoplastHct)
sdplastHct<-sd(chem$Hct.perctPCL[plasticbirds],na.rm = T)
sdnoplastHct<-sd(chem$Hct.perctPC[-plasticbirds],na.rm = T)
y.sdHct<-c(sdplastHct,sdnoplastHct)

midHct<-barplot(Hct,col = c("dark grey", "white"), plot = F)
midHct
barplot(Hct,col = c("dark grey", "white"),ylim=range(0,45),yaxp=c(0, 180, 18),names.arg = c("+","-"),ylab = "Hematocrit (%)",las=1) + arrows(x0=midHct, y0=Hct-y.sdHct, x1=midHct, y1=Hct+y.sdHct, code=3, angle=90, length=0.1)

## Bar Plot for Hb-g-dL
#now find the mean and s.d. of Hb.g.dL for birds in these rows and in rows outside this set:
mean(chem$Hb.g.dL[plasticbirds],na.rm = T) 
mean(chem$Hb.g.dL[-plasticbirds],na.rm = T)
sd(chem$Hb.g.dL[plasticbirds],na.rm = T)
sd(chem$Hb.g.dL[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
Hg_pos<-(chem$Hb.g.dL[plasticbirds])
Hg_neg<-(chem$Hb.g.dL[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(Hg_pos,Hg_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastHg<-mean(chem$Hb.g.dL[plasticbirds],na.rm = T)
meannoplastHg<-mean(chem$Hb.g.dL[-plasticbirds],na.rm = T)
Hg<-c(meanplastHg,meannoplastHg)
sdplastHg<-sd(chem$Hb.g.dL[plasticbirds],na.rm = T)
sdnoplastHg<-sd(chem$Hb.g.dL[-plasticbirds],na.rm = T)
y.sdHg<-c(sdplastHg,sdnoplastHg)

midHg<-barplot(Hg,col = c("dark grey", "white"), plot = F)
midHg
barplot(Hg,col = c("dark grey", "white"),ylim=range(0,15),yaxp=c(0, 180, 18),names.arg = c("+","-"),ylab = "Hemoglobin (g/dL)",las=1) + arrows(x0=midHg, y0=Hg-y.sdHg, x1=midHg, y1=Hg+y.sdHg, code=3, angle=90, length=0.1)

## Bar plot for AnGap mmol/L 
#now find the mean and s.d. of Clmmol/L for birds in these rows and in rows outside this set:
mean(chem$AnGapmmol.L[plasticbirds],na.rm = T)
mean(chem$AnGapmmol.L[-plasticbirds],na.rm = T)
sd(chem$AnGapmmol.L[plasticbirds],na.rm = T)
sd(chem$AnGapmmol.L[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
An_pos<-(chem$AnGapmmol.L[plasticbirds])
An_neg<-(chem$AnGapmmol.L[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(An_pos,An_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastAn<-mean(chem$AnGapmmol.L[plasticbirds],na.rm = T)
meannoplastAn<-mean(chem$AnGapmmol.L[-plasticbirds],na.rm = T)
An<-c(meanplastAn,meannoplastAn)
sdplastAn<-sd(chem$AnGapmmol.L[plasticbirds],na.rm = T)
sdnoplastAn<-sd(chem$AnGapmmol.L[-plasticbirds],na.rm = T)
y.sdAn<-c(sdplastAn,sdnoplastAn)

midAn<-barplot(An,col = c("dark grey", "white"), plot = F)
midAn
barplot(An,col = c("dark grey", "white"),ylim=range(0,25),yaxp=c(0, 180, 18),names.arg = c("+","-"),ylab = "Anion Gap (mmol)",las=1) + arrows(x0=midAn, y0=An-y.sdAn, x1=midAn, y1=An+y.sdAn, code=3, angle=90, length=0.1)

par()
## Bar plot for Crea 
#now find the mean and s.d. of Crea for birds in these rows and in rows outside this set:
mean(chem$Creamg.dL[plasticbirds],na.rm = T)
mean(chem$Creamg.dL[-plasticbirds],na.rm = T)
sd(chem$Creamg.dL[plasticbirds],na.rm = T)
sd(chem$Creamg.dL[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
Cr_pos<-(chem$Creamg.dL[plasticbirds])
Cr_neg<-(chem$Creamg.dL[-plasticbirds])
##y1 <- c(Cr_pos, Cr_neg)
t.test(Cr_pos,Cr_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastCr<-mean(chem$Creamg.dL[plasticbirds],na.rm = T)
meannoplastCr<-mean(chem$Creamg.dL[-plasticbirds],na.rm = T)
Cr<-c(meanplastCr,meannoplastCr)
sdplastCr<-sd(chem$Creamg.dL[plasticbirds],na.rm = T)
sdnoplastCr<-sd(chem$Creamg.dL[-plasticbirds],na.rm = T)
y.sdCr<-c(sdplastCr,sdnoplastCr)

midCr<-barplot(Cr,col = c("dark grey", "white"), plot = F)
midCr
barplot(Cr,col = c("dark grey", "white"),ylim=range(0,1),yaxp=c(0, 1, 1),names.arg = c("+","-"),ylab = "Creatine (mg/dL)",las=1) + arrows(x0=midCr, y0=Cr-y.sdCr, x1=midCr, y1=Cr+y.sdCr, code=3, angle=90, length=0.1)
help(attach)


## Bar Plot for weight and plastic
#now find the mean and s.d. of weight for birds in these rows and in rows outside this set:

## mean(MOR$Weight.24[plasticbirds],na.rm = T)
plasticbirds2<-which(MOR$Plastic == "+")
mean(MOR$Weight.24[plasticbirds2])
mean(MOR$Weight.24[-plasticbirds2])
#Find the standard deviation
sd(MOR$Weight.24[plasticbirds2])
sd(MOR$Weight.24[-plasticbirds2])
#t-test of differences in mean between these groups
W_pos<-(MOR$Weight.24[plasticbirds2])
W_neg<-(MOR$Weight.24[-plasticbirds2])
##y1 <- c(Na_pos, Na_neg)
t.test(W_pos,W_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastW<-mean(MOR$Weight.24[plasticbirds2],na.rm = T)
meannoplastW<-mean(MOR$Weight.24[-plasticbirds2],na.rm = T)
W<-c(meanplastW,meannoplastW)
sdplastW<-sd(MOR$Weight.24[plasticbirds2],na.rm = T)
sdnoplastW<-sd(MOR$Weight.24[-plasticbirds2],na.rm = T)
y.sdW<-c(sdplastW,sdnoplastW)

midW<-barplot(W,col = c("dark grey", "white"), plot = F)
midW
barplot(W,col = c("dark grey", "white"),ylim=range(0,600),yaxp=c(0, 500, 10),names.arg = c("+","-"),ylab = "Body Weight (grams)",las=1) + arrows(x0=midW, y0=W-y.sdW, x1=midW, y1=W+y.sdW, code=3, angle=90, length=0.1)
### together barplots ###

par(mfrow = c(4, 3))


##Shapiro-Wilk normality test
shapiro.test(mannwhitweight$weight)

## Mann-Whitney U Test
#create data frame 
mannwhitweight <- data.frame(weight = c(juntosOR$Weight.24), category = c(juntosOR$wfact), plastic = c(juntosOR$Plastic.y))
head(mannwhitweight)
#run Mann-Whitney U test
test_result <- wilcox.test(weight ~ plastic, data = mannwhitweight)
print(test_result)
