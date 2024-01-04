chem<-read.table("/Users/scottedwards/OneDrive\ -\ Harvard\ University/Students/Nicole/Thesis/WTSH_chem_Sheet1_v2b.txt",sep="\t",header = T)


head(chem)
setwd("~/Documents/Thesis_Docs/WTSH_CHEM.txt")


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

#now find the mean and s.d. of Nammol/L for birds in these rows and in rows outside this set:
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

#to make error bars, first create objects with mean and sd for a given variable
meanplast<-mean(chem$Nammol.L[plasticbirds],na.rm = T)
meannoplast<-mean(chem$Nammol.L[-plasticbirds],na.rm = T)
y<-c(meanplast,meannoplast)
sdplast<-sd(chem$Nammol.L[plasticbirds],na.rm = T)
sdnoplast<-sd(chem$Nammol.L[-plasticbirds],na.rm = T)
y.sd<-c(sdplast,sdnoplast)

#now create midpoints for where you want the bars and error bars to go
mid<-barplot(y,col = rainbow(c(2,3)), plot = F)
mid

#now plot the barplot first and then use the arrows function and 90 degree arrows to create standard error bars
barplot(y,col = rainbow(c(2,3)),ylim=range(0,180),yaxp=c(0, 180, 18),names.arg = c("with plastic","without plastic"),ylab = "Mmol Na",las=1)
arrows(x0=mid, y0=y-y.sd, x1=mid, y1=y+y.sd, code=3, angle=90, length=0.1)

## Bar Plot for Kmmol
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

midk<-barplot(k,col = rainbow(c(2,3)), plot = F)
midk
barplot(k,col = rainbow(c(2,3)),ylim=range(0,7),yaxp=c(0, 7, 2),names.arg = c("with plastic","without plastic"),ylab = "Mmol K",las=1)
arrows(x0=midk, y0=k-y.sdk, x1=midk, y1=k+y.sdk, code=3, angle=90, length=0.1)


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

midcl<-barplot(cl,col = rainbow(c(2,3)), plot = F)
midcl
barplot(cl,col = rainbow(c(2,3)),ylim=range(0,140),yaxp=c(0, 180, 18),names.arg = c("with plastic","without plastic"),ylab = "Mmol Cl",las=1)
arrows(x0=midcl, y0=cl-y.sdcl, x1=midcl, y1=cl+y.sdcl, code=3, angle=90, length=0.1)

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

midIC<-barplot(ic,col = rainbow(c(2,3)), plot = F)
midIC
barplot(ic,col = rainbow(c(2,3)),ylim=range(0,2),yaxp=c(0, 2, 1),names.arg = c("with plastic","without plastic"),ylab = "Mmol iCa",las=1)
arrows(x0=midIC, y0=ic-y.sdic, x1=midIC, y1=ic+y.sdic, code=3, angle=90, length=0.1)

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

midCO2<-barplot(CO2,col = rainbow(c(2,3)), plot = F)
midCO2
barplot(CO2,col = rainbow(c(2,3)),ylim=range(0,17),yaxp=c(0, 180, 18),names.arg = c("with plastic","without plastic"),ylab = "Mmol TCO2",las=1)
arrows(x0=midCO2, y0=CO2-y.sdCO2, x1=midCO2, y1=CO2+y.sdCO2, code=3, angle=90, length=0.1)

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

midGl<-barplot(Gl,col = rainbow(c(2,3)), plot = F)
midGl
help(barplot)
barplot(Gl,col = rainbow(c(2,3)),ylim=range(0,300), yaxp=c(0, 300, 10), names.arg = c("with plastic","without plastic"),ylab = "Mmol Glu mg/dL",las=1)
arrows(x0=midGl, y0=Gl-y.sdGl, x1=midGl, y1=Gl+y.sdGl, code=3, angle=90, length=0.1)

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

midBUN<-barplot(BUN,col = rainbow(c(2,3)), plot = F)
midBUN
barplot(BUN,col = rainbow(c(2,3)),ylim=range(0,20),yaxp=c(0, 180, 18),names.arg = c("with plastic","without plastic"),ylab = "mg/dL BUN",las=1)
arrows(x0=midBUN, y0=BUN-y.sdBUN, x1=midBUN, y1=BUN+y.sdBUN, code=3, angle=90, length=0.1)

## Bar Plot for Hct
#now find the mean and s.d. of Hct for birds in these rows and in rows outside this set:
mean(chem$Hct.perctPCU[plasticbirds],na.rm = T)
mean(chem$Hct.perctPC[-plasticbirds],na.rm = T)
sd(chem$Hct.perctPC[plasticbirds],na.rm = T)
sd(chem$Hct.perctPC[-plasticbirds],na.rm = T)
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

midHct<-barplot(Hct,col = rainbow(c(2,3)), plot = F)
midHct
barplot(Hct,col = rainbow(c(2,3)),ylim=range(0,45),yaxp=c(0, 180, 18),names.arg = c("with plastic","without plastic"),ylab = "%pcu Hct",las=1)
arrows(x0=midHct, y0=Hct-y.sdHct, x1=midHct, y1=Hct+y.sdHct, code=3, angle=90, length=0.1)

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

midHg<-barplot(Hg,col = rainbow(c(2,3)), plot = F)
midHg
barplot(Hg,col = rainbow(c(2,3)),ylim=range(0,15),yaxp=c(0, 180, 18),names.arg = c("with plastic","without plastic"),ylab = "Hb g/dL",las=1)
arrows(x0=midHg, y0=Hg-y.sdHg, x1=midHg, y1=Hg+y.sdHg, code=3, angle=90, length=0.1)

## Bar plot for AnGao mmol/L 
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

midAn<-barplot(An,col = rainbow(c(2,3)), plot = F)
midAn
barplot(An,col = rainbow(c(2,3)),ylim=range(0,25),yaxp=c(0, 180, 18),names.arg = c("with plastic","without plastic"),ylab = "Mmol AnGap",las=1)
arrows(x0=midAn, y0=An-y.sdAn, x1=midAn, y1=An+y.sdAn, code=3, angle=90, length=0.1)

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

midCr<-barplot(Cr,col = rainbow(c(2,3)), plot = F)
midCr
barplot(Cr,col = rainbow(c(2,3)),ylim=range(0,1),yaxp=c(0, 1, 1),names.arg = c("with plastic","without plastic"),ylab = "Crea mg/dL",las=1)
arrows(x0=midCr, y0=Cr-y.sdCr, x1=midCr, y1=Cr+y.sdCr, code=3, angle=90, length=0.1)
help(attach)


## Bar Plot for weight and plastic
#now find the mean and s.d. of iCammol/L for birds in these rows and in rows outside this set:
mean(MOR$Weight.24[plasticbirds],na.rm = T)
mean(MOR$Weight.24[-plasticbirds],na.rm = T)
sd(MOR$Weight.24[plasticbirds],na.rm = T)
sd(MOR$Weight.24[-plasticbirds],na.rm = T)
#t-test of differences in mean between these groups
W_pos<-(MOR$Weight.24[plasticbirds])
W_neg<-(MOR$Weight.24[-plasticbirds])
##y1 <- c(Na_pos, Na_neg)
t.test(iCa_pos,iCa_neg)
#now create midpoints for where you want the bars and error bars to go
meanplastIC<-mean(MOR$Weight.24[plasticbirds],na.rm = T)
meannoplastIC<-mean(MOR$Weight.24[-plasticbirds],na.rm = T)
W<-c(meanplastIC,meannoplastIC)
sdplastW<-sd(MOR$Weight.24[plasticbirds],na.rm = T)
sdnoplastW<-sd(MOR$Weight.24[-plasticbirds],na.rm = T)
y.sdW<-c(sdplastW,sdnoplastW)

midW<-barplot(W,col = rainbow(c(2,3)), plot = F)
midW
barplot(W,col = rainbow(c(2,3)),ylim=range(0,600),yaxp=c(0, 500, 10),names.arg = c("with plastic","without plastic"),ylab = "Weight",las=1)
arrows(x0=midW, y0=W-y.sdW, x1=midW, y1=W+y.sdW, code=3, angle=90, length=0.1)

