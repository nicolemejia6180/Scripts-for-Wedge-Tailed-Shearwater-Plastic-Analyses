##barplot using ggplot2##
#create dataframe from mean values calculated in "script_complete_figures.R##
##significant TCO2
dfco2 <- data.frame(Plastic=c("+", "-"), TCO2 =c(11.90909, 13.85714), sd=c(1.700267, 3.08488))
dfco2
#create barplot
TCO2 <- ggplot(data=dfco2, aes(x=Plastic, y=TCO2, fill=Plastic)) + 
  geom_bar(stat = "identity", color="black") + labs(y = "TCO2 (mmol)") + theme_classic() + scale_fill_grey(start=.5, end=1) + geom_errorbar(aes(x=Plastic, ymin=TCO2-sd, ymax=TCO2+sd), colour="black", alpha=0.9, width=0.4) + theme(legend.position="none") + theme(axis.title=element_text(size=7))
TCO2
#Nammol barplot
dfna <- data.frame(Plastic=c("+", "-"), NAMMOL =c(155.8182, 154.8571), sd=c(5.473905, 4.588567))
dfna
#create barplot
NAmmol <- ggplot(data=dfna, aes(x=Plastic, y=NAMMOL, fill=Plastic)) + 
  geom_bar(stat = "identity", color="black") + labs(y = "Sodium (mmol)") + theme_classic() + scale_fill_grey(start=.5, end=1) + geom_errorbar(aes(x=Plastic, ymin=NAMMOL-sd, ymax=NAMMOL+sd), colour="black", alpha=0.9, width=0.4) + theme(legend.position="none") + theme(axis.title=element_text(size=7))     
NAmmol                        

#Potassium barplot 
dfk <- data.frame(Plastic=c("+", "-"), K=c(3.272727, 3.314286), sd=c(0.7044017, 0.4943638))
dfk
Kmmol <- ggplot(data=dfk, aes(x=Plastic, y=K, fill=Plastic)) + geom_bar(stat = "identity", color="black") + labs(y = "Potassium (mmol)") + theme_classic() + scale_fill_grey(start=.5, end=1) + geom_errorbar(aes(x=Plastic, ymin=K-sd, ymax=K+sd), colour="black", alpha=0.9, width=0.4) + theme(legend.position="none") + theme(axis.title=element_text(size=7))                              
Kmmol

#Chloride
dfcl <- data.frame(Plastic=c("+", "-"), CL=c(132.2727, 128.6429), sd=c(4.900835, 5.583039))
dfcl
Clmmol <- ggplot(data=dfcl, aes(x=Plastic, y=CL, fill=Plastic)) + geom_bar(stat = "identity", color="black") + labs(y = "Chloride (mmol)") + theme_classic() + scale_fill_grey(start=.5, end=1) + geom_errorbar(aes(x=Plastic, ymin=CL-sd, ymax=CL+sd), colour="black", alpha=0.9, width=0.4) + theme(legend.position="none") + theme(axis.title=element_text(size=7)) 
Clmmol

#Ionized Calcium 
dfic <- data.frame(Plastic=c("+", "-"), IC=c(1.127273, 1.145), sd=c(0.09360458, 0.09920841))
dfic
ICmmol <- ggplot(data=dfic, aes(x=Plastic, y=IC, fill=Plastic)) + geom_bar(stat = "identity", color="black") + labs(y = "Ionized Calcium (mmol)") + theme_classic() + scale_fill_grey(start=.5, end=1) + geom_errorbar(aes(x=Plastic, ymin=IC-sd, ymax=IC+sd), colour="black", alpha=0.9, width=0.4) + theme(legend.position="none") + theme(axis.title=element_text(size=7))  
ICmmol

#Glucose
dfglu <- data.frame(Plastic=c("+", "-"), Glu=c(256, 261.3571), sd=c(24.62113, 33.69685))
dfglu
Glumgdl <- ggplot(data=dfglu, aes(x=Plastic, y=Glu, fill=Plastic)) + geom_bar(stat = "identity", color="black") + labs(y = "Glucose (mg/dL)") + theme_classic() + scale_fill_grey(start=.5, end=1) + geom_errorbar(aes(x=Plastic, ymin=Glu-sd, ymax=Glu+sd), colour="black", alpha=0.9, width=0.4) + theme(legend.position="none") + theme(axis.title=element_text(size=7)) 
Glumgdl

##BUN 
dfbun <- data.frame(Plastic=c("+", "-"), BUN=c(9.998182, 8.497857), sd=c(9.727693, 6.82659))
dfbun
BUNmgdl <- ggplot(data=dfbun, aes(x=Plastic, y=BUN, fill=Plastic)) + geom_bar(stat = "identity", color="black") + labs(y = "BUN (mg/dL)") + theme_classic() + scale_fill_grey(start=.5, end=1) + geom_errorbar(aes(x=Plastic, ymin=BUN-sd, ymax=BUN+sd), colour="black", alpha=0.9, width=0.4) + theme(legend.position="none") + theme(axis.title=element_text(size=7)) 
BUNmgdl

#Creatine 
dfcr <- data.frame(Plastic=c("+", "-"), CR=c(0.199, 0.2205), sd=c(0, 0.08044563))
dfcr
CR <- ggplot(data=dfcr, aes(x=Plastic, y=CR, fill=Plastic)) + geom_bar(stat = "identity", color="black") + labs(y = "Creatine (mg/dL)") + theme_classic() + scale_fill_grey(start=.5, end=1) + geom_errorbar(aes(x=Plastic, ymin=CR-sd, ymax=CR+sd), colour="black", alpha=0.9, width=0.4) + theme(legend.position="none") + theme(axis.title=element_text(size=7)) 
CR

##HCT
dfhct <- data.frame(Plastic=c("+", "-"), HCT=c(34.36364, 35.21429), sd=c(4.500505, 3.745327))
dfhct
HCTpcu <- ggplot(data=dfhct, aes(x=Plastic, y=HCT, fill=Plastic)) + geom_bar(stat = "identity", color="black") + labs(y = "Hematocrit (%)") + theme_classic() + scale_fill_grey(start=.5, end=1) + theme(legend.position="none") + geom_errorbar(aes(x=Plastic, ymin=HCT-sd, ymax=HCT+sd), colour="black", alpha=0.9, width=0.4) + theme(axis.title=element_text(size=7)) 
HCTpcu

#Hemoglobin
dfhg <- data.frame(Plastic=c("+", "-"), HG=c(11.7, 12.20714), sd=c(1.524467, 0.9110132))
dfhg
HG <- ggplot(data=dfhg, aes(x=Plastic, y=HG, fill=Plastic)) + geom_bar(stat = "identity", color="black") + labs(y = "Hemoglobin (g/dL)") + theme_classic() + scale_fill_grey(start=.5, end=1) + theme(legend.position="none") + geom_errorbar(aes(x=Plastic, ymin=HG-sd, ymax=HG+sd), colour="black", alpha=0.9, width=0.4) + theme(axis.title=element_text(size=7)) 
HG

#AnionGap
dfan <- data.frame(Plastic=c("+", "-"), AN=c(13.63636, 15.71429), sd=c(5.463931, 5.150483))
dfan
AN <- ggplot(data=dfan, aes(x=Plastic, y=AN, fill=Plastic)) + geom_bar(stat = "identity", color="black") + labs(y = "Anion Gap (mmol)") + theme_classic() + scale_fill_grey(start=.5, end=1) + theme(legend.position="none") + geom_errorbar(aes(x=Plastic, ymin=AN-sd, ymax=AN+sd), colour="black", alpha=0.9, width=0.4) + theme(axis.title=element_text(size=7)) 
AN

##Body weight
dfbw <- data.frame(Plastic=c("+", "-"), BW=c(377.0833, 407.1250), sd=c(46.1806236423892, 41.05865922588))
dfbw
BW <- ggplot(data=dfbw, aes(x=Plastic, y=BW, fill=Plastic)) + geom_bar(stat = "identity", color="black") + labs(y = "Body weight (g)") + theme_classic() + scale_fill_grey(start=.5, end=1) + theme(legend.position="none") + geom_errorbar(aes(x=Plastic, ymin=BW-sd, ymax=BW+sd), colour="black", alpha=0.9, width=0.4) + theme(axis.title=element_text(size=7)) 
BW
##together?
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,3)))
print(NAmmol, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(Kmmol, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(Clmmol, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(ICmmol, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(TCO2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(Glumgdl, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
print(BUNmgdl, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(CR, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
print(HCTpcu, vp = viewport(layout.pos.row = 3, layout.pos.col = 3))
print(HG, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
print(AN, vp = viewport(layout.pos.row = 4, layout.pos.col = 2))
print(BW, vp = viewport(layout.pos.row = 4, layout.pos.col = 3))



