############             Library            #####
library(readr)
library(agricolae)
library(readr)
library(agricolae)
library(lme4)
library(lmerTest)
library(glm2)
library(lme4)
library(lmerTest)
library(glm2)
library(ggplot2)
library(stats)
library(ade4)
library(gclus)
library(ape)
library(vegan)
library(grid)
library(gridExtra)
library(mosaic)
library(tigerstats)
library(multcomp)
library(ggrepel)
library(grid)
library(gridExtra)
library(car)
library(ggrepel)
library(emmeans)
library(ggrepel)
library(emmeans)
library(PerformanceAnalytics)
library(openxlsx)
library(EnvStats)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(ggridges)
library(dplyr)
library(gridExtra)
library(Rmisc)
theme_set(theme_pubr())
library(emmeans)
library(ggrepel)
library(emmeans)
library(PerformanceAnalytics)
library(openxlsx)
library(EnvStats)
library("ggpubr")

############################## Current Bio script #############################
###### COrrelations NUCLEI AND SPORE SIZE
Nuclearnumber <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/DATA/confocal/Nuclearnumber.csv")

scP<-ggscatter(Nuclearnumber, x = "SporeSize", y = "Nuclei", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "Spore Size (um) in vitro", ylab = "Number of Nuclei")
scP


scP<-ggscatter(Nuclearnumber, x = "SporeSize", y = "Nuclei", color="Karyosis",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "Spore Size (um)", ylab = "Number of Nuclei")
scP

SCP2<- ggscatter(Nuclearnumber, x = "SporeSize", y = "Nuclei",
                 color = "Karyosis", shape = "Karyosis",
                 ellipse = TRUE, mean.point = TRUE,
                 star.plot = TRUE,xlab = "Spore Size (um)", ylab = "Number of Nuclei")
SCP2
library(gridExtra)
grid.arrange(scP, SCP2, nrow = 1)


NuclearnumberHOM <- Nuclearnumber[1:125,]
NuclearnumberDIK <- Nuclearnumber[126:173,]

scP<-ggscatter(NuclearnumberDIK, x = "SporeSize", y = "Nuclei", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "Spore Size (um) in vitro", ylab = "Number of Nuclei")
scP

scP<-ggscatter(NuclearnumberHOM, x = "SporeSize", y = "Nuclei", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "Spore Size (um) in vitro", ylab = "Number of Nuclei")
scP


####### ANCOVA
mod1 <- aov(SporeSize~Nuclei*Karyosis, data=Nuclearnumber)
summary(mod1) 

#Df Sum Sq Mean Sq F value Pr(>F)    
# Nuclei            1  55468   55468 222.280 <2e-16 ***
# Karyosis          1   1137    1137   4.556 0.0342 *  
# Nuclei:Karyosis   1     58      58   0.233 0.6298    
# Residuals       169  42173     250  

#### Since interaction is not significant then the slopes of the model 
# are not different and the correlation of nuclei~spore size is similar 
# for DIK and HOM

reg1 <- lm(Nuclei~SporeSize, data=NuclearnumberDIK); summary(reg1)
reg2 <- lm(Nuclei~SporeSize, data=NuclearnumberHOM); summary(reg1)

plot(SporeSize~Nuclei, data=Nuclearnumber, type='n')
points(NuclearnumberDIK$SporeSize,NuclearnumberDIK$Nuclei, pch=20)
points(NuclearnumberHOM$SporeSize,NuclearnumberHOM$Nuclei, pch=1)
abline(reg1, lty=1)
abline(reg2, lty=2)
legend("bottomright", c("DIK","HOM"), lty=c(1,2), pch=c(20,1) )


mod1b <- aov(Nuclei~Isolate, data=Nuclearnumber)
summary(mod1b)
tukey.test <- TukeyHSD(mod1b)
tukey.test
TukeyHSD(mod1b)
Tukey = HSD.test(mod1b, "Isolate")
Tukey

a1<-ggplot(data = Nuclearnumber, aes( x = Karyosis, y = Nuclei, fill = Isolate ) ) +    # print bar chart
  geom_boxplot()+
  xlab("") + ylab("Number of nuclei")
a1

######## Number of nuclei confocal vs ddPCR ##############
Nuclearnumber_A4 <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Nuclear data/Nuclearnumber_A4.csv")
Nuclearnumber_A5 <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Nuclear data/Nuclearnumber_A5.csv")
Nuclearnumber_SL1 <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Nuclear data/Nuclearnumber_SL1.csv")
Nuclearnumber_G1 <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Nuclear data/Nuclearnumber_G1.csv")

par(mfrow=c(1,2))
boxplot(Nuclei~Isolate,data=Nuclearnumber_A4, ylab= "Number ot nuclei ",
        xlab=" Method of quantification", names=c("Confocal", "ddPCR"),ylim=c(0, 800 ), frame = FALSE, col=terrain.colors(3) )
boxplot(Nuclei~Strain,data=Nuclearnumber_A4, names=c("A4"),tittle = "Compiled data",ylab= "", ylim=c(0, 800 ), xlab="Compiled data",col=terrain.colors(1) )

boxplot(Nuclei~Isolate,data=Nuclearnumber_A5, ylab= "Number ot nuclei ",
        xlab=" Method of quantification", names=c("Confocal", "ddPCR"),ylim=c(0, 1800 ), frame = FALSE, col=terrain.colors(3) )
boxplot(Nuclei~Strain,data=Nuclearnumber_A5, names=c("A4"),tittle = "Compiled data",ylab= "", ylim=c(0, 1800 ), xlab="Compiled data",col=terrain.colors(1) )

boxplot(Nuclei~Isolate,data=Nuclearnumber_SL1, ylab= "Number ot nuclei ",
        xlab=" Method of quantification", names=c("Confocal", "ddPCR"),ylim=c(0, 700 ), frame = FALSE, col=terrain.colors(3) )
boxplot(Nuclei~Strain,data=Nuclearnumber_SL1, names=c("A4"),tittle = "Compiled data",ylab= "", ylim=c(0, 700 ), xlab="Compiled data",col=terrain.colors(1) )


boxplot(Nuclei~Isolate,data=Nuclearnumber_G1, ylab= "Number ot nuclei ",
        xlab=" Method of quantification", names=c("Confocal", "ddPCR"),ylim=c(0, 800 ), frame = FALSE, col=terrain.colors(3) )
boxplot(Nuclei~Strain,data=Nuclearnumber_G1, names=c("A4"),tittle = "Compiled data",ylab= "", ylim=c(0, 800 ), xlab="Compiled data",col=terrain.colors(1) )


kruskal.test(Nuclei~Isolate, data = Nuclearnumber_A4)
#OUTCOME:Kruskal-Wallis chi-squared = 0.3493, df = 1, p-value = 0.5545
kruskal.test(Nuclei~Isolate, data = Nuclearnumber_A5)
#Kruskal-Wallis chi-squared = 1.0836, df = 1, p-value = 0.2979
kruskal.test(Nuclei~Isolate, data = Nuclearnumber_SL1)
#Kruskal-Wallis chi-squared = 1.5536, df = 1, p-value = 0.2126
kruskal.test(Nuclei~Isolate, data = Nuclearnumber_G1)
#Kruskal-Wallis chi-squared = 2.572, df = 1, p-value = 0.1088


### Ratio examination
#### eXAMPLE FIGURE  (G1), HOW TO READ DDpcr AND ANALYZE DATA
Ratio_data <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Ratio_data.csv")
boxplot(log1p(Value)~Target, outline = FALSE, data=Ratio_data)
boxplot(log1p(Value)~Target + Culture, outline = FALSE, frame = FALSE,xlab="",ylab= "copies per ??l (log)", col = c("chartreuse4", "blue", "lightsalmon"), data=Ratio_data)
abline(v=3.5, col="blue",lwd=3, lty=2)

Ratio_data_G1 <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Ratio_data_G1.csv")
boxplot(Value ~  Culture, outline = FALSE, frame = FALSE,xlab="",ylab= "Ratio (MAT-1 : MAT-5)", col = c("lightsalmon"), data=Ratio_data_G1)
leveneTest(Value ~ Culture, data = Ratio_data_G1)
ANOVA<-aov(Value~Culture, data = Ratio_data_G1)
summary(ANOVA)

### ALL RATIOS WITH CAROT (old cultivar)
Ratio_data_all_carrot <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Ratio_data_all_carrot_FINAL.csv")
boxplot(Value~Strain, outline = FALSE,frame = FALSE, col = c("lightsalmon"), xlab="",ylab= "Ratio (MAT-A : MAT-B)", data=Ratio_data_all_carrot)
abline(h=1, col="blue")
leveneTest(Value ~ Strain, data = Ratio_data_all_carrot)
#### Non parametric test due to inequality of n ###
kruskal.test(Value~Strain, data = Ratio_data_all_carrot)
pairwise.wilcox.test(Ratio_data_all_carrot$Value, Ratio_data_all_carrot$Strain,
                     p.adjust.method = "BH")

#### Subset_data to test variation between sub-cultures #####
boxplot(Value~Culture+Strain,frame = FALSE, col = c("lightsalmon","lightsalmon4"), xlab="",ylab= "Ratio (MAT-A : MAT-B)", outline = FALSE, data=Ratio_data_all_carrot)
abline(h=1, col="blue")
G1_subcultures <- Ratio_data_all_carrot[1:87,]
A5_subcultures <- Ratio_data_all_carrot[88:220,]
SL1_subcultures <- Ratio_data_all_carrot[221:298,]
A4_subcultures <- Ratio_data_all_carrot[299:373,]

ANOVA<-aov(Value~Culture, data = G1_subcultures)
summary
#             Df  Sum Sq  Mean Sq F value Pr(>F)
#Culture      1 0.00171 0.001714   0.529  0.469
#Residuals   85 0.27527 0.003238
leveneTest(Value ~ Culture, data = G1_subcultures)

ANOVA<-aov(Value~Culture, data = A5_subcultures)
summary(ANOVA)
#             Df Sum Sq  Mean Sq F value Pr(>F)
#Culture       1 0.0056 0.005564    0.46  0.499
#Residuals   131 1.5860 0.012107 
leveneTest(Value ~ Culture, data = A5_subcultures)

ANOVA<-aov(log1p(Value)~Culture, data = SL1_subcultures)
summary(ANOVA)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#Culture      1  16.26  16.263   36.52 5.22e-08 ***
# Residuals   76  33.85   0.445 
leveneTest(log1p(Value) ~ Culture, data = SL1_subcultures)

ANOVA<-aov(Value~Culture, data = A4_subcultures)
summary(ANOVA)
#           Df Sum Sq Mean Sq F value Pr(>F)
#Culture      1  0.006  0.0057   0.051  0.822
#Residuals   73  8.190  0.1122 
leveneTest(Value ~ Culture, data = A4_subcultures)


### aLL RATIOS WITH CAROT CV. p68
Ratio_all_P68 <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Ratio_all_P68_FINAL.csv")
boxplot(Value~Strain, outline = FALSE,frame = FALSE, ylim=c(0.3, 1.2 ),col = c("lightsalmon"), xlab="",ylab= "Ratio (MAT-A : MAT-B)", data=Ratio_all_P68)
abline(h=1, col="blue")
ANOVA<-aov(Value~Strain, data = Ratio_all_P68)
summary(ANOVA)
tukey.test <- TukeyHSD(ANOVA)
tukey.test
TukeyHSD(ANOVA)
Tukey = HSD.test(ANOVA, "Strain")
Tukey
shapiro.test(resid(ANOVA))

#### Subset_data to test variation between sub-cultures #####
boxplot(Value~Culture+Strain,frame = FALSE, ylim=c(0.3, 1.2 ),col = c("lightsalmon","lightsalmon4"), xlab="",ylab= "Ratio (MAT-A : MAT-B)", outline = FALSE, data=Ratio_all_P68)
abline(h=1, col="blue")

G1_subcultures <- Ratio_all_P68[85:132,]
A5_subcultures <- Ratio_all_P68[1:44,]
SL1_subcultures <- Ratio_all_P68[45:84,]
A4_subcultures <- Ratio_all_P68[133:180,]

ANOVA<-aov(Value~Culture, data = G1_subcultures)
summary(ANOVA)
#             Df Sum Sq  Mean Sq F value Pr(>F)
#Culture      1 0.0111 0.011078   1.428  0.238
#Residuals   46 0.3568 0.007756
ANOVA<-aov(Value~Culture, data = A5_subcultures)
summary(ANOVA)
#            Df Sum Sq  Mean Sq F value Pr(>F)
#Culture      1 0.0000 0.000003       0  0.983
#Residuals   42 0.3235 0.007703 
ANOVA<-aov(Value~Culture, data = SL1_subcultures)
summary(ANOVA)
#            Df  Sum Sq  Mean Sq F value Pr(>F)
#Culture      1 0.00528 0.005280   2.029  0.163
#Residuals   38 0.09891 0.002603
ANOVA<-aov(Value~Culture, data = A4_subcultures)
summary(ANOVA)
#           Df  Sum Sq  Mean Sq F value Pr(>F)
#Culture      1 0.00137 0.001370   0.231  0.633
#Residuals   46 0.27346 0.005945

##### rATIO MEDICAGO
Ratio_all_medicago <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Ratio_all_medicago_FINAL.csv")
boxplot(Value~Strain, outline = FALSE,frame = FALSE, ylim=c(0.3, 1.6),col = c("lightsalmon"), xlab="",ylab= "Ratio (MAT-A : MAT-B)", data=Ratio_all_medicago)
abline(h=1, col="blue")
boxplot(Value~Culture+Strain,frame = FALSE, ylim=c(0.3, 1.6 ),col = c("lightsalmon","lightsalmon4"), xlab="",ylab= "Ratio (MAT-A : MAT-B)", outline = FALSE, data=Ratio_all_medicago)
abline(h=1, col="blue")
ANOVA<-aov(log1p(Value)~Strain, data = Ratio_all_medicago)
summary(ANOVA)
tukey.test <- TukeyHSD(ANOVA)
tukey.test
TukeyHSD(ANOVA)
Tukey = HSD.test(ANOVA, "Strain")
Tukey
shapiro.test(resid(ANOVA))
### Non parametric##
kruskal.test(Value~Strain, data = Ratio_all_medicago)
pairwise.wilcox.test(Ratio_all_medicago$Value, Ratio_all_medicago$Strain,
                     p.adjust.method = "BH")



#### Subset_data to test variation between sub-cultures #####
G1_subcultures <- Ratio_all_medicago[91:136,]
A5_subcultures <- Ratio_all_medicago[45:90,]
SL1_subcultures <- Ratio_all_medicago[137:180,]
A4_subcultures <- Ratio_all_medicago[1:44,]

ANOVA<-aov(Value~Culture, data = G1_subcultures)
summary(ANOVA)
##            Df Sum Sq Mean Sq F value Pr(>F)
#Culture      1 0.0618 0.06179   1.825  0.184
#Residuals   44 1.4900 0.03386 
ANOVA<-aov(Value~Culture, data = A5_subcultures)
summary(ANOVA)
#             Df Sum Sq  Mean Sq F value Pr(>F)
#Culture      1 0.0044 0.004402   0.667  0.419
#Residuals   44 0.2906 0.006604 
ANOVA<-aov(Value~Culture, data = SL1_subcultures)
summary(ANOVA)
#             Df  Sum Sq  Mean Sq F value Pr(>F)
#Culture      1 0.00359 0.003590   0.983  0.327
#Residuals   40 0.14604 0.003651
ANOVA<-aov(Value~Culture, data = A4_subcultures)
summary(ANOVA)
#             Df Sum Sq  Mean Sq F value Pr(>F)
#Culture      1 0.0006 0.000582   0.038  0.846
#Residuals   42 0.6358 0.015137#


##### rATIO CHICORY
Ratio_data_all_chicory <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Ratio_data_all_chicory_FINAL.csv")
boxplot(Value~Strain, outline = FALSE,frame = FALSE, ylim=c(0.1, 2.3),col = c("lightsalmon"), xlab="",ylab= "Ratio (MAT-A : MAT-B)", data=Ratio_data_all_chicory)
abline(h=1, col="blue")
boxplot(Value~Culture+Strain,frame = FALSE, ylim=c(0.1, 2.3 ),col = c("lightsalmon","lightsalmon4"), xlab="",ylab= "Ratio (MAT-A : MAT-B)", outline = FALSE, data=Ratio_data_all_chicory)
abline(h=1, col="blue")
ANOVA<-aov(log1p(Value)~Strain, data = Ratio_data_all_chicory)
summary(ANOVA)
tukey.test <- TukeyHSD(ANOVA)
tukey.test
TukeyHSD(ANOVA)
Tukey = HSD.test(ANOVA, "Strain")
Tukey
shapiro.test(resid(ANOVA))
### Non parametric##
kruskal.test(Value~Strain, data = Ratio_data_all_chicory)
pairwise.wilcox.test(Ratio_data_all_chicory$Value, Ratio_data_all_chicory$Strain,
                     p.adjust.method = "BH")


#### Subset_data to test variation between sub-cultures #####
G1_subcultures <- Ratio_data_all_chicory[45:85,]
SL1_subcultures <- Ratio_data_all_chicory[86:134,]
A4_subcultures <- Ratio_data_all_chicory[1:44,]

ANOVA<-aov(Value~Culture, data = G1_subcultures)
summary(ANOVA)
#Df  Sum Sq  Mean Sq F value Pr(>F)
#Culture      1 0.00413 0.004128   0.622  0.435
#Residuals   39 0.25877 0.006635
ANOVA<-aov(Value~Culture, data = SL1_subcultures)
summary(ANOVA)
#Df Sum Sq Mean Sq F value Pr(>F)
#Culture      1   1.92   1.922   0.557  0.459
#Residuals   47 162.18   3.451

ANOVA<-aov(Value~Culture, data = A4_subcultures)
summary(ANOVA)
#Df Sum Sq Mean Sq F value Pr(>F)
#Culture      1    516   515.5   0.986  0.326
#Residuals   42  21962   522.9  

##### Ratio Nicotiana
Ratio_nicotiana <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Ratio_nicotiana.csv")
boxplot(Value~Strain, outline = FALSE,frame = FALSE, ylim=c(0.1, 2.3),col = c("lightsalmon"), xlab="",ylab= "Ratio (MAT-A : MAT-B)", data=Ratio_nicotiana)
abline(h=1, col="blue")
boxplot(Value~Culture+Strain,frame = FALSE, ylim=c(0.1, 2.3 ),col = c("lightsalmon","lightsalmon4"), xlab="",ylab= "Ratio (MAT-A : MAT-B)", outline = FALSE, data=Ratio_nicotiana)
abline(h=1, col="blue")
ANOVA<-aov(log1p(Value)~Culture, data = Ratio_nicotiana)
summary(ANOVA)
leveneTest(Value ~ Culture, data = Ratio_nicotiana)
shapiro.test(resid(ANOVA))


########### All Hosts (Figure 3 alternative)
Ratio_G1 <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Ratio_G1_FINAL.csv")
Ratio_A4 <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Ratio_A4_FINAL.csv")
Ratio_A5 <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Ratio_A5_FINAL.csv")
Ratio_SL1 <- read.csv("~/PDF/PDF_Experiments/Nuclear ratios/Final manuscript/Revised/CB/revisions/Data/Ratio_SL1_FINAL.csv")


par(mfrow=c(1,2))
boxplot(Value~Host, main = "G1", outline = FALSE,frame = FALSE, ylim=c(0.1, 1.2),col = c("lightsalmon"), xlab="",ylab= "Ratio (MAT-1 : MAT-5)", data=Ratio_G1)
abline(h=1, col="blue")
ANOVA<-aov(Value~Host, data = Ratio_G1)
summary(ANOVA)
tukey.test <- TukeyHSD(ANOVA)
tukey.test
TukeyHSD(ANOVA)
Tukey = HSD.test(ANOVA, "Host")
Tukey
shapiro.test(resid(ANOVA))
###non parametric
kruskal.test(Value~Host, data = Ratio_G1)
pairwise.wilcox.test(Ratio_G1$Value, Ratio_G1$Host,
                     p.adjust.method = "BH")


boxplot(Value~Host,main = "A4", outline = FALSE,frame = FALSE, ylim=c(0.1, 2.3),col = c("lightsalmon"), xlab="",ylab= "Ratio (MAT-1 : MAT-2)", data=Ratio_A4)
abline(h=1, col="blue")
ANOVA<-aov(log1p(Value)~Host, data = Ratio_A4)
summary(ANOVA)
tukey.test <- TukeyHSD(ANOVA)
tukey.test
TukeyHSD(ANOVA)
Tukey = HSD.test(ANOVA, "Host")
Tukey
shapiro.test(resid(ANOVA))
###non parametric
kruskal.test(Value~Host, data = Ratio_A4)
pairwise.wilcox.test(Ratio_A4$Value, Ratio_A4$Host,
                     p.adjust.method = "BH")



boxplot(Value~Host, main = "A5",outline = FALSE,frame = FALSE, ylim=c(0.1, 1.5),col = c("lightsalmon"), xlab="",ylab= "Ratio (MAT-6 : MAT-3)", data=Ratio_A5)
abline(h=1, col="blue")
ANOVA<-aov(log1p(Value)~Host, data = Ratio_A5)
summary(ANOVA)
tukey.test <- TukeyHSD(ANOVA)
tukey.test
TukeyHSD(ANOVA)
Tukey = HSD.test(ANOVA, "Host")
Tukey
shapiro.test(resid(ANOVA))
###non parametric
kruskal.test(Value~Host, data = Ratio_A5)
pairwise.wilcox.test(Ratio_A5$Value, Ratio_A5$Host,
                     p.adjust.method = "BH")

boxplot(Value~Host,main = "SL1", outline = FALSE,frame = FALSE, ylim=c(0.1, 2.9),col = c("lightsalmon"), xlab="",ylab= "Ratio (MAT-1 : MAT-5)", data=Ratio_SL1)
abline(h=1, col="blue")
ANOVA<-aov(Value~Host, data = Ratio_SL1)
summary(ANOVA)
tukey.test <- TukeyHSD(ANOVA)
tukey.test
TukeyHSD(ANOVA)
Tukey = HSD.test(ANOVA, "Host")
Tukey
shapiro.test(resid(ANOVA))
###non parametric
kruskal.test(Value~Host, data = Ratio_SL1)
pairwise.wilcox.test(Ratio_SL1$Value, Ratio_SL1$Host,
                     p.adjust.method = "BH")