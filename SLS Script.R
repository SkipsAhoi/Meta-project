#FINAL RAWSCORE MODELS
#R Packages - some unused.
library(kinship2) ; library(MCMCglmm) ; library(pedigree) ; library(pedigreemm) ; library(MasterBayes) ; library(MuMIn) ; library(AICcmodavg) ; library(sjstats); library(coda)
library(lme4);library(arm);library(sjPlot);library(glmulti);library(bblme);library(Rcpp);library(coefplot);library(ggplot2); library(multcomp); library(blmeco)
library(gridExtra);library(cowplot);library(ggplot2); library(lattice); library(coefplot2); library(mcmcplots); library(coda); library(mfx); library(spatialprobit); library(MASS); library(boot)
library(broom); library(tidyr)

#The Data
#datafile setup and variable correction
d <- as.matrix(read.table("~/Meta-project/MCMC_GLMM.txt", header=TRUE, na.strings="NA", sep = "\t")); d <- as.data.frame(d); d[is.na(d)] <- NA;str(d)
d$RawScore<-as.numeric(as.character(d$RawScore)); d$AdjScore<-as.numeric(as.character(d$AdjScore)); d$TotalExp<-as.numeric(as.character(d$TotalExp))
d$Age<-as.numeric(as.character(d$Age));  d$CurExp<-as.numeric(as.character(d$CurExp)); str(d)
d$ID<- d$animal
#d$SLBinary<-as.numeric(as.character(d$SLBinary))
d$Age.group<- as.numeric(as.character(d$Age.group))
#eliminate NA rows to streamline the dataset
d<-subset(d, !is.na(RawScore)); str(d)
# Let's subset out all of the ghost conditions, which don't have social learning
D<-d; D<-subset(D,!(Context == "Ghost")) ; D$Context<-drop.levels(D$Context) ; str(D)
# And then let's remove individuals with unknown rearing histories
D<-subset(D,!(Rearing == "U")) ; D$Rearing<-drop.levels(D$Rearing) ; str(D)
D$ID<-drop.levels(D$ID) ; str(D)
D$animal<-drop.levels(D$animal) ; str(D)
#drop individuals of age 1 
#D<-subset(D,!(Age.group == "1")) ; D$Age.group<-drop.levels(D$Age.group) ; str(D)
#subsetting individuals according to participation. Not sure about this.
#D<-subset(D, TotalExp > 2) ; D$TotalExp<-drop.levels(D$TotalExp) ; str(D)
str(D)

#This is our pedigree file. One column for animals, one for sires, one for dams. All sires and dams must be in animal column.
p <- read.csv("~/Meta-project/kinship2.CSV", header = T, na.strings="NA"); p<- as.data.frame(p); p[is.na(p)] <- NA; str(p)
#This command orders the pedigree file so that sires+dams come first in the list. Necessary for MCMCglmm to read it properly.
p<-orderPed(p);p[is.na(p)] <- NA; str(p)
# Let's check to make sure all Sires and Dams are in the animal list
a <- p$animal; b <- p$dam; c <- p$sire; e <- d$animal
difs <- setdiff(b,a); length(difs); list(difs); difs <- setdiff(c,a);length(difs); list(difs); difs <- setdiff(e,a);length(difs); list(difs)
#for manual checking use:
#b %in% a; c %in% a; e %in% a; e %in% a; b %in% c; a %in% e

#The Priors. Weakly informative.
prior3 <- list(R = list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.02), G3 = list(V=1, nu=0.02)))

#The models
#Null model
m1 <- MCMCglmm(RawScore ~ 1, random = ~animal + ID + study2, family = "gaussian",
               prior = prior3, pedigree = p, verbose = F, data = D, nitt = 250000, burnin = 10000, thin = 250)
#Rearing1
m2 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID + study2, family = "gaussian",
                prior = prior3, pedigree = p, data = D, nitt = 250000, verbose = F, burnin = 50000, thin = 250)
#Rearing 2
m3 <- MCMCglmm(RawScore ~ 1 + sex + Rearing2 + Age + CurExp, random = ~animal + ID + study2, family = "gaussian",
                prior = prior3, pedigree = p, data = D, nitt = 250000, verbose = F, burnin = 50000, thin = 250)
m3sing <- MCMCglmm(RawScore ~ 1 + sex + Rearing2 + Age + CurExp, random = ~animal + ID + study2, family = "gaussian",
                   prior = prior3, pedigree = p, data = D, nitt = 250000, verbose = F, burnin = 50000, thin = 250, singular.ok=T)
summary(m3sing)
#Rearing 3
m4 <- MCMCglmm(RawScore ~ 1 + sex + Rearing3 + Age + CurExp, random = ~animal + ID + study2, family = "gaussian",
                prior = prior3, pedigree = p, data = D, nitt = 250000, verbose = F, burnin = 50000, thin = 250, singular.ok=T)
summary(m4)
summary(m1); summary(m2); summary(m3); summary(m4)
summary(m3)

#full model summary list
#summary(m1); summary(m2); summary(m2.5) ; summary(m3); summary(m3.5) ; summary(m4) ;  summary(m4.5)
#mlistDIC<-list(m1$DIC, m2$DIC,  m2.5$DIC,  m3$DIC,  m3.5$DIC,  m4$DIC,  m4.5$DIC); names(mlistDIC)<-paste(1:7) ;  mlistDIC
summary(m1); summary(m2); summary(m3); summary(m4)
mlistDIC<-list(m1$DIC, m2$DIC,  m3$DIC,  m4$DIC); names(mlistDIC)<-paste(1:4) ;  mlistDIC
#Check for Autocorrelation
autocorr.diag(m1$Sol) ; autocorr.diag(m1$VCV) ;autocorr.diag(m2$Sol) ; autocorr.diag(m2$VCV) ;autocorr.diag(m3$Sol) ; autocorr.diag(m3$VCV) ; autocorr.diag(m4$Sol) ; autocorr.diag(m4$VCV)
#No autocorrelation (lag 250 = <0.1), so let's check the mixing plots
plot(m1$Sol); plot(m1$VCV)
plot(m2$Sol); plot(m2$VCV)
plot(m3$Sol); plot(m3$VCV)
plot(m4$Sol); plot(m4$VCV)
#Everything looks good! No need to run chains for longer or alter thinning interval.
# Now let's check the effective sample sizes
#Check effective sample sizes
effectiveSize(m1$Sol);effectiveSize(m1$VCV);effectiveSize(m2$Sol);effectiveSize(m2$VCV);effectiveSize(m2.5$Sol);effectiveSize(m2.5$VCV);effectiveSize(m3$Sol);effectiveSize(m3$VCV); effectiveSize(m3.5$Sol);effectiveSize(m3.5$VCV);effectiveSize(m4$Sol) ; effectiveSize(m4$VCV);effectiveSize(m4.5$Sol) ; effectiveSize(m4.5$VCV)
#Everything is close to the full sample size, and well over 1000, so I'm happy with this.

#phenotypic variance conditional on fixed effects
#heritability
herit1 <- m1$VCV[, "animal"]/(m1$VCV[, "animal"] + m1$VCV[, "units"]); effectiveSize(herit1); mean(herit1); HPDinterval(herit1) 
herit2 <- m2$VCV[, "animal"]/(m2$VCV[, "animal"] + m2$VCV[, "units"]); effectiveSize(herit2); mean(herit2); HPDinterval(herit2) 
herit3 <- m3$VCV[, "animal"]/(m3$VCV[, "animal"] + m3$VCV[, "units"]); effectiveSize(herit3); mean(herit3); HPDinterval(herit3) 
herit4 <- m4$VCV[, "animal"]/(m4$VCV[, "animal"] + m4$VCV[, "units"]); effectiveSize(herit4); mean(herit4); HPDinterval(herit4)

#rep
rep1 <- m1$VCV[, "ID"]/(m1$VCV[, "ID"] + m1$VCV[, "units"]); effectiveSize(rep1); mean(rep1); HPDinterval(rep1) 
rep2 <- m2$VCV[, "ID"]/(m2$VCV[, "ID"] + m2$VCV[, "units"]); effectiveSize(rep2); mean(rep2); HPDinterval(rep2) 
rep3 <- m3$VCV[, "ID"]/(m3$VCV[, "ID"] + m3$VCV[, "units"]); effectiveSize(rep3); mean(rep3); HPDinterval(rep3) 
rep4 <- m4$VCV[, "ID"]/(m4$VCV[, "ID"] + m4$VCV[, "units"]); effectiveSize(rep4); mean(rep4); HPDinterval(rep4)



newdfr = expand.grid(sex = c("Female", "Male"),
                     Rearing = c("W", "C", "N")
)
p.mod<-predict(m2)

ggplot(D, aes(x = sex, y = RawScore, color = Rearing)) +
  ylab("Score") + xlab("Sex") +
  geom_point() 
geom_line(data = newdfr, aes(x = sex, y = p.mod)) 
+
  facet_wrap(~ Subject, ncol = 3) +
  scale_x_continuous(breaks= seq(30, 130, by = 30))

plot(predict(m2),D$RawScore,
     xlab="predicted",ylab="actual")
abline(a=0,b=1)
par(mfrow = c(2, 2))
plot(m1)
