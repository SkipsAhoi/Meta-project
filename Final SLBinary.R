#FINAL SL BINARY 

#FINAL SLBinary MODELS
#R Packages - some unused.
MetaPacks<-c("MCMCglmm", "MuMIn", "ggplot2", "lattice", "cowplot", "mcmcplots")

library(kinship2) ; library(MCMCglmm) ; library(pedigree) ; library(pedigreemm) ; library(MasterBayes) ; library(MuMIn) ; library(AICcmodavg) ; library(sjstats); library(coda)
library(lme4);library(arm);library(sjPlot);library(glmulti);library(bblme);library(Rcpp);library(coefplot);library(ggplot2); library(multcomp); library(blmeco)
library(gridExtra);library(cowplot);library(ggplot2); library(lattice); library(coefplot2); library(mcmcplots); library(coda); library(mfx); library(spatialprobit); library(MASS); library(boot)

#The Data
#datafile setup and variable correction
d <- as.matrix(read.table("~/Meta-project/MCMC_GLMM.txt", header=TRUE, na.strings="NA", sep = "\t")); d <- as.data.frame(d); d[is.na(d)] <- NA;str(d)
d$SLBinary<-as.numeric(as.character(d$SLBinary)); d$AdjScore<-as.numeric(as.character(d$AdjScore)); d$TotalExp<-as.numeric(as.character(d$TotalExp))
d$Age<-as.numeric(as.character(d$Age)); d$SLBinary<-as.numeric(as.character(d$SLBinary)); d$CurExp<-as.numeric(as.character(d$CurExp)); str(d)
d$ID<- d$animal
d$Age.group<- as.numeric(as.character(d$Age.group))
#eliminate NA rows to streamline the dataset
d<-subset(d, !is.na(SLBinary)); str(d)
# Let's subset out all of the ghost conditions, which don't have social learning
D<-d; D<-subset(D,!(Context == "Ghost")) ; D$Context<-drop.levels(D$Context) ; str(D)
# And then let's remove individuals with unknown rearing histories
D<-subset(D,!(Rearing == "U")) ; D$Rearing<-drop.levels(D$Rearing) ; str(D)
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
a <- p2$animal; b <- p2$dam; c <- p2$sire; e <- d$animal
difs <- setdiff(b,a); length(difs); list(difs); difs <- setdiff(c,a);length(difs); list(difs); difs <- setdiff(e,a);length(difs); list(difs)
#for manual checking use:
#b %in% a; c %in% a; e %in% a; e %in% a; b %in% c; a %in% e

#The Priors
prior3 <- list(R = list(V=1, nu=0.2), G = list(G1 = list(V=1, nu=0.2), G2 = list(V=1, nu=0.2), G3 = list(V=1, nu=0.2)))

#The models
#Null model
m1b <- MCMCglmm(SLBinary ~ 1, random = ~animal + ID + study2, family = "categorical",
               prior = prior3, pedigree = p, verbose = F, data = D, nitt = 500000, burnin = 100000, thin = 250)
summary(m2b)
#Rearing1
m2b <- MCMCglmm(SLBinary ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID + study2, family = "categorical",
                prior = prior3, pedigree = p, data = D, nitt = 500000, verbose = T, burnin = 100000, thin = 250)
#Rearing 2
m3b<- MCMCglmm(SLBinary ~ 1 + sex + Rearing2 + Age + CurExp, random = ~animal + ID + study2, family = "categorical",
                prior = prior3, pedigree = p, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250)
summary(m3b)
m4b<- MCMCglmm(SLBinary ~ 1 + sex + Rearing2 + Age + CurExp, random = ~animal + ID + study2, family = "categorical",
               prior = prior3, pedigree = p, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250, singular.ok=T)
summary(m4b)
#Rearing 3
m5b<- MCMCglmm(SLBinary ~ 1 + sex + Rearing3 + Age + CurExp, random = ~animal + ID + study2, family = "categorical",
                prior = prior3, pedigree = p2, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250)
summary(m1b); summary(m2b);
summary(m3b); summary(m4b);  summary(m4.5b)
summary(m4b)
#In each case, fit new model with only significant predictors.
#Rearing1
m2.5b <- MCMCglmm(SLBinary ~ 1 + sex + Rearing, random = ~animal + ID + study2, family = "categorical",
                  prior = prior3, pedigree = p2, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250)
#Rearing2
m3.5b <- MCMCglmm(SLBinary ~ 1 + sex + Rearing2, random = ~animal + ID + study2, family = "categorical",
                  prior = prior3, pedigree = p2, data = D, nitt = 10000, verbose = F, burnin = 100000, thin = 250)
#Rearing3
m4.5b <- MCMCglmm(SLBinary ~ 1 + sex + Rearing3, random = ~animal + ID + study2, family = "categorical",
                prior = prior3, pedigree = p2, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250)
summary(m1b); summary(m2b); summary(m2.5b) ;summary(m3b); summary(m3.5b) ; summary(m4b);  summary(m4.5b)

#full model summary list
summary(c(m1b, m2, m2.5b, m3, m3.5b, m4, m4.5b))

#Check for Autocorrelation
autocorr.diag(m1b$Sol) ; autocorr.diag(m1b$VCV) ;
autocorr.diag(m2b$Sol) ; autocorr.diag(m2b$VCV) 
autocorr.diag(m3b$Sol) ; autocorr.diag(m3b$VCV) ;
autocorr.diag(m4b$Sol) ; autocorr.diag(m4b$VCV)
#No autocorrelation (lag 250 = <0.1), so let's check the mixing plots
plot(m1b$Sol); plot(m1b$VCV)
plot(m2b$Sol); plot(m2b$VCV)
plot(m2.5b$Sol); plot(m2.5b$VCV)
plot(m3b$Sol); plot(m3b$VCV)
plot(m3.5b$Sol); plot(m3.5b$VCV)
plot(m4b$Sol); plot(m4b$VCV)
plot(m4.5b$Sol); plot(m4.5b$VCV)
#Everything looks good! No need to run chains for longer or alter thinning interval.
# Now let's check the effective sample sizes
#Check effective sample sizes
effectiveSize(m1b$Sol);effectiveSize(m1b$VCV);effectiveSize(m2b$Sol);effectiveSize(m2b$VCV);effectiveSize(m2.5b$Sol);effectiveSize(m2.5b$VCV);effectiveSize(m3b$Sol);effectiveSize(m3b$VCV); effectiveSize(m3.5b$Sol);effectiveSize(m3.5b$VCV);effectiveSize(m4b$Sol) ; effectiveSize(m4b$VCV);effectiveSize(m4.5b$Sol) ; effectiveSize(m4.5b$VCV)
#Everything is close to the full sample size, and well over 1000, so I'm happy with this.

#phenotypic variance conditional on fixed effects
mVP<-m2.5b$VCV[,"animal"]+m2.5b$VCV[,"ID"]+m2.5b$VCV[,"units"]
#between individual variance
mIP<-m2.5b$VCV[,"animal"]+m2.5b$VCV[,"ID"]
#repeatability
rep<-posterior.mode(mIP/mVP); rep
#heritability
her<-posterior.mode(m2.5b$VCV[,"ID"]/mVP); her
her<-posterior.mode(m2.5b$VCV[,"animal"]/mVP); her
her<-posterior.mode(m2.5b$VCV[,"study2"]/mVP); her
