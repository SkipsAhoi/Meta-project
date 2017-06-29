#FINAL meta-anlaysis R scripts

#R Packages - some unused.

library(kinship2) ; library(MCMCglmm) ; library(pedigree) ; library(pedigreemm) ; library(MasterBayes) ; library(MuMIn) ; library(AICcmodavg) ; library(sjstats); library(coda)
library(lme4);library(arm);library(sjPlot);library(glmulti);library(bblme);library(Rcpp);library(coefplot);library(ggplot2); library(multcomp); library(blmeco)
library(gridExtra);library(cowplot);library(ggplot2); library(lattice); library(coefplot2); library(mcmcplots); library(coda); library(mfx); library(spatialprobit); library(MASS); library(boot)

pkgs_CRAN <- c("lme4","MCMCglmm","blme",
               "pbkrtest","coda","aods3","bbmle","ggplot2",
               "reshape2","plyr","numDeriv","Hmisc",
               "plotMCMC","gridExtra","R2admb")
install.packages(pkgs_CRAN)
rr <- "http://www.math.mcmaster.ca/bolker/R"
install.packages(c("glmmADMB","coefplot2"),type="source",
                 repos=rr)
library("coefplot2"); library("reshape2"); library("biomod2")

#The Data
#datafile setup and variable correction
d <- as.matrix(read.table("~/Meta-project/MCMC_GLMM.txt", header=TRUE, na.strings="NA", sep = "\t")); d <- as.data.frame(d); d[is.na(d)] <- NA;str(d)
d$RawScore<-as.numeric(as.character(d$RawScore)); d$AdjScore<-as.numeric(as.character(d$AdjScore)); d$TotalExp<-as.numeric(as.character(d$TotalExp))
d$Age<-as.numeric(as.character(d$Age)); d$SLBinary<-as.numeric(as.character(d$SLBinary)); d$CurExp<-as.numeric(as.character(d$CurExp)); str(d)
d$ID<- d$animal
d$Age.group<- as.numeric(as.character(d$Age.group))
#eliminate NA rows to streamline the dataset
d<-subset(d, !is.na(RawScore)); str(d)
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
p2<-orderPed(p);p2[is.na(p2)] <- NA; str(p2)
# Let's check to make sure all Sires and Dams are in the animal list
a <- p2$animal; b <- p2$dam; c <- p2$sire; e <- d$animal
difs <- setdiff(b,a); length(difs); list(difs); difs <- setdiff(c,a);length(difs); list(difs); difs <- setdiff(e,a);length(difs); list(difs)
#for manual checking use:
#b %in% a; c %in% a; e %in% a; e %in% a; b %in% c; a %in% e
UMCMCglmm<- updateable(MCMCglmm)
#The Priors
prior2 <- list(R = list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.02))); prior2.5 <- list(R = list(V=2, nu=0.002), G = list(G1 = list(V=1, nu=0.002), G2 = list(V=1, nu=0.002)))
prior3 <- list(R = list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.02), G3 = list(V=1, nu=0.02)))
prior2 <- list(R = list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.02)));
#without pedigree. Used for testing repeatability of outcome
m1NP <- MCMCglmm(RawScore ~ 1+ Age, random = ~ID + study, family = "gaussian",
               prior = prior2, pedigree = p2, verbose = FALSE, data = D, nitt = 150000, burnin = 15000, thin = 150)
summary(m1NP)
#The Models. Outcome = RawScore
#Final models: Null. Full with Rearing. Full with Rearing 2. Full with Rearing 3 (captive/wild). Study1. Study2.
#Try with Random = Publication
m1t <- MCMCglmm(RawScore ~ 1, random = ~animal + ID + Task, family = "gaussian",
               prior = prior3, pedigree = p2, verbose = FALSE, data = D, nitt = 150000, burnin = 15000, thin = 150)
m2t <- UMCMCglmm(RawScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID + Task, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250)
m3t <- UMCMCglmm(RawScore ~ 1 + sex + Rearing2 + Age + CurExp, random = ~animal + ID + Task, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250)
m4t <- UMCMCglmm(RawScore ~ 1 + sex + Rearing3 + Age + CurExp, random = ~animal + ID + Task, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250)
summary(m1t); summary(m2t); summary(m3t); summary(m4t)

m1 <- MCMCglmm(RawScore ~ 1, random = ~animal + ID + study2, family = "gaussian",
               prior = prior3, pedigree = p2, verbose = F, data = D, nitt = 500000, burnin = 100000, thin = 250)
m2 <- UMCMCglmm(RawScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID + study2, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250)
m3 <- UMCMCglmm(RawScore ~ 1 + sex + Rearing2 + Age + CurExp, random = ~animal + ID + study2, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250)
m4 <- UMCMCglmm(RawScore ~ 1 + sex + Rearing3 + Age + CurExp, random = ~animal + ID + study2, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 500000, verbose = F, burnin = 100000, thin = 250)

summary(m1); summary(m2); summary(m3); summary(m4)
DICdiff<-m5$DIC - m1$DIC


#Diagnostics
#Check for Autocorrelation
autocorr.diag(m1$Sol) ; autocorr.diag(m1$VCV)
autocorr.diag(m2$Sol) ; autocorr.diag(m2$VCV)
autocorr.diag(m3$Sol) ; autocorr.diag(m3$VCV)
autocorr.diag(m4$Sol) ; autocorr.diag(m4$VCV)
autocorr.diag(m5$Sol) ; autocorr.diag(m5$VCV)
#Check for mixing for each model
plot(m1$Sol)
plot(m1$VCV)
plot(m2$Sol)
plot(m2$VCV)
plot(m3$Sol)
plot(m3$VCV)
plot(m4$Sol)
plot(m4$VCV)
plot(m5$Sol)
plot(m5$VCV)
#Check effective sample sizes
effectiveSize(m1$Sol) ; effectiveSize(m1$VCV) ; effectiveSize(m2$Sol) ; effectiveSize(m2$VCV) ; effectiveSize(m3$Sol) ; effectiveSize(m3$VCV) ; effectiveSize(m4$Sol) ; effectiveSize(m4$VCV) ; effectiveSize(m5$Sol) ; effectiveSize(m5$VCV) 

#MEASURING REPEATABILITY
IC.1 <-(m4$VCV)/(rowSums(m4$VCV) + pi^2/3); plot(IC.1); summary(IC.1)

#Dredging for models. dm = dredge model. R = rearing. Hashed to prevent accidental running.
#dmR1<-dredge(m2, rank="DIC") 
#dmR2<-dredge(m3, rank="DIC") 
#dmR3<-dredge(m4, rank="DIC") 
dmR1; dmR2; dmR3
library(MuMIn)
dmR1s<- subset(dmR1, delta < 4)
dmR2s<- subset(dmR2, delta < 4)
dmR3s<- subset(dmR3, delta < 4)

dmR1s; dmR2s; dmR3s
mavg<-model.avg(dredgeModel);mavg
#withouth subset
dmR1avg<- model.avg(dmR1);dmR1avg
dmR2avg<- model.avg(dmR2);dmR2avg
dmR3avg<- model.avg(dmR3);dmR3avg
#with subset
dmR1avg<- model.avg(dmR1s);dmR1avg
dmR2avg<- model.avg(dmR2s);dmR2avg
dmR3avg<- model.avg(dmR3s);dmR3avg

summary(dmR1avg); summary(dmR2avg); summary(dmR3avg)
confint(dmR1avg, full = T); confint(dmR2avg, full = T); confint(dmR3avg, full = T);  

#phenotypic variance conditional on fixed effects

mVP<-m4$VCV[,"animal"]+m4$VCV[,"ID"]+m4$VCV[,"units"]
#between individual variance
mIP<-m4$VCV[,"animal"]+m4$VCV[,"ID"]
#repeatability
rep<-posterior.mode(mIP/mVP); rep
#heritability
her<-posterior.mode(m4$VCV[,"animal"]/mVP); her
her<-posterior.mode(m4$VCV[,"study"]/mVP); her

herit1 <- m1NP$VCV[, "ID"]/(m1NP$VCV[, "ID"] + m1NP$VCV[, "units"]); effectiveSize(herit1); mean(herit1); HPDinterval(herit1) 
herit1 <- m1NP$VCV[, "study"]/(m1NP$VCV[, "study"] + m1NP$VCV[, "units"]); effectiveSize(herit1); mean(herit1); HPDinterval(herit1) 

herit1 <- m5$VCV[, "ID"]/(m5$VCV[, "ID"] + m5$VCV[, "units"]); effectiveSize(herit1); mean(herit1); HPDinterval(herit1) 
herit1 <- m5$VCV[, "study"]/(m5$VCV[, "study"] + m5$VCV[, "units"]); effectiveSize(herit1); mean(herit1); HPDinterval(herit1)

#heritability
herit1 <- m1$VCV[, "animal"]/(m1$VCV[, "animal"] + m1$VCV[, "units"]); effectiveSize(herit1); mean(herit1); HPDinterval(herit1) 
herit2 <- m2$VCV[, "animal"]/(m2$VCV[, "animal"] + m2$VCV[, "units"]); effectiveSize(herit2); mean(herit2); HPDinterval(herit2) 
herit3 <- m3$VCV[, "animal"]/(m3$VCV[, "animal"] + m3$VCV[, "units"]); effectiveSize(herit3); mean(herit3); HPDinterval(herit3) 
herit4 <- m4$VCV[, "animal"]/(m4$VCV[, "animal"] + m4$VCV[, "units"]); effectiveSize(herit4); mean(herit4); HPDinterval(herit4)
herit5 <- m5$VCV[, "animal"]/(m5$VCV[, "animal"] + m5$VCV[, "units"]); effectiveSize(herit5); mean(herit5); HPDinterval(herit5)

#plotting coefficients
coefplot2(m2$VCV)
coefplot2(m2$Sol)
caterplot(m2$Sol)
xyplot(as.mcmc(m2$Sol))


#undone
sex * tool use?
