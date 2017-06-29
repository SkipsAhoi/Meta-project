#FINAL meta-anlaysis R scripts

#R Packages - some unused.

library(kinship2) ; library(MCMCglmm) ; library(pedigree) ; library(pedigreemm) ; library(MasterBayes) ; library(MuMIn) ; library(AICcmodavg) ; library(sjstats); library(coda)
library(lme4);library(arm);library(sjPlot);library(glmulti);library(bblme);library(Rcpp);library(coefplot);library(ggplot2); library(multcomp); library(blmeco)
library(gridExtra);library(cowplot);library(ggplot2); library(lattice); library(coefplot2); library(mcmcplots); library(coda); library(mfx); library(spatialprobit); library(MASS); library(boot)

#The Data
#datafile setup and variable correction
d <- as.matrix(read.table("~/Meta-project/MCMC_GLMM.txt", header=TRUE, na.strings="NA", sep = "\t")); d <- as.data.frame(d); d[is.na(d)] <- NA;str(d)
d$AdjScore<-as.numeric(as.character(d$AdjScore)); d$AdjScore<-as.numeric(as.character(d$AdjScore)); d$TotalExp<-as.numeric(as.character(d$TotalExp))
d$Age<-as.numeric(as.character(d$Age)); d$SLBinary<-as.numeric(as.character(d$SLBinary)); d$CurExp<-as.numeric(as.character(d$CurExp)); str(d)
d$ID<- d$animal
#eliminate NA rows to streamline the dataset
d<-subset(d, !is.na(AdjScore)); str(d)
# Let's subset out all of the ghost conditions, which don't have social learning
D<-d; D<-subset(D,!(Context == "Ghost")) ; D$Context<-drop.levels(D$Context) ; str(D)
# And then let's remove individuals with unknown rearing histories
D<-subset(D,!(Rearing == "U")) ; D$Rearing<-drop.levels(D$Rearing) ; str(D)
#drop individuals of age 1 
D<-subset(D,!(Age.group == "1")) ; D$Age.group<-drop.levels(D$Age.group) ; str(D)
#subsetting individuals according to participation. Not sure about this.
#D<-subset(D, TotalExp > 1) ; D$TotalExp<-drop.levels(D$TotalExp) ; str(D)
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

#The Priors
prior3 <- list(R = list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.02), G3 = list(V=1, nu=0.02)))
prior2 <- list(R = list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.02)))

#The Models. Outcome = AdjScore. No ranom=study because each outcome is already centered on study means.
mA1 <- MCMCglmm(AdjScore ~ 1, random = ~animal + ID, family = "gaussian",
               prior = prior2, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)
mA2 <- MCMCglmm(AdjScore ~ 1 + sex, random = ~animal + ID, family = "gaussian",
               prior = prior2, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)
mA3 <- MCMCglmm(AdjScore ~ 1 + sex + Rearing, random = ~animal + ID, family = "gaussian",
               prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA4 <- MCMCglmm(AdjScore ~ 1 + sex + Rearing + Age, random = ~animal + ID, family = "gaussian",
               prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA5 <- MCMCglmm(AdjScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID, family = "gaussian",
               prior = prior2, pedigree = p2, data = D, nitt = 150000, verbose = F, burnin = 15000, thin = 150)
summary(mA1); summary(mA2); summary(mA3); summary(mA4); summary(mA5)
summary(mA5)
#Diagnostics
#Check for Autocorrelation
autocorr.diag(mA1$Sol) ; autocorr.diag(mA1$VCV)
autocorr.diag(mA2$Sol) ; autocorr.diag(mA2$VCV)
autocorr.diag(mA3$Sol) ; autocorr.diag(mA3$VCV)
autocorr.diag(mA4$Sol) ; autocorr.diag(mA4$VCV)
autocorr.diag(mA5$Sol) ; autocorr.diag(mA5$VCV)
#Check for mixing for each model
plot(mA1$Sol)
plot(mA1$VCV)
plot(mA2$Sol)
plot(mA2$VCV)
plot(mA3$Sol)
plot(mA3$VCV)
plot(mA4$Sol)
plot(mA4$VCV)
plot(mA5$Sol)
plot(mA5$VCV)
#Check effective sample sizes
effectiveSize(mA1$Sol) ; effectiveSize(mA1$VCV) ; effectiveSize(mA2$Sol) ; effectiveSize(mA2$VCV) ; effectiveSize(mA3$Sol) ; effectiveSize(mA3$VCV) ; effectiveSize(mA4$Sol) ; effectiveSize(mA4$VCV) ; effectiveSize(mA5$Sol) ; effectiveSize(mA5$VCV) 

#MEASURING REPEATABILITY
IC.1 <-mA5$VCV[, 2]/(rowSums(mA5$VCV) + pi^2/3); plot(IC.1); summary(IC.1)

#Inference
#Model Comparison
ms<-model.sel(mA1, mA2, mA3, mA4, mA5, rank = DIC); ms
ms2<-model.sel(mA1, mA2, mA3, mA4, mA5, m6, m7, m8, m9, mA10, rank = DIC)
ms; ms2

summary(mA3)
#phenotypic variance conditional on fixed effects
mVP<-mA4$VCV[,"animal"]+mA4$VCV[,"ID"]+mA4$VCV[,"units"]
#between individual variance
mIP<-mA4$VCV[,"animal"]+mA4$VCV[,"ID"]
#repeatability
rep<-posterior.mode(mAIP/mVP); rep
#heritability
her<-posterior.mode(mA5$VCV[,"animal"]/mVP); her

#heritability
herit1 <- mA1$VCV[, "animal"]/(mA1$VCV[, "animal"] + mA1$VCV[, "units"]); effectiveSize(herit1); mean(herit1); HPDinterval(herit1) 
herit2 <- mA2$VCV[, "animal"]/(mA2$VCV[, "animal"] + mA2$VCV[, "units"]); effectiveSize(herit2); mean(herit2); HPDinterval(herit2) 
herit3 <- mA3$VCV[, "animal"]/(mA3$VCV[, "animal"] + mA3$VCV[, "units"]); effectiveSize(herit3); mean(herit3); HPDinterval(herit3) 
herit4 <- mA4$VCV[, "animal"]/(mA4$VCV[, "animal"] + mA4$VCV[, "units"]); effectiveSize(herit4); mean(herit4); HPDinterval(herit4)
herit5 <- mA5$VCV[, "animal"]/(mA5$VCV[, "animal"] + mA5$VCV[, "units"]); effectiveSize(herit5); mean(herit5); HPDinterval(herit5)

#plotting coefficients
coefplot2(mA5$VCV)
coefplot2(mA5$Sol)
caterplot#FINAL meta-anlaysis R scripts

#R Packages - some unused.

library(kinship2) ; library(MCMCglmm) ; library(pedigree) ; library(pedigreemm) ; library(MasterBayes) ; library(MuMIn) ; library(AICcmodavg) ; library(sjstats); library(coda)
library(lme4);library(arm);library(sjPlot);library(glmulti);library(bblme);library(Rcpp);library(coefplot);library(ggplot2); library(multcomp); library(blmeco)
library(gridExtra);library(cowplot);library(ggplot2); library(lattice); library(coefplot2); library(mcmcplots); library(coda); library(mfx); library(spatialprobit); library(MASS); library(boot)

#The Data
#datafile setup and variable correction
d <- as.matrix(read.table("~/Meta-project/MCMC_GLMM.txt", header=TRUE, na.strings="NA", sep = "\t")); d <- as.data.frame(d); d[is.na(d)] <- NA;str(d)
d$RawScore<-as.numeric(as.character(d$RawScore)); d$AdjScore<-as.numeric(as.character(d$AdjScore)); d$TotalExp<-as.numeric(as.character(d$TotalExp))
d$Age<-as.numeric(as.character(d$Age)); d$SLBinary<-as.numeric(as.character(d$SLBinary)); d$CurExp<-as.numeric(as.character(d$CurExp)); str(d)
d$ID<- d$animal
#eliminate NA rows to streamline the dataset
d<-subset(d, !is.na(RawScore)); str(d)
# Let's subset out all of the ghost conditions, which don't have social learning
D<-d; D<-subset(D,!(Context == "Ghost")) ; D$Context<-drop.levels(D$Context) ; str(D)
# And then let's remove individuals with unknown rearing histories
D<-subset(D,!(Rearing == "U")) ; D$Rearing<-drop.levels(D$Rearing) ; str(D)
#drop individuals of age 1 
D<-subset(D,!(Age.group == "1")) ; D$Age.group<-drop.levels(D$Age.group) ; str(D)
#subsetting individuals according to participation. Not sure about this.
#D<-subset(D, TotalExp > 1) ; D$TotalExp<-drop.levels(D$TotalExp) ; str(D)
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

#The Priors
prior3 <- list(R = list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.02), G3 = list(V=1, nu=0.02)))

#The Models. Outcome = RawScore
mA1 <- MCMCglmm(RawScore ~ 1, random = ~animal + ID + Task, family = "gaussian",
               prior = prior3, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)
mA2 <- MCMCglmm(RawScore ~ 1 + sex, random = ~animal + ID + Task, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)
mA3 <- MCMCglmm(RawScore ~ 1 + sex + Rearing, random = ~animal + ID + Task, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA4 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age, random = ~animal + ID + Task, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA5 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID+ Task, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
#Diagnostics
#Check for Autocorrelation
autocorr.diag(mA1$Sol) ; autocorr.diag(mA1$VCV)
autocorr.diag(mA2$Sol) ; autocorr.diag(mA2$VCV)
autocorr.diag(mA3$Sol) ; autocorr.diag(mA3$VCV)
autocorr.diag(mA4$Sol) ; autocorr.diag(mA4$VCV)
autocorr.diag(mA5$Sol) ; autocorr.diag(mA5$VCV)
#Check for mixing for each model
plot(mA1$Sol)
plot(mA1$VCV)
plot(mA2$Sol)
plot(mA2$VCV)
plot(mA3$Sol)
plot(mA3$VCV)
plot(mA4$Sol)
plot(mA4$VCV)
plot(mA5$Sol)
plot(mA5$VCV)
#Check effective sample sizes
effectiveSize(mA1$Sol) ; effectiveSize(mA1$VCV) ; effectiveSize(mA2$Sol) ; effectiveSize(mA2$VCV) ; effectiveSize(mA3$Sol) ; effectiveSize(mA3$VCV) ; effectiveSize(mA4$Sol) ; effectiveSize(mA4$VCV) ; effectiveSize(mA5$Sol) ; effectiveSize(mA5$VCV) 

#MEASURING REPEATABILITY
IC.1 <-mA5$VCV[, 2]/(rowSums(mA5$VCV) + pi^2/3); plot(IC.1); summary(IC.1)

#Inference
#Model Comparison
ms<-model.sel(mA1, mA2, mA3, mA4, mA5, rank = DIC)
ms2<-model.sel(mA1, mA2, mA3, mA4, mA5, m6, m7, m8, m9, mA10, rank = DIC)
ms; ms2

summary(m8)
#phenotypic variance conditional on fixed effects
mVP<-mA4$VCV[,"animal"]+mA4$VCV[,"ID"]+mA4$VCV[,"units"]
#between individual variance
mIP<-mA4$VCV[,"animal"]+mA4$VCV[,"ID"]
#repeatability
rep<-posterior.mode(mIP/mVP); rep
#heritability
her<-posterior.mode(mA5$VCV[,"animal"]/mVP); her

#heritability
herit1 <- mA1$VCV[, "animal"]/(mA1$VCV[, "animal"] + mA1$VCV[, "units"]); effectiveSize(herit1); mean(herit1); HPDinterval(herit1) 
herit2 <- mA2$VCV[, "animal"]/(mA2$VCV[, "animal"] + mA2$VCV[, "units"]); effectiveSize(herit2); mean(herit2); HPDinterval(herit2) 
herit3 <- mA3$VCV[, "animal"]/(mA3$VCV[, "animal"] + mA3$VCV[, "units"]); effectiveSize(herit3); mean(herit3); HPDinterval(herit3) 
herit4 <- mA4$VCV[, "animal"]/(mA4$VCV[, "animal"] + mA4$VCV[, "units"]); effectiveSize(herit4); mean(herit4); HPDinterval(herit4)
herit5 <- mA5$VCV[, "animal"]/(mA5$VCV[, "animal"] + mA5$VCV[, "units"]); effectiveSize(herit5); mean(herit5); HPDinterval(herit5)

#plotting coefficients
coefplot2(mA5$VCV)
coefplot2(mA5$Sol)
caterplot(mA5$Sol)
xyplot(as.mcmc(mA5$Sol))



#EXTRAS



#with Context instead of Task
m6 <- MCMCglmm(RawScore ~ 1, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)
m7 <- MCMCglmm(RawScore ~ 1 + sex, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)
m8 <- MCMCglmm(RawScore ~ 1 + sex + Rearing, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m9 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA10 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID+ Context, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

#with Study instead of Task
m6 <- MCMCglmm(RawScore ~ 1, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)
m7 <- MCMCglmm(RawScore ~ 1 + sex, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)
m8 <- MCMCglmm(RawScore ~ 1 + sex + Rearing, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m9 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA10 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID+ study, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)#FINAL meta-anlaysis R scripts

#R Packages - some unused.

library(kinship2) ; library(MCMCglmm) ; library(pedigree) ; library(pedigreemm) ; library(MasterBayes) ; library(MuMIn) ; library(AICcmodavg) ; library(sjstats); library(coda)
library(lme4);library(arm);library(sjPlot);library(glmulti);library(bblme);library(Rcpp);library(coefplot);library(ggplot2); library(multcomp); library(blmeco)
library(gridExtra);library(cowplot);library(ggplot2); library(lattice); library(coefplot2); library(mcmcplots); library(coda); library(mfx); library(spatialprobit); library(MASS); library(boot)

#The Data
#datafile setup and variable correction
d <- as.matrix(read.table("~/Meta-project/MCMC_GLMM.txt", header=TRUE, na.strings="NA", sep = "\t")); d <- as.data.frame(d); d[is.na(d)] <- NA;str(d)
d$RawScore<-as.numeric(as.character(d$RawScore)); d$AdjScore<-as.numeric(as.character(d$AdjScore)); d$TotalExp<-as.numeric(as.character(d$TotalExp))
d$Age<-as.numeric(as.character(d$Age)); d$SLBinary<-as.numeric(as.character(d$SLBinary)); d$CurExp<-as.numeric(as.character(d$CurExp)); str(d)
d$ID<- d$animal
#eliminate NA rows to streamline the dataset
d<-subset(d, !is.na(RawScore)); str(d)
# Let's subset out all of the ghost conditions, which don't have social learning
D<-d; D<-subset(D,!(Context == "Ghost")) ; D$Context<-drop.levels(D$Context) ; str(D)
# And then let's remove individuals with unknown rearing histories
D<-subset(D,!(Rearing == "U")) ; D$Rearing<-drop.levels(D$Rearing) ; str(D)
#drop individuals of age 1 
D<-subset(D,!(Age.group == "1")) ; D$Age.group<-drop.levels(D$Age.group) ; str(D)
#subsetting individuals according to participation. Not sure about this.
#D<-subset(D, TotalExp > 1) ; D$TotalExp<-drop.levels(D$TotalExp) ; str(D)
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

#The Priors
prior3 <- list(R = list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.02), G3 = list(V=1, nu=0.02)))

#The Models. Outcome = RawScore
mA1 <- MCMCglmm(RawScore ~ 1, random = ~animal + ID + Task, family = "gaussian",
               prior = prior3, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)
mA2 <- MCMCglmm(RawScore ~ 1 + sex, random = ~animal + ID + Task, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)
mA3 <- MCMCglmm(RawScore ~ 1 + sex + Rearing, random = ~animal + ID + Task, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA4 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age, random = ~animal + ID + Task, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA5 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID+ Task, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
#Diagnostics
#Check for Autocorrelation
autocorr.diag(mA1$Sol) ; autocorr.diag(mA1$VCV)
autocorr.diag(mA2$Sol) ; autocorr.diag(mA2$VCV)
autocorr.diag(mA3$Sol) ; autocorr.diag(mA3$VCV)
autocorr.diag(mA4$Sol) ; autocorr.diag(mA4$VCV)
autocorr.diag(mA5$Sol) ; autocorr.diag(mA5$VCV)
#Check for mixing for each model
plot(mA1$Sol)
plot(mA1$VCV)
plot(mA2$Sol)
plot(mA2$VCV)
plot(mA3$Sol)
plot(mA3$VCV)
plot(mA4$Sol)
plot(mA4$VCV)
plot(mA5$Sol)
plot(mA5$VCV)
#Check effective sample sizes
effectiveSize(mA1$Sol) ; effectiveSize(mA1$VCV) ; effectiveSize(mA2$Sol) ; effectiveSize(mA2$VCV) ; effectiveSize(mA3$Sol) ; effectiveSize(mA3$VCV) ; effectiveSize(mA4$Sol) ; effectiveSize(mA4$VCV) ; effectiveSize(mA5$Sol) ; effectiveSize(mA5$VCV) 

#MEASURING REPEATABILITY
IC.1 <-mA5$VCV[, 2]/(rowSums(mA5$VCV) + pi^2/3); plot(IC.1); summary(IC.1)

#Inference
#Model Comparison
ms<-model.sel(mA1, mA2, mA3, mA4, mA5, rank = DIC)
ms2<-model.sel(mA1, mA2, mA3, mA4, mA5, m6, m7, m8, m9, mA10, rank = DIC)
ms; ms2

summary(m8)
#phenotypic variance conditional on fixed effects
mVP<-mA4$VCV[,"animal"]+mA4$VCV[,"ID"]+mA4$VCV[,"units"]
#between individual variance
mIP<-mA4$VCV[,"animal"]+mA4$VCV[,"ID"]
#repeatability
rep<-posterior.mode(mIP/mVP); rep
#heritability
her<-posterior.mode(mA5$VCV[,"animal"]/mVP); her

#heritability
herit1 <- mA1$VCV[, "animal"]/(mA1$VCV[, "animal"] + mA1$VCV[, "units"]); effectiveSize(herit1); mean(herit1); HPDinterval(herit1) 
herit2 <- mA2$VCV[, "animal"]/(mA2$VCV[, "animal"] + mA2$VCV[, "units"]); effectiveSize(herit2); mean(herit2); HPDinterval(herit2) 
herit3 <- mA3$VCV[, "animal"]/(mA3$VCV[, "animal"] + mA3$VCV[, "units"]); effectiveSize(herit3); mean(herit3); HPDinterval(herit3) 
herit4 <- mA4$VCV[, "animal"]/(mA4$VCV[, "animal"] + mA4$VCV[, "units"]); effectiveSize(herit4); mean(herit4); HPDinterval(herit4)
herit5 <- mA5$VCV[, "animal"]/(mA5$VCV[, "animal"] + mA5$VCV[, "units"]); effectiveSize(herit5); mean(herit5); HPDinterval(herit5)

#plotting coefficients
coefplot2(mA5$VCV)
coefplot2(mA5$Sol)
caterplot(mA5$Sol)
xyplot(as.mcmc(mA5$Sol))



#EXTRAS



#with Context instead of Task
m6 <- MCMCglmm(RawScore ~ 1, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)
m7 <- MCMCglmm(RawScore ~ 1 + sex, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)
m8 <- MCMCglmm(RawScore ~ 1 + sex + Rearing, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m9 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA10 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID+ Context, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

#with Study instead of Task
m6 <- MCMCglmm(RawScore ~ 1, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)
m7 <- MCMCglmm(RawScore ~ 1 + sex, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)
m8 <- MCMCglmm(RawScore ~ 1 + sex + Rearing, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m9 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA10 <- MCMCglmm(RawScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID+ study, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)5$Sol)
xyplot(as.mcmc(mA5$Sol))



#EXTRAS



#with Context instead of study
m6 <- MCMCglmm(AdjScore ~ 1, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)
m7 <- MCMCglmm(AdjScore ~ 1 + sex, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)
m8 <- MCMCglmm(AdjScore ~ 1 + sex + Rearing, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m9 <- MCMCglmm(AdjScore ~ 1 + sex + Rearing + Age, random = ~animal + ID + Context, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA10 <- MCMCglmm(AdjScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID+ Context, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

#with Study instead of study
m6 <- MCMCglmm(AdjScore ~ 1, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)
m7 <- MCMCglmm(AdjScore ~ 1 + sex, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)
m8 <- MCMCglmm(AdjScore ~ 1 + sex + Rearing, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m9 <- MCMCglmm(AdjScore ~ 1 + sex + Rearing + Age, random = ~animal + ID + study, family = "gaussian",
               prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
mA10 <- MCMCglmm(AdjScore ~ 1 + sex + Rearing + Age + CurExp, random = ~animal + ID+ study, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)