# Packages required
library(kinship2) ; library(MCMCglmm) ; library(pedigree) ; library(pedigreemm) ; library(MasterBayes) ; library(rethinking) ; library(MuMIn) ; library(AICcmodavg) ; library(sjstats); library(coda)

#datafile setup and variable correction
d <- as.matrix(read.table("~/Meta-project/MCMC_GLMM.txt", header=TRUE, na.strings="NA", sep = "\t")); d <- as.data.frame(d); d[is.na(d)] <- NA;str(d)
d$RawScore<-as.numeric(as.character(d$RawScore)); d$AdjScore<-as.numeric(as.character(d$AdjScore)); d$TotalExp<-as.numeric(as.character(d$TotalExp))
d$Age<-as.numeric(as.character(d$Age)); d$SLBinary<-as.numeric(as.character(d$SLBinary)); d$CurExp<-as.numeric(as.character(d$CurExp)); str(d)
d$ID<- d$animal
#eliminate NA lines to streamline the dataset
d<-subset(d, !is.na(RawScore))
str(d)

# Let's subset out all of the ghost conditions
D<-d
D<-subset(D,!(Context == "Ghost")) ; D$Context<-drop.levels(D$Context) ; str(D)
# And then let's remove unknown rearing histories
D<-subset(D,!(Rearing == "U")) ; D$Rearing<-drop.levels(D$Rearing) ; str(D)
D<-subset(D,!(Age.group == "1")) ; D$Age.group<-drop.levels(D$Age.group) ; str(D)
#subsetting individuals
#D<-subset(D, TotalExp > 1) ; D$TotalExp<-drop.levels(D$TotalExp) ; str(D)


#This is our pedigree file. One column for animals, one for sires, one for dams. All sires and dams must be in animal column.
p <- read.csv("~/Meta-project/kinship2.CSV", header = T, na.strings="NA"); p<- as.data.frame(p); p[is.na(p)] <- NA; str(p)
#This command orders the pedigree file so that sires+dams come first in the list. Necessary for MCMCglmm to read it properly.
p2<-orderPed(p);p2[is.na(p2)] <- NA; str(p2)

# Let's check to make sure all Sires and Dams are in the animal list
a <- p2$animal; b <- p2$dam; c <- p2$sire; e <- d$animal
difs <- setdiff(b,a); length(difs); list(difs); difs <- setdiff(c,a);length(difs); list(difs); difs <- setdiff(e,a);length(difs); list(difs)
#for manual checking use:
#b %in% a; c %in% a; e %in% a; e %in% a; b %in% c; a %in% e

# Priors: We're going to use weakly informative priors here
#For one random effect
prior1 <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002))); prior1.5 <- list(R = list(V=2, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))
#for two random effects use:
prior2 <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002), G2 = list(V=1, nu=0.002))); prior2.5 <- list(R = list(V=2, nu=0.002), G = list(G1 = list(V=1, nu=0.002), G2 = list(V=1, nu=0.002)))

#Null model
m1 <- MCMCglmm(AdjScore ~ 1, random = ~animal + study, family = "gaussian",
               prior = prior2, pedigree = p2, verbose = FALSE, data = d, nitt = 100000, burnin = 10000, thin = 100, pr = T)
summary(m1)

autocorr.diag(m1$Sol) ; autocorr.diag(m1$VCV)
plot(m1$Sol)
plot(m1$VCV)
effectiveSize(m1$Sol) ; effectiveSize(m1$VCV)

m2 <- MCMCglmm(AdjScore ~ sex + 1, random = ~animal + study, family = "gaussian",
               prior = prior2, pedigree = p2, data = d, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100, pr = T)
summary(m2)

autocorr.diag(m2$Sol);autocorr.diag(m2$VCV)
plot(m2$Sol)
plot(m2$VCV)
effectiveSize(m2$Sol) ; effectiveSize(m2$VCV)

m3 <- MCMCglmm(AdjScore ~ sex + Rearing + 1, random = ~animal + study, family = "gaussian",
               prior = prior2, pedigree = p2, data = d, nitt = 100000, burnin = 10000, thin = 100, pr = T)
summary(m3)
autocorr.diag(m3$Sol) ; autocorr.diag(m3$VCV)
plot(m3$Sol)
plot(m3$VCV)
effectiveSize(m3$Sol) ; effectiveSize(m3$VCV)

m4 <- MCMCglmm(AdjScore ~ sex + Rearing + Age + 1, random = ~animal + study, family = "gaussian",
               prior = prior2, pedigree = p2, data = d, nitt = 100000, burnin = 10000, thin = 100, pr = T)
summary(m4)

autocorr.diag(m4$Sol) ; autocorr.diag(m4$VCV)
plot(m4$Sol)
plot(m4$VCV)
effectiveSize(m4$Sol) ; effectiveSize(m4$VCV)

m5 <- MCMCglmm(RawScore ~ sex + Rearing + Age + CurExp, random = ~animal, family = "gaussian",
               prior = prior1, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
summary(m5)

m5.5 <- MCMCglmm(RawScore ~ sex + Rearing + Age, random = ~animal, family = "gaussian",
                prior = prior1, pedigree = p2, data = d, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
summary(m5.5)

IC.1 <-m5.5$VCV[, 1]/(rowSums(m5.5$VCV) + pi^2/3); plot(IC.1); summary(IC.1)

plot(IC.1)
autocorr.diag(m5$Sol) ; autocorr.diag(m5$VCV)
plot(m5$Sol)
plot(m5$VCV)
effectiveSize(m5$Sol) ; effectiveSize(m5$VCV)

m6 <- MCMCglmm(SLBinary ~ sex + Rearing + Age + CurExp, random = ~animal, family = "categorial",
                prior = prior1, pedigree = p2, data = d, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
summary(m6)

IC.2 <-m6$VCV[, 1]/(rowSums(m6$VCV) + pi^2/3); summary(IC.2); plot(IC.2)
summary(IC.2)
IC.2
plot(IC.2)

c2 <- ((16 * sqrt(3))/(15 * pi))^2
Int.1 <- m6$Sol/sqrt(1 + c2 * m6$VCV[, 2]); summary (Int.1); plot(Int.1)

autocorr.diag(m6$Sol) ; autocorr.diag(m6$VCV)
plot(m6$Sol)
plot(m6$VCV)
effectiveSize(m6$Sol) ; effectiveSize(m6$VCV)

m7 <- glmer(SLBinary ~ sex + Rearing + Age +  (1|animal) + (1|CurExp), data = d, family = binomial)
summary(m7)
print(m7)
icc(m7)

m8 <-   rptGaussian(RawScore ~ sex + Rearing + Age +  (1|animal) + (1|), data = d, grname = "animal", nboot= 10, npermut= 0, adjusted=F)
summary(m8)
print(m8)
icc(m7)

dred <- dredge(m5)
avg <-model.avg(dred)
summary(avg)

m6 <- uMCMCglmm(SLBinary ~ sex + Rearing + Age + CurExp, random = ~animal + study, family = "gaussian",
                prior = prior2, pedigree = p2, data = d, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
summary(m6)

#heritability
herit <- m1$VCV[, "animal"]/(m1$VCV[, "animal"] + m1$VCV[, "units"])
effectiveSize(herit)
mean(herit)
HPDinterval(herit) 

herit <- m2$VCV[, "animal"]/(m2$VCV[, "animal"] + m2$VCV[, "units"])
effectiveSize(herit)
mean(herit)
HPDinterval(herit) 

herit <- m3$VCV[, "animal"]/(m3$VCV[, "animal"] + m3$VCV[, "units"])
effectiveSize(herit)
mean(herit)
HPDinterval(herit) 

herit <- m4$VCV[, "animal"]/(m4$VCV[, "animal"] + m4$VCV[, "units"])
effectiveSize(herit)
mean(herit)
HPDinterval(herit)

herit <- m5$VCV[, "animal"]/(m5$VCV[, "animal"] + m5$VCV[, "units"])
effectiveSize(herit)
mean(herit)
HPDinterval(herit)



