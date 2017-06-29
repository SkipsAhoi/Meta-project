# MCMCglmm animal model

#CHECK FOR APRIL
#CHECK FOR SAMMY
#CHECK FOR X'S

d <- as.matrix(read.table("~/Meta-project/MCMC_GLMM.txt", header=TRUE, na.strings="NA", sep = "\t"))
d <- as.data.frame(d)
d[is.na(d)] <- NA
str(d)

d$RawScore<-as.numeric(as.character(d$RawScore))
d$AdjScore<-as.numeric(as.character(d$AdjScore))
d$TotalExp<-as.numeric(as.character(d$TotalExp))
d$Age<-as.numeric(as.character(d$Age))
d$SLBinary<-as.numeric(as.character(d$SLBinary))

library(kinship2)
library(MCMCglmm)
library(pedigree)
library(pedigreemm)
library(MasterBayes)
#pedigree
p <- read.csv("kinship2.CSV", header = T, na.strings="NA")
p<- as.data.frame(p)
p[is.na(p)] <- NA
p
p2<-orderPed(p)
p2
p2[is.na(p2)] <- NA

str(p)
str(p2)
a<- p2$animal
b<- p2$dam
c<- p2$sire

e<- d$animal

b %in% a
c %in% a
e %in% a
e %in% a
b %in% c
a %in% e

difs <- setdiff(b,a)
difs <- setdiff(c,a)
difs <- setdiff(e,a)
length(difs)
list(difs)



#ID may have to be renamed 'animal'
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))

#prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002),
#            G2 = list(V=1, nu=0.002), G3=list(V=1, nu=0.002)))

#lets fit a simple model
m1 <- MCMCglmm(AdjScore ~ 1, random = ~animal, family = "gaussian",
                  prior = prior, pedigree = p2, data = d, nitt = 10000,
                  burnin = 1000, thin = 100)
summary(m1)

autocorr.diag(m1$Sol)
autocorr.diag(m1$VCV)

plot(m1$Sol)
plot(m1$VCV)
effectiveSize(m1$Sol)
effectiveSize(m1$VCV)

m2 <- MCMCglmm(AdjScore ~ sex + 1, random = ~animal, family = "gaussian",
                  prior = prior, pedigree = p2, data = d, nitt = 10000,
                  burnin = 1000, thin = 100)
summary(m2)

autocorr.diag(m2$Sol)
autocorr.diag(m2$VCV)

plot(m2$Sol)
plot(m2$VCV)
effectiveSize(m2$Sol)
effectiveSize(m2$VCV)

m3 <- MCMCglmm(AdjScore ~ sex + Rearing + 1, random = ~animal, family = "gaussian",
               prior = prior, pedigree = p2, data = d, nitt = 10000,
               burnin = 1000, thin = 100)
summary(m3)
autocorr.diag(m3$Sol)
autocorr.diag(m3$VCV)

plot(m3$Sol)
plot(m3$VCV)
effectiveSize(m3$Sol)
effectiveSize(m3$VCV)

m4 <- MCMCglmm(AdjScore ~ sex + Context + 1, random = ~animal, family = "gaussian",
               prior = prior, pedigree = p2, data = d, nitt = 10000,
               burnin = 1000, thin = 100)
summary(m4)

m4 <- MCMCglmm(AdjScore ~ sex + Context + Task + 1, random = ~animal, family = "gaussian",
               prior = prior, pedigree = p2, data = d, nitt = 10000,
               burnin = 1000, thin = 100)
summary(m4)

autocorr.diag(m3$Sol)
autocorr.diag(m23VCV)

plot(m3$Sol)
plot(m3$VCV)
effectiveSize(m3$Sol)
effectiveSize(m3$VCV)



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

