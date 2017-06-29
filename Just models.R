library(gridExtra);library(cowplot);library(ggplot2); library(lattice); library
#install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R",type="source")
library(coefplot2); library(mcmcplots); library(coda); library(mfx); library(spatialprobit)
library(rethinking); library(MuMIn)
library(kinship2) ; library(MCMCglmm) ; library(pedigree) ; library(pedigreemm) ; library(MasterBayes) ; library(rethinking) ; library(MuMIn) ; library(AICcmodavg) ; library(sjstats); library(coda)
library(MuMIn);
library(MASS); library(boot)

prior1 <- list(R = list(V=1, nu = 0.02), G = list(G1 = list(V=1, nu = 0.02))); prior1.5 <- list(R = list(V=2, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))
#for two random effects use:
prior2 <- list(R = list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.002))); prior2.5 <- list(R = list(V=2, nu=0.002), G = list(G1 = list(V=1, nu=0.002), G2 = list(V=1, nu=0.002)))
prior3 <- list(R = list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.02), G3 = list(V=1, nu=0.02)))


#phenotypic variance conditional on fixed effects
mVP<-m1x$VCV[,"animal"]+m1x$VCV[,"ID"]+m1x$VCV[,"units"]
#between individual variance
mIP<-m1x$VCV[,"animal"]+m1x$VCV[,"ID"]
#repeatability
posterior.mode(mIP/mVP)

m1x <- MCMCglmm(RawScore ~ 1, random = ~animal + ID , family = "gaussian",
                prior = prior2, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)

m2x <- MCMCglmm(RawScore ~ sex, random = ~animal + ID , family = "gaussian",
               prior = prior2, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)

m3x <- MCMCglmm(RawScore ~ sex + Rearing, random = ~animal + ID , family = "gaussian",
               prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

m4x <- MCMCglmm(RawScore ~ sex + Rearing + Age, random = ~animal + ID , family = "gaussian",
               prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

m5x <- MCMCglmm(RawScore ~ sex + Rearing + Age + CurExp, random = ~animal + ID, family = "gaussian",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

m6x <- MCMCglmm(RawScore ~ sex + Age, random = ~animal + ID, family = "gaussian",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m7x <- MCMCglmm(RawScore ~ sex + Age  +CurExp, random = ~animal + ID, family = "gaussian",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

m8x <- MCMCglmm(RawScore ~ Rearing, random = ~animal + ID, family = "gaussian",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

m9x <- MCMCglmm(RawScore ~ Rearing + Age, random = ~animal + ID, family = "gaussian",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

m9x <- MCMCglmm(RawScore ~ Rearing + Age + CurExp, random = ~animal + ID, family = "gaussian",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

m9x <- MCMCglmm(RawScore ~ Rearing + CurExp, random = ~animal + ID, family = "gaussian",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

m10x <- MCMCglmm(RawScore ~ sex+CurExp, random = ~animal + ID, family = "gaussian",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

m11x <- MCMCglmm(RawScore ~ Age+CurExp, random = ~animal + ID, family = "gaussian",
                 prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

m12x <- MCMCglmm(RawScore ~ Age, random = ~animal + ID, family = "gaussian",
                 prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)

(m5x)

model.avg(m1x, m2x, m3x, m4x, m5x, m6x, m7x, m8x, m9x, m10x, m11x, m12x, subset = delta < 4)
ms<-model.sel(m1x, m2x, m3x, m4x, m5x, m6x, m7x, m8x, m9x, m10x, m11x, m12x, rank = DIC)
ms
mavg<-model.avg(ms,  subset = delta < 4)
mavg
confint(mavg)
mavg$coef.shrinkage

summary(m5x)

plot(plogis(.55-0.10*(0:1000)))
inv.pr
plogis(0.55-0.1)
fixef(m5x)

pframe <- data.frame(sex=factor(levels(D$sex),
                                levels=levels(D$sex)))
pframe
cpred1 <- predict(m5x,re.form=NA,newdata=pframe,type="response")
predict(m5x)

W.1<-cBind(m5x$X, m5x$Z) # note X and Z are sparse so use cBind
prediction.1<-W.1%*%posterior.mode(m5x$Sol)
exp(m5x$coef["(Intercept)"])

exp(mean(m5x$Sol[, "(Intercept)"])) + 0.5 * m5x$VCV[,+ 1]
m<-exp(mean(m5x$Sol[, "sexM"])) + 0.5 * m5x$VCV[,+ 1]
m<-inv.logit(mean(m5x$Sol[, "sexM"])) + 0.5 * m5x$VCV[,+ 1]
mean(m)
inv.logit(0.8)


Int <- t(apply(m5x$Sol[, 1:2], 1, function(x) {
  + D %*% (x/sqrt(1 + c2 * diag(IJ)))
  + }))
summary(mcmc(exp(Int)/rowSums(exp(Int))))


exp(0.82)
exp(0.82-0.105)

mfx(m5x)

W.1<-cBind(m6x$X, m6x$Z) # note X and Z are sparse so use cBind
prediction.1<-W.1%*%posterior.mode(m6x$Sol)
plot1<-xyplot(RawScore+prediction.1@x~Age|ID, data=D)
plot1

xyplot(RawScore+ predict(m5x, marginal = NULL) ~
          sex | ID, data = D)


pop.int <- posterior.mode(m5x$Sol[, 1])
pop.slope <- posterior.mode(m5x$Sol[, 5])
pop.quad <- posterior.mode(m5x$Sol[, 6])
chick.int <- posterior.mode(m5x$Sol[, c(7:56)])


#DIC list
model_list<-list(m1,m2,m3,m4,m5)
sapply(model_list,"[[","DIC")
coefplot2(m5$VCV)
coefplot2(m5$Sol)
caterplot(m5$Sol)
xyplot(as.mcmc(m5x$Sol))

plot(D$RawScore, predict(m5x, marginal=m5x$ID))

#plotting predictions

p1 <- ggplot(dat = D, aes(x = age, y = RawScore)) +
  geom_smooth(aes(x = age, ymin = lower, ymax = upper), stat = "identity") + facet_wrap(~ Interaction)

x1<-predict(m5x, newdata=NULL, marginal=m5$ID,
        type="response", interval="none", level=0.95)
plot(x1)

summary(m5x)
plot(m5x)

plot(m5x$VCV)
HPDinterval(m5x$Sol)

pred1<- predict(m5x, newdata=NULL, marginal= NULL,
        type="response", interval="confidence", level=0.95)
plot(pred1)
colnames(m5x$Sol)

#phenotypic variance conditional on fixed effects
mVP<-m5x$VCV[,"animal"]+m5x$VCV[,"ID"]+m5x$VCV[,"units"]
#between individual variance
mIP<-m5x$VCV[,"animal"]+m5x$VCV[,"ID"]
#repeatability
posterior.mode(mIP/mVP)
#heritability
posterior.mode(m5x$VCV[,"animal"]/mVP)

gelman.diag(mcmc.list(m4$Sol, m4$Sol))
gelman.plot(mcmc.list(m5$Sol, m5$Sol))

m1c <- MCMCglmm(SLBinary ~ 1, random = ~animal + ID, family = "categorical",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, burnin = 10000, thin = 100)
m2c <- MCMCglmm(SLBinary ~ sex, random = ~animal+ID, family = "categorical",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m3c <- MCMCglmm(SLBinary ~ sex + Rearing, random = ~animal+ID, family = "categorical",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m4c <- MCMCglmm(SLBinary ~ sex + Rearing + Age, random = ~animal + ID, family = "categorical",
               prior = prior2, pedigree = p2, data = D, nitt = 100000, burnin = 10000, thin = 100)
m5c <- MCMCglmm(SLBinary ~ sex + Rearing + Age + CurExp, random = ~animal+ID, family = "ordinal",
               prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
summary(m5c)
inv.logit(2.95-0.44)

model_list<-list(m1c,m2c,m3c,m4c,m5c)
sapply(model_list,"[[","DIC")

m1x <- MCMCglmm(RawScore ~ 1, random = ~animal + ID + Task, family = "gaussian",
                prior = prior3, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)

m2x <- MCMCglmm(RawScore ~ 1 + sex, random = ~animal + ID + Task, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, verbose = FALSE, nitt = 100000, burnin = 10000, thin = 100)
m3x <- MCMCglmm(RawScore ~ 1 + sex + Rearing2, random = ~animal + ID + Task, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m4x <- MCMCglmm(RawScore ~ 1 + sex + Rearing2 + Age, random = ~animal + ID + Task, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
m5x <- MCMCglmm(RawScore ~ 1 + sex + Rearing2 + Age + CurExp, random = ~animal + ID+ Task, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)


summary(m5x)
model.avg(m1x, m2x, m3x, m4x, m5x)
ms<-model.sel(m1x, m2x, m3x, m4x, m5x, rank = DIC)
ms
mavg<-model.avg(ms)
mavg
confint(mavg)
HPDinterval(m4x)

# THIS MODEL GIVES MEASURE OF WITHIN-INDIVIDUAL CONSISTENCY
rep1 <- MCMCglmm(RawScore ~ 1, random = ~animal, family = "gaussian",
                prior = prior1, pedigree = p2, verbose = FALSE, data = D, nitt = 100000, burnin = 10000, thin = 100)
summary(rep1)
rep.2<-(rep1$VCV[,"animal"]/ (rep1$VCV[,"animal"]+rep1$VCV[,"units"]))
posterior.mode(rep.2)

rep2 <- MCMCglmm(RawScore ~ 1, random = ~animal + ID, family = "gaussian",
                 prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)
summary(rep2)
rep.1<-(rep2$VCV[,"animal"]/ (rep2$VCV[,"animal"]+rep2$VCV[,"units"]))
rep.2<-(rep2$VCV[,"ID"]/ (rep2$VCV[,"ID"]+rep2$VCV[,"units"]))
posterior.mode(rep.2)

#phenotypic variance conditional on fixed effects
mVP<-rep2$VCV[,"animal"]+rep2$VCV[,"ID"]+rep2$VCV[,"units"]
#between individual variance
mIP<-rep2$VCV[,"animal"]+rep2$VCV[,"ID"]
#repeatability
posterior.mode(mIP/mVP)
#heritability
posterior.mode(rep2$VCV[,"animal"]/mVP)

posterior.mode(rep.1)
posterior.mode(rep.2)
mean(rep.1)
mean(rep.2)
HPDinterval(rep.1)

# calculate phenotypic variance conditional fixed effects
modelVP<-rep1$VCV[,"animal"]+rep1$VCV[,"ID"]+rep1$VCV[,"units"]
# calculate between individual variance
modelVA<-rep1$VCV[,"animal"]+rep1$VCV[,"ID"]
#repeatability
posterior.mode(modelVA/modelVP)



#with pr pi savex
m5.1x <- MCMCglmm(RawScore ~ sex + Rearing + Age + CurExp, random = ~animal + ID, family = "gaussian",
                  prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100, , pl=T, pr=T, saveX=T, saveZ=T)

m6x <- MCMCglmm(AdjScore ~ sex + Rearing + Age + CurExp, random = ~animal + ID, family = "gaussian",
                prior = prior2, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100, pl=T, pr=T, saveX=T, saveZ=T)
summary(m5x)

m5x <- MCMCglmm(AdjScore ~ sex + Rearing2 + Age + CurExp, random = ~animal + ID + study, family = "gaussian",
                prior = prior3, pedigree = p2, data = D, nitt = 100000, verbose = F, burnin = 10000, thin = 100)