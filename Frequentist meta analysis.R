#library(devtools); install_github("lme4","lme4")
install.packages('sp')
install.packages("AICcmodavg", dependencies = TRUE)


library(AICcmodavg); library(sjstats); library(blmeco)
library(lme4);library(arm);library(sjPlot);library(glmulti);library(bblme);library(MuMIn);;library(Rcpp);library(AICcmodavg);library(coefplot);library(ggplot2); library(multcomp)
library(kinship2) ; library(MCMCglmm) ; library(pedigree) ; library(pedigreemm) ; library(MasterBayes) ; library(rethinking) ; library(MuMIn) ; library(AICcmodavg) ; library(sjstats)
#datafile setup and variable correction
d <- as.matrix(read.table("~/Meta-project/MCMC_GLMM.txt", header=TRUE, na.strings="NA", sep = "\t")); d <- as.data.frame(d); d[is.na(d)] <- NA;str(d)
d$RawScore<-as.numeric(as.character(d$RawScore)); d$AdjScore<-as.numeric(as.character(d$AdjScore)); d$TotalExp<-as.numeric(as.character(d$TotalExp))
d$Age<-as.numeric(as.character(d$Age)); d$SLBinary<-as.numeric(as.character(d$SLBinary)); d$CurExp<-as.numeric(as.character(d$CurExp)); str(d)

#eliminate NA lines to streamline the dataset
d<-subset(d, !is.na(RawScore))
str(d)
# Let's subset out all of the ghost conditions
D<-d
D<-subset(D,!(Context == "Ghost")) ; D$Context<-drop.levels(D$Context) ; str(D)
# And then let's remove unknown rearing histories
D<-subset(D,!(Rearing == "U")) ; D$Rearing<-drop.levels(D$Rearing) ; str(D)
#D<-subset(D,!(Age.group == "1")) ; D$Age.group<-drop.levels(D$Age.group) ; str(D) ; D$Age.group<-as.numeric(as.character(D$Age.group))
#subsetting individuals
#D<-subset(D, TotalExp > 2) ; D$TotalExp<-drop.levels(D$TotalExp) ; str(D)
str(D)

m1<- lmer (RawScore ~ 1 + (1 | ID) , data = D, REML = F)
m2<- lmer (RawScore ~ sex + (1 | ID:study) , data = D, REML = F)
m3<- lmer (RawScore ~ 1 +sex + Rearing + (1 | ID) , data = D, REML = F)
m4<- lmer (RawScore ~ 1 +Rearing +sex + Age + (1 | ID) , data = D, REML = F)
m5<- lmer (RawScore ~ 1 +Rearing + sex + Age + CurExp + (1 | ID), data = D, REML = F)
m7<- lmer (RawScore ~ 1 +Rearing2 + sex + Age + CurExp + (1 | ID) + (1|study2), data = D, REML = F)
summary(m7)


mb<- glmer (RawScore ~ 1 +Rearing2 + sex + Age + (1 | ID) + (1|study2), family=binomial, data = D,weights = CurExp)
weights = size
summary(mb)
summary(m7)
confint(m5)
dm7<-dredge(m7)
dm7
dm7s<-subset(dm7, delta < 4)
dmavg<-model.avg(dm7s);  dmavg
dmavg<-model.avg(dm7);  dmavg
model.sel(dm7)
summary(dmavg)
confint(dmavg)
HPDinterval(m5)

icc(m7)

m1<- lmer (RawScore ~ (1 | ID) + (1|study) , data = D, REML = F)
m2<- lmer (RawScore ~ sex + (1 | ID) + (1|study) , data = D, REML = F)
m3<- lmer (RawScore ~ sex + Rearing + (1 | ID) + (1|study) , data = D, REML = F)
m4<- lmer (RawScore ~ Rearing +sex + Age + (1 | ID) + (1|study) , data = D, REML = F)
m5<- lmer (RawScore ~ Rearing + sex + Age + CurExp + (1 | ID) + (1|study), data = D, REML = F)
m6<- lmer (RawScore ~ Rearing2 + sex + Age + CurExp + (1 | ID) + (1|study), data = D, REML = F)
m7<- lmer (RawScore ~ Rearing2 + sex + Age + CurExp + Context + (1 | ID) + (1|study), data = D, REML = F)
m8<- lmer (RawScore ~ Rearing2 + sex + Age + CurExp + Context + Task + (1 | ID) + (1|study), data = D, REML = F)

m1<- lmer (RawScore ~ (1 | ID) + (1|study) , data = D, REML = F)
m2<- lmer (RawScore ~ sex + (1 | ID) + (1|study) , data = D, REML = F)
m3<- lmer (RawScore ~ sex + Rearing + (1 | ID) + (1|study) , data = D, REML = F)
m4<- lmer (RawScore ~ Rearing +sex + Age + (1 | ID) + (1|study) , data = D, REML = F)
m5<- lmer (RawScore ~ Rearing + sex + Age + CurExp + (1 | ID) + (1|study), data = D, REML = F)
m6<- lmer (RawScore ~ Rearing2 + sex + Age + CurExp + (1 | ID) + (1|study) + (1|Task) + (1| Context), data = D, REML = F)
coef(m6)

dispersion_glmer(m7)

m1<- lmer (RawScore ~ (1 | ID) + (1|Context) + (1|Task) , data = D, REML = F)
m2<- lmer (RawScore ~ sex + (1 | ID) + (1|Context) + (1|Task) , data = D, REML = F)
m3<- lmer (RawScore ~ sex + Rearing + (1 | ID) + (1|Context) + (1|Task) , data = D, REML = F)
m4<- lmer (RawScore ~ Rearing +sex + Age + (1 | ID) + (1|Context) + (1|Task) , data = D, REML = F)
m5<- lmer (RawScore ~ Rearing + sex + Age + CurExp + (1 | ID) + (1|Context) + (1|Task), data = D, REML = F)
m6<- lmer (RawScore ~ sex + Rearing2 + (1 | ID) + (1|Context) + (1|Task) , data = D, REML = F)
m7<- lmer (RawScore ~ Rearing2 +sex + Age + (1 | ID) + (1|Context) + (1|Task) , data = D, REML = F)
m8<- lmer (RawScore ~ Rearing2 + sex + Age + CurExp + (1 | ID) + (1|Context) + (1|Task), data = D, REML = F)


m1<- lmer (AdjScore ~ (1 | ID) + (1|study) , data = D, REML = F)
m2<- lmer (AdjScore ~ sex + (1 | ID) + (1|study) , data = D, REML = F)
m3<- lmer (AdjScore ~ sex + Rearing + (1 | ID) + (1|study) , data = D, REML = F)
m4<- lmer (AdjScore ~ Rearing +sex + Age + (1 | ID) + (1|study) , data = D, REML = F)
m5<- lmer (AdjScore ~ Rearing + sex + Age + CurExp + (1 | ID) + (1|study), data = D, REML = F)
m6<- lmer (AdjScore ~ Rearing2 + sex + Age + CurExp + study + (1 | ID), data = D, REML = F)
summary(m6)
summary(glht(m4, mcp(Rearing="Tukey")))
con.models<-list(m1,m2,m3,m4, m5)
con.models.names<-c("m1", "m2", "m3", "m4", "m5")
con.models<-list(m1,m2,m3,m4, m5, m6, m7, m8)
con.models.names<-c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8")
#table for AICc
aictab(cand.set=con.models, modnames=con.models.names)
model.sel(con.models)
mavgf<-model.avg(con.models)
mavgf
confint(mavgf)
summary(con.models)
plot(m4)

#anova(m1, m1.1, m1.2, m1.3, m1.4,m1.5,m1.6,m1.7,m1.8)
#model averaged estimates to take into account model uncertainty - also CI's!!!
m.avg.sex<-modavg(cand.set=con.models, modnames=con.models.names, parm="sexM")
m.avg.age2<-modavg(cand.set=con.models, modnames=con.models.names, parm="Age")
m.avg.age3<-modavg(cand.set=con.models, modnames=con.models.names, parm="Age.group")
m.avg.age4<-modavg(cand.set=con.models, modnames=con.models.names, parm="Age.group4")
m.avg.rearingN<-modavg(cand.set=con.models, modnames=con.models.names, parm="RearingN")
m.avg.rearingW<-modavg(cand.set=con.models, modnames=con.models.names, parm="RearingW")
m.avg.rearing2<-modavg(cand.set=con.models, modnames=con.models.names, parm="Rearing2N")
m.avg.curexp<-modavg(cand.set=con.models, modnames=con.models.names, parm="CurExp")

m.avg.sex
m.avg.age2
m.avg.age3
m.avg.rearingN
m.avg.rearingW
m.avg.curexp

con.models<-list(m1,m2,m3,m4, m5, m6)
con.models.names<-c("m1", "m2", "m3", "m4", "m5", "m6")
#table for AICc
aictab(cand.set=con.models, modnames=con.models.names)

summary(con.models)
plot(m4)

#odds of intercept
OR<-logistic(c(-14.8, -23.07, -6.52))
round(OR, digits=3)
#odds of main effect
logistic(c((-14.8 + 15.28), (-6.52+6.61), (-23.07+23.94)))

m1<- glmer (SLBinary ~ 1+ (1 | ID/study), data = D,  family=binomial)
m2<- glmer (SLBinary ~ 1+ sex + (1 | ID/study), data = D,  family=binomial)
m3<- glmer (SLBinary ~ 1+ sex + Rearing + (1 | ID/study), data = D, family=binomial)
m4<- glmer (SLBinary ~ 1+ Rearing +sex + Age + (1 | ID/study), data = D,  family=binomial)
m5<- glmer (SLBinary ~ 1+ Rearing + sex + Age + CurExp + (1 | ID/study), data = D, family=binomial)
m6<- glmer (SLBinary ~ 1+ Rearing + sex + Age +  (1 | ID/study), data = D, family=binomial)

m1<- glmer (SLBinary ~ 1+ (1 | ID)+ (1|study), data = D,  family=binomial)
m2<- glmer (SLBinary ~ 1+ sex + (1 | ID)+ (1|study), data = D,  family=binomial)
m3<- glmer (SLBinary ~ 1+ sex + Rearing + (sex | ID)+ (1|study), data = D, family=binomial)
m4<- glmer (SLBinary ~ 1+ Rearing +sex + Age + (1 | ID)+ (1|study), data = D,  family=binomial)
m5<- glmer (SLBinary ~ 1+ Rearing + sex + Age + CurExp +  (1 | ID)+ (1|study), data = D, family=binomial)
m6<- glmer (SLBinary ~ 1+ Rearing2 + sex + Age +  (1 | ID) + (1|study2), data = D, family=binomial)
summary(m6)

ranef(m3)
m1<- glmer (SLBinary ~ (1 | ID), data = D,  family=binomial)
m2<- glmer (SLBinary ~ sex + (1 | ID), data = D,  family=binomial)
m3<- glmer (SLBinary ~ sex + Rearing + (1 | ID), data = D, family=binomial)
m4<- glmer (SLBinary ~  Rearing +sex + Age + (1 | ID), data = D,  family=binomial)
m5<- glmer (SLBinary ~ Rearing + sex + Age + CurExp +  (1 | ID), data = D, family=binomial)

summary(m6)
icc(m6)
library(AICcmodavg)
con.models<-list(m1,m2,m3,m4, m5)
con.models.names<-c("m1", "m2", "m3", "m4", "m5")
#table for AICc
aictab(cand.set=con.models, modnames=con.models.names)
model.sel(con.models)
mavgf<-model.avg(con.models)
mavgf
confint(mavgf)
summary(con.models)
dispersion_glmer(m4)

#to get the fitted average reaction time per subject
reaction_subject <- fixef(m4) + ranef(m4)$ID 
reaction_subject$ID<-rownames(reaction_subject)
names(reaction_subject)[1]<-"Intercept"
reaction_subject <- reaction_subject[c(3,1)]
#plot
ggplot(reaction_subject,aes(x=ID,y=Intercept))+geom_point()


#the next line put all the estimated intercept and slope per
#subject into a dataframe
reaction_slp <- as.data.frame(t(apply(ranef(m4)$ID,
                                      1,function(x) fixef(m4) + x)))
#to get the predicted regression lines we need one further
#step, writing the linear equation: Intercept + Slope*Days
#with different coefficient for each subject
pred_slp <- melt(apply(reaction_slp,1,function(x) x[1] + x[2]*0:1),
                 value.name = "SLBinary")
#some re-formatting for the plot
names(pred_slp)[1:2] <- c("sex","ID")
pred_slp$sex <- pred_slp$sex - 1
pred_slp$ID <- as.factor(pred_slp$ID)
pred_slp$SLBinary<- as.Error: Discrete value supplied to continuous scale(pred_slp$value)

ggplot(pred_slp,aes(x=sex,y=value,color=ID))+
  geom_line()+
  geom_point(data=D,aes(x=sex,y=SLBinary))
facet_wrap(~ID,nrow=3)


#anova(m1, m1.1, m1.2, m1.3, m1.4,m1.5,m1.6,m1.7,m1.8)
#model averaged estimates to take into account model uncertainty - also CI's!!!
m.avg.sex<-modavg(cand.set=con.models, modnames=con.models.names, parm="sexM")
m.avg.age<-modavg(cand.set=con.models, modnames=con.models.names, parm="Age")
m.avg.rearingN<-modavg(cand.set=con.models, modnames=con.models.names, parm="RearingN")
m.avg.rearingW<-modavg(cand.set=con.models, modnames=con.models.names, parm="RearingW")
m.avg.curexp<-modavg(cand.set=con.models, modnames=con.models.names, parm="CurExp")

m.avg.sex
m.avg.age
m.avg.rearingN
m.avg.rearingW
m.avg.curexp
