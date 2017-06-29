str(D)
#create a new dataframe: names should match those of the originary dataframe
# Create 'dummy' data frame for prediction
# - all IDs, all opponent sizes,
# - block set at 0 (as blocks are coded as -0.5, 0.5)
# - mean body size and assay repeats
df_ri_ind <- expand(D,
                    ID, RawScore, sex)
df<-expand(D)
# Get predicted values based on RR model and dummy data frame
# and using random effects structure as in the model
df_ri_ind$fit <- predict(lmer_a, newdata = df_ri_ind, re.form = NULL)
# Plot predictions and overlay original data points
ggplot(df_ri_ind, aes(x = opp_size, y = fit, group = ID)) +
  geom_line() +
  geom_point(data = df_plast,
             aes(y = aggression),
             alpha = 0.3) +
  scale_x_continuous(breaks = c(-1,0,1)) +
  theme_classic() +
  facet_wrap(~ID)





ggplot(D, aes(sex, RawScore)) + geom_point()
curve(predict(m2,data.frame(bodysize=x),type="resp"),add=TRUE) # draws a curve based on prediction from logistic regression model
?curve

newdat <- expand.grid(RawScore = plotsc, sex = c("Female", "Male"))
newdat$pred = predict(m2, newdata = newdat)


#META PLOTS
D2<-D
D2<-subset(D, TotalExp > 3) ; D$TotalExp<-drop.levels(D$TotalExp) ; str(D)

aggmedian<-aggregate(D2[,"RawScore"], list(D2$ID, D2$sex, D2$Rearing, D2$Rearing2), median); aggmedian$ID<-aggmedian$Group.1; aggmedian$MedianScore<-aggmedian$x; aggmedian$sex<-aggmedian$Group.2plotagg<-ggplot(aggmedian, aes(sex, MedianScore)) 
aggmedian$Rearing <- aggmedian$Group.3 ; aggmedian$Rearing2<- aggmedian$Group.4 ; aggmedian$Rearing3<- aggmedian$Group.5
  scale_y_continuous(limits=c(-0.05,1))+
  ylab("Score") +
  xlab("sex") +
  geom_line(aes(colour=ID), size=2 ) + 
  scale_shape_manual(values=c(22,21)) + 
  geom_jitter(width = 0.3, height = 0.01, size = 3) + 
  theme(legend.position="none", text = element_text(size=20))
plotagg

plotagg<-ggplot(aggmedian, aes(, MedianScore)) +
  scale_y_continuous(limits=c(-0.05,1))+
  ylab("Score") +
  xlab("sex") +
  geom_line(aes(colour=ID), size=2 ) + 
  scale_shape_manual(values=c(22,21)) + 
  geom_jitter(width = 0.3, height = 0.01, size = 3) + 
  theme(legend.position="none", text = element_text(size=20))
plotagg

plot(aggmedian$Rearing, aggmedian$MedianScore)
plot(D2$Rearing, D2$RawScore)
plot(D2$Rearing2, D2$RawScore)
plot(D2$Rearing3, D2$RawScore)
boxplot(sex ~ MedianScore, data = aggmedian)

aggmedian<-aggregate(D[,"RawScore"], list(D$ID, D$sex), median)
aggmedian2<-aggregate(D[,"RawScore"], list(D$ID, D$sex), sd)
aggmedian$sd<-aggmedian2$sd
aggmedian$ID<-aggmedian$Group.1; aggmedian$MedianScore<-aggmedian$x; aggmedian$sex<-aggmedian$Group.2
plotagg<-ggplot(aggmedian, aes(sex, MedianScore)) +
  scale_y_continuous(limits=c(-0.05,1))+
  ylab("Score") +
  xlab("sex") +
  geom_line(aes(colour=ID), size=2 ) + 
  scale_shape_manual(values=c(22,21)) + 
  geom_jitter(width = 0.3, height = 0.01, size = 3) + 
  theme(legend.position="none", text = element_text(size=20))
plotagg

#plotting coefficients
cm1<-caterplot(m1$Sol); cm1v<-caterplot(m1$VCV)
cm2<-caterplot(m2$Sol); cm2v<-caterplot(m2$VCV)
cm2.5<-caterplot(m2.5$Sol); cm2.5v<-caterplot(m2.5$Sol)
cm3<-caterplot(m3$VCV); cm3v<-caterplot(m3$VCV)
cm3.5<-caterplot(m3.5$Sol); cm3.5v<-caterplot(m3.5$VCV)
cm4<-caterplot(m2$Sol); cm4v<-caterplot(m4$VCV)
cm4.5<-caterplot(m4.5$Sol);cm4.5v<-caterplot(m4.5$VCV)
cm1; cm1v; cm2; cm2v; cm2.5; cm2.5v; cm3; cm3v; cm3.5; cm3.5v; cm4; cm4v; cm4.5; cm4.5v
caterplot(m2$VCV)
#plots diagnostics
xyplot(as.mcmc(m2$Sol))

#plotting predictions?
predPlot + geom_point(data = d.predNew, stat="identity", position = position_dodge(width=0.3), size = 2.8) + 
  geom_errorbar(limits, width = 0.08, position = position_dodge(width=0.3)) +
  geom_hline(aes(yintercept=0.5), linetype="dashed", show.legend=FALSE) + 
  theme_bw() + theme(text = element_text(size=12), axis.title.x=element_blank(), axis.title.y=element_text(margin=margin(0,12,0,0))) + 
  ylab("Proportion Chose Social Source") +
  scale_y_continuous(limits=c(0,1), expand = c(0,0)) +
  scale_x_discrete(limits=c("Control", "Social Risky","Asocial Risky")) 


personalityPlot <- ggplot(D, aes(D$RawScore, fill = sex)) 
personalityPlot + scale_fill_grey(start = 0.1, end = 0.9) + geom_density(alpha = 0.2) + theme_bw() + 
  theme(text = element_text(size=12), axis.title.y=element_text(margin=margin(0,12,0,0))) +
  scale_y_continuous(limits=c(0,0.09), expand = c(0,0)) +
  scale_x_continuous(limits=c(1,2.5), expand= c(0,0)) +
  xlab("\nRawScore") + ylab("Density") 
personalityPlot

personalityMean = tapply(D$RawScore, list(D$sex),mean); personalityMean
personalitySD = tapply(D$RawScore, list(D$sex),sd); personalitySD

tryingPlot <- ggplot(D, aes(Rearing, RawScore, colour = sex))
tryingPlot + geom_point(data = D, size = 3.5) + 
  theme_bw() + theme(text = element_text(size=26)) + ylab("Score") + ylim(0,1) +
  scale_x_discrete(limits=c("M", "N", "W")) +
  ggtitle("Raw scores")
tryingPlot
