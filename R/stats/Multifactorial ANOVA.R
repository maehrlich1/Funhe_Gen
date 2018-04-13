###################
MULTI-FACTORIAL ANOVA
###################

###Housekeeping###

str(Madeira)
dotchart(Madeira$SppRichness)

dotchart(Madeira$SppRichness,groups=factor(Madeira$Productivity)) #Structures data by groups also

Madeira$Productivity<-factor(Madeira$Productivity, c("Ambient", "Low", "High")) #Sorts the factors, otherwise it would be alphabetically sorted which does not make sense.

boxplot(Madeira$SppRichness~Madeira$Productivity, ylim=c(0,15), xlab="Nutrient Availability", ylab='Species Richness', main="Productivity")

boxplot(Madeira$SppRichness~Madeira$Disturbance, ylim=c(0,15), xlab="Disturbance Level", ylab="Species Richness", main="Disturbance")

Means<-tapply(Madeira$SppRichness, interaction(Madeira$Disturbance,Madeira$Productivity), mean, na.rm=TRUE)#With this me can extract group means, the last term removes missing values. The interaction part tells R to calc means at the combination level, we do not want means of only one factor.

round(Means, digits=0)

###Check for Interactions###

library(lattice)

bwplot(Madeira$SppRichness~Madeira$Disturbance|Madeira$Productivity)#Can plot several boxplots so see whether there are interactions.

bwplot(Madeira$SppRichness~Madeira$Productivity|Madeira$Disturbance)#Can also reverse the independent variables.

#Can see some very slight interactions but they do not look significant.

interaction.plot(Madeira$Disturbance, Madeira$Productivity, Madeira$SppRichness)#The last one is the response variable. The first one the x-axis and the second the trace factor.
#Not recommended to use this plot since it does not show CI.

###Diagnostics###

plot(resid(Model1)~fitted(Model1))#Remember the fitted values are nothing but the group means. The residuals the difference between the single data points to their group mean.
abline(h=0, lwd=2, lty=2, col="black")
#Variances do not look perfectly homogeneous but it doesnt seem too bad.

fligner.test(Madeira$SppRichness~interaction(Madeira$Disturbance,Madeira$Productivity))
#Can test using FK but we have to look at the interaction, not at single factors.

hist(resid(Model1))#Looks ok!

shapiro.test(resid(Model1))#Looks good but remember Shapiro becomes very conservative/critical the bigger the number of observations. Works well up to 70 obs.

plot(cooks.distance(Model1), type="h") #No influential data points.

###Modelling###

Model1<-aov(Madeira$SppRichness~Madeira$Productivity*Madeira$Disturbance)#Same as one-way only we indicate interaction with an asterisk.

###Results###

summary(Model1)
#Check interaction first, if there is no significant interaction term we can simplify things and look at the other two factor individually.

Model2<-aov(Madeira$SppRichness~Madeira$Productivity+Madeira$Disturbance)#Can use a plus instead of asterisk. We are forcing R not to consider interactions.
#Can exclude the interaction term now that we know it is insignificant. Like this we have a more powerful test!

summary(Model2)
#Less degrees of freedom but only slightly higher sums of resids.

###Non-parametric multi-factorial ANOVA###

sheirer-ray-hare test #Not really used since we can use PERMANOVA which is more powerful.
