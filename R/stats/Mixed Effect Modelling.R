###MIXED-EFFECT MODELLING###

str(Madeira)

library(lattice)

bwplot(SppRichness~Disturbance|Block, data=Madeira) #Can only check for crossed factors.

###Modelling###

library(nlme)

Model1<-lme(fixed=SppRichness~Disturbance*Productivity, random=~1|Block, data=Madeira) #This is the modelling of a mixed effect model. Splitting fixed and random effects.

###Diagnostics###

plot(resid(Model1)~fitted(Model1))
abline(h=0, lwd=2, lty=2, col="black")

#IN MEM NORMALITY OF ERRORS CANNOT BE CHECKED USING THE NORMALITY OF THE MODEL RESIDUALS. YOU HAVE TO CHECK NORMALITY OF ERRORS FOR EACH LEVEL OF RANDOM FACTOR!

qqnorm(Model1,~resid(.)|Block) #Gives a qqplot for every block (Random factor level).

###Results###

summary(Model1)
#This gives lots of result data. Its actually a regression output.

anova(Model1) #Can express output in anova style.

#To extract the amount of variance explained bythe random effect we have to use another package.

library(lme4)
Model2<-lmer(SppRichness~Disturbance+Productivity+(1|Block), data=Madeira)

summary(Model2)
#Looking at the "Random effects" section you can see how much of the unexplained variance is covered by "Blocks" and how much is still unexplained (residuals). In this case about 40% are covered (looking at the proportion of variances)