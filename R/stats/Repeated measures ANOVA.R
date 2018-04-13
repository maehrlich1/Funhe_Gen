###Repeated-measures ANOVA###

###Housekeeping###

str(Littorina)

dotchart(Littorina$Eggs)

interaction.plot(Littorina$Time, factor(Littorina$Group), Littorina$Eggs, ylab="Mean of Eggs", xlab="Time")
#Not much interaction. Not much of an effect. Rhythmic egg release over time.

library(lattice)
library(nlme)

Snails<-groupedData(Eggs~Time|ID, data=Littorina, outer=~Group)#Prepare data for plotting.
plot(Snails)

plot(Snails, outer=TRUE)
#Just different plots to visualise the data

###Modelling###

library(nlme)

Model1<-lme(Eggs~factor(Group)*factor(Time), random=~1|ID, data=Littorina) #To add more than 1 random factor you can use the syntax: list(~1|whatever, ~1|whatever)

###Diagnostics###

plot(resid(Model1)~fitted(Model1))
abline(h=0, lwd=2, lty=2, col="black")
#Doesnt look too bad.

qqnorm(Model1, ~resid(.)|ID)
#Also looks good.

##RESULTS###

summary(Model1)
#This gives the regression output.

anova(Model1)
#This gives the anova output. The interaction is insignificant and so is the effect of group. Nonetheless the seems to be an effect of time on egg number.