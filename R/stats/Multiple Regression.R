###Multiple Regression###

###Housekeeping###

str(Eelgrass)

par(mfrow=c(2,2)) #This splits the graph output window into 4 parts. Easy visualisation.

dotchart(Eelgrass$Growth, main="Growth")
dotchart(Eelgrass$Light, main="Light")
dotchart(Eelgrass$Phosphate, main="Phosphate")
dotchart(Eelgrass$Epibionts, main="Epibionts")

par(mfrow=c(1,1))

pairs(Eelgrass, panel=panel.smooth) #Produces pair plots to see relationships between all variables. The panel command adds smoothers to all of them.

###Diagnostics###

#Check for colinearity of predictors. Need to compare all with all.

plot(Eelgrass$Phosphate~Eelgrass$Light)
abline(lm(Eelgrass$Phosphate~Eelgrass$Light))
cor.test(Eelgrass$Phosphate,Eelgrass$Light) #Pearson test to check for correlation, if there is this is problematic for modelling!

plot(Eelgrass$Phosphate~Eelgrass$Epibionts)
abline(lm(Eelgrass$Phosphate~Eelgrass$Epibionts))
cor.test(Eelgrass$Phosphate,Eelgrass$Epibionts) #Pearson test to check for correlation, if there is this is problematic for modelling!

plot(Eelgrass$Light~Eelgrass$Epibionts)
abline(lm(Eelgrass$Light~Eelgrass$Epibionts))
cor.test(Eelgrass$Light,Eelgrass$Epibionts) #Pearson test to check for correlation, if there is this is problematic for modelling!

#Can do pair plots and attach correlation coefficients! Correlation coefficients are plotted according to their absolute value. Does NOT give p-values and therefore can NOT check for colinearity.

panel.cor<-function(x,y, ...)
{
  par(usr=c(0,1,0,1))
  txt<-as.character(format(cor(x,y), digits=2))
  text(0.5, 0.5, txt, cex=6*abs(cor(x,y)))
}
pairs(Eelgrass, upper.panel=panel.cor)

#Check for interactions using coplots. Need to compare all with all. The panels command is producing the trendlines.
library(lattice)

coplot(Growth~Light|Phosphate, data=Eelgrass, panel=function(x,y,...){tmp<-lm(y~x)
                                                                      abline(tmp)
                                                                      points(x,y)})

#The grey bars at the top show which plot covers which Phosphate range. R tries to even out the number of samples per plot.

coplot(Growth~Light|Epibionts, data=Eelgrass, panel=function(x,y,...){tmp<-lm(y~x)
                                                                      abline(tmp)
                                                                      points(x,y)})

#The same for Epibionts.

coplot(Growth~Phosphate|Epibionts, data=Eelgrass, panel=function(x,y,...){tmp<-lm(y~x)
                                                                      abline(tmp)
                                                                      points(x,y)})

#The last combination.

#Can also check for interactions using a tree plot.

library(tree)

Model1<-tree(Growth~., data=Eelgrass)

plot(Model1)
text(Model1) #Adds text to the graph.

#Graph shows light is the most relevant predictor. At high light there is an interaction with Phosphate, in low light there isnt. In high light and low phosphate the epibionts also play a role.

###Modelling###

#Start with the most complex/maximal model.
Model1<-lm(Growth~Light*Phosphate*Epibionts+I(Light^2)+I(Phosphate^2)+I(Epibionts^2), data=Eelgrass)

summary(Model1)
#Almost 60% explained variation. P-value is highly significant.

#We need to eliminate predictors which with turn our model into a more powerful one. Follow set of rules:
#1 Eliminate the highest order interaction if insignificant. TRUE!

Model2<-update(Model1, ~.-Light:Phosphate:Epibionts, data=Eelgrass)# The update command allows you to modify an existing model. Inthis case it removes the 3 way interaction.

summary(Model2)
#The R value is slightly better but its not very reliable to compare between models. Use an anova approach.

anova(Model1, Model2)
#The p-value suggests the two models are not significantly different. Therefore you can safely eliminate the three-way interaction.

#Now we want to eliminate more predictors. We have to take out the interaction terms in parallel since they are at the same hierarchical level.

Model3a<-update(Model2, ~.-Phosphate:Epibionts, data=Eelgrass)#This removes the P:E interaction
Model3b<-update(Model2, ~.-Light:Epibionts, data=Eelgrass)#So on...
Model3c<-update(Model2, ~.-Light:Phosphate, data=Eelgrass)

#We now compare all of them to the previous model (Model2).

anova(Model2, Model3a)
anova(Model2, Model3b)
anova(Model2, Model3c)

#They are all insignificant so we are allowed to remove the interaction terms!

Model4<-update(Model3a,~.-Light:Epibionts-Light:Phosphate, data=Eelgrass) #This creates our new model where all of the interaction terms are removed.

summary(Model4) #Remaining predictors become more significant since model aquires higher resolution/power!

#now we can do the same for the quadratic terms.

Model5a<-update(Model4,~.-I(Light^2), data=Eelgrass)
Model5b<-update(Model4,~.-I(Phosphate^2), data=Eelgrass)
Model5c<-update(Model4,~.-I(Epibionts^2), data=Eelgrass)

anova(Model4, Model5a)
anova(Model4, Model5b)
anova(Model4, Model5c)

#The quadratic light term is significant! We cannot take this out! Nonetheless we can eliminate the other two!

Model6<-update(Model5b,~.-I(Epibionts^2), data=Eelgrass)

summary(Model6)
#Light is definitely the driver. This is the end point regarding removal of predictors. Could try to remove E and P but they are very close to significance.

#Can now see whether our simplification really improved the model.

AIC(Model1,Model6) #The Akaike's Information Criterion. The lower the number the better.
#It assesses both the explanatory power and the complexity of the models. Something is strange, the numbers are similar...

###Post-diagnostics###

#Checking for colinearity
library(car)
vif(Model6) #This gives the variance inflation factor. If they are larger 5 there is a problem! In this case light is a problem, might be useful to remodel without phosphate!

#Homoscedasticity, normality of errors and influential data points tested as in simple regression.

#Autocorrelation between residuals

library(car)
durbinWatsonTest(Model6)
#Should not be smaller than 1.