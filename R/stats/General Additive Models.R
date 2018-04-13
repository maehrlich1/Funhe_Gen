###########################
General Additive Models
###########################

###Housekeeping###

str(Plantgrowth)

plot(Plantgrowth$Biomass~Plantgrowth$Nutrients)

abline(lm(Plantgrowth$Biomass~Plantgrowth$Nutrients))

#Actually the data looks like a log function/saturation curve. Might try some nonlinear model...

###Nonlinear Modelling###

library(mgcv)

Model1<-gam(Biomass~s(Nutrients, k=3), data=Plantgrowth) #this gives us a non-linear model. s indicates a smoothing of the independent variable. k indicates the complexity of the nonlinear model, try to minimise. data calls the dataframe we want to use.

###Plotting###

plot(Model1) #gives the model and the CI. For some reason the y-axis is strange?? Use ggplot...

library(ggplot2)

Figure1<-ggplot(Plantgrowth, aes(x=Nutrients, y=Biomass))

Figure1+
  geom_point(colour="red")+
  geom_smooth(method=lm)+
  theme(panel.background=element_blank(), axis.line=element_line(colour="black"), panel.grid=element_blank())

#Graph of 3 layers: axes, data points and smoother. The rest is formatting.
#But we want to use a non-linear smoother!Just change the method...

Figure2<-ggplot(Plantgrowth, aes(x=Nutrients, y=Biomass))

Figure2+
  geom_point(colour="red")+
  geom_smooth(method=loess)+
  theme(panel.background=element_blank(), axis.line=element_line(colour="black"), panel.grid=element_blank())

###Diagnostics###

plot(resid(Model1)~fitted(Model1)) #This is the residual plot also given in the plotting part of the linear modelling. Checks homoscedasticity.

hist(resid(Model1)) #Normality

plot(cooks.distance(Model2), type="h") #No highly influential data points

summary(Model1) #Significant model!

coef(Model1) #Can get the coefficients of the nonlinear function? Unsure...

###Compare to a LINEAR MODEL###

Model2<-lm(Plantgrowth$Biomass~Plantgrowth$Nutrients)

summary(Model2)

plot(Model2) #Homoscedacity

hist(resid(Model2)) #Normality of resid

durbinWatsonTest(Model2) #No autocorrelation

plot(cooks.distance(Model2), type = "h") #No data points are extremely influential.

#Linear model could also be used! Need to check which one is better.

AIC(Model1, Model2) #This index judges the complexity vs. explanatory power of a model. The smaller the better! Non-linear seems better!
