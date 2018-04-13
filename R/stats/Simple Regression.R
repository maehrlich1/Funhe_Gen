########Regression Analysis#########

#######Housekeeping########

str(Predator)

####Plotting#####

dotchart(Predator$Predators)
dotchart(Predator$Prey)

plot(Predator$Prey~Predator$Predators)

abline(lm(Predator$Prey~Predator$Predators))#This draws a line through the datapoints. lm stands for linear model. Just visual, no modelling yet

########Modelling###########

Model1<-lm(Predator$Prey~Predator$Predators)#Creates the linear regression model.

Fitted<-fitted(Model1) #This makes a vector with the fitted values of our model.

Fitted<-as.integer(Fitted) #This rounds the numbers in the vector to integers.

plot(Fitted~Predator$Predators) #This is essentially the fitted line of our model.

######Assumptions##########

plot(Model1)
#Residuals do not seem very homogeneous in variance but its not too bad.
#Residuals seem normally distributed

hist(resid(Model1))#Residuals seem normally distributed
shapiro.test(resid(Model1)) #Shapiro agrees...

plot(cooks.distance(Model1), type="h")#Cooks distances not >1 so it looks like there are no highly influential data points...

library(car)
durbinWatsonTest(Model1) #This is the test for autocorrelation of residuals. If the D-W statistic is between 1-3 there is no problem with autocorrelation.

###Model Results####

summary(Model1)
#Model seems to fit the data very well (High R-squared)
#Intercept at 30 with standard error. t-test checks whether intercept is sig. different from zero, comes with a p-value.
#Gradient is -1.2 with standard error. t-test checks whether gradient is sig. different from zero, comes with a p-value.
#F stat to test whether our model is sig different to the null model
#F stat comes with a p-value