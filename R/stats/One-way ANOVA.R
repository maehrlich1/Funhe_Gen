#########One-Way ANOVA#######

#############
#Housekeeping
#############

str(Plantgrowth)#Structure of data set incl. data quality levels

Plantgrowth$Nutrients<-as.factor(Plantgrowth$Nutrients)#Since R did not recognise the Nutrients as an integer/continuous variable. Need to change the format of the Nutrients column.

levels(Plantgrowth$Nutrients)#Now that R changed the structure it can show the levels of that column

table(Plantgrowth$Nutrients)#Shows number of replicates per category. Can be used to check for a balanced design

tapply(Plantgrowth$Biomass, Plantgrowth$Nutrients, mean) #Shows you the mean of the individual factor levels

############
#Plotting Graphs
################

dotchart(Plantgrowth$Biomass, groups=factor(Plantgrowth$Nutrients))#Second part tells R to group the data points according to the groups in our independent variable

boxplot(Plantgrowth$Biomass~Plantgrowth$Nutrients)#Normal boxplot of Biomass as a function of Nutrients

plot.design(Plantgrowth$Biomass~Plantgrowth$Nutrients)

#Can also do a mean-error plot using ggplot2

library(ggplot2)

library(Hmisc)

Figure1<-ggplot(Plantgrowth, aes(Nutrients, Biomass)) #This layer defines the axes

Figure1+
  stat_summary(fun.y=mean, geom="bar", colour="gray", size=1, fill="grey60")+
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", position=position_dodge(width=0.9), width=0.0, colour="black", linetype=1, size=0.8)+
  geom_point(colour="red")+
  labs(x="Nutrient addition in %", y="Biomass(g/m2)")+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(colour="black"),
        axis.title=element_text(size=20, colour="black"),
        axis.text=element_text(size=15))

#1st layer: Axes
#2nd layer: Adds the mean as an errorbar
#3rd layer: Adds the original data points
#4th layer: Adds labels
#5th layer: Format

############
#Modelling##

Model1<-aov(Plantgrowth$Biomass~Plantgrowth$Nutrients)#Creates the ANOVA model, same syntax as a graph.

######################
#Diagnostics##

plot(Model1)#Automatically produces important plots of the Model1

plot(resid(Model1)~fitted(Model1))
abline(h=0, lwd=2, lty=2, col="black") #Can also produce the residual plot manually. In ANOVA the fitte values are the sample means!

fligner.test(Plantgrowth$Biomass~Plantgrowth$Nutrients)#Fligner-Kileen test for homogeneity of variances. Same syntax as graphs.

hist(resid(Model1))#creates a histogram of the residuals of the ANOVA! To check for normality of errors

shapiro.test(resid(Model1))#normality test for the residuals. checking for normality. Only works for 10<n<70

plot(cooks.distance(Model1), type="h")#Takes out individual data points from the model and reruns. The discrepancy between the original and the modified is the Cook's Distance! Should be <1. "h" is for histogram style

cooks.distance(Model1)#Can ask for Cook's disatnce directly but its easier to visualise in graph form above.

#################
####Output##

summary(Model1)#Gives the test stats of the ANOVA

summary.lm(Model1)#Gives regression stats of the ANOVA incl. R-squared (percentage of explained variation)

#######################
#Retro-Power analysis##

power.anova.test(groups=6, n=10, between.var=2371534, within.var=50537, power=NULL, sig.level=0.05)#Retrospective power analysis to estimate variables given other parameters. One of them needs to be set to NULL to be calculated.

######################
#Post-hoc Test#######

#Essentially a multiple T-test

TukeyHSD(Model1)#Lists all possible 1-1 comparisons to check for significant difference between each and every group

plot(TukeyHSD(Model1), las=2) #Shows all of the different group comparisons. If the CI of of the group means include the zero lines, we conclude they are not significantly different. The las turns the labes by 90 deg.

########################
###NON-PARAMETRIC OW-ANOVA####

kruskal.test(Plantgrowth$Biomass~Plantgrowth$Nutrients)#Non-parametric equivalent of OW-ANOVA. Less powerful due to non-parametrics

####################
#Welch adjustment for heteroscedacity###

oneway.test(Plantgrowth$Biomass~Plantgrowth$Nutrients)#OW-ANOVA + Welch adjustment for heterogeneous variances. Otherwise syntax same as ANOVA