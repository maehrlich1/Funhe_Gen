############
#Plotting Graphs
################

dotchart(Plantgrowth$Biomass, groups=factor(Plantgrowth$Nutrients))#Second part tells R to group the data points according to the groups in our independent variable

boxplot(Plantgrowth$Biomass~Plantgrowth$Nutrients)#Normal boxplot of Biomass as a function of Nutrients

par(bty="l")#Removes boxing around graph. Needs to be run before the graphing command!

boxplot(Plantgrowth$Biomass~Plantgrowth$Nutrients,
        xlab="Nutrients", ylab="Biomass", main="ANOVA Exmp", col="turquoise",
        ylim=c(2000,4500), boxwex=0.5, whisklty=1, staplecol=FALSE, outline=FALSE,
        yaxp=c(2000,4500,5))#Same as above +Labels +Colours +Deliminations +Box Size +Line Type +Line ends +Remove Outliers +No. of ticks

colors()#Lists all of the available colours

plot(Predator$Prey~Predator$Predators, pch=16, col="red")#Can change the point character

plot.design(Plantgrowth$Biomass~Plantgrowth$Nutrients) #Can show the design of the experiment

#####Error Bars######

library(Hmisc)

errbar(df$X, df$Y, df$Y + df$SD, df$Y - df$SD) # Syntax is (x,y,yplus,yminus) where yplus and yminus are the upper and lower limits of the error bars.

###GGPLOT2###

library(ggplot2)

library(Hmisc)

Figure1<-ggplot(Plantgrowth, aes(Nutrients, Biomass)) #This layer defines the axes

Figure1+
  stat_summary(fun.y=mean, geom="bar", colour="gray", size=1, fill="grey60")+
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", position=position_dodge(width=0.9), width=0.0, colour="black", linetype=1, size=0.8)+
  geom_point(colour="red")+
  geom_smooth(method=lm)+
  theme_bw()+
  labs(x="Nutrient addition in %", y="Biomass(g/m2)")+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line=element_line(colour="black"),
        axis.title=element_text(size=20, colour="black"),
        axis.text=element_text(size=15))

#First layer: Axes
#Second layer: Adds the mean as an errorbar
#Third layer: Adds the original data points
#Fourth layer: Adds a smoother
#Fifth layer: Changes theme to BW
#6th layer: Labels
#7th layer: Formatting