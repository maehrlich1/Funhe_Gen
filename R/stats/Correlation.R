#######Correlations########

#######Housekeeping########

str(Gannets)
length(Gannets$Food)#Should always indicate n in correlations

########Plotting########

dotchart(Gannets$Growth)
dotchart(Gannets$Food)

plot(Gannets$Food,Gannets$Growth)#This is a simple scatterplot of Food vs Growth. Use the , instead of the ~ because the ~ indicates a model! (Does work though...)

#####Assumptions#########

hist(Gannets$Growth)#Does not look normal.
shapiro.test(Gannets$Growth)#Shapiro suggests everything is fine 

hist(Gannets$Food)
shapiro.test(Gannets$Food)#Shapiro suggests non-normality. DATA NOT NORMAL DISTR!

hist(sqrt(Gannets$Food))#Try to transform to achieve better result...doesnt really improve
hist(log(Gannets$Food))#Also not very good. Nonetheless had we done a transformation we would have had to do the same to both variables!!

####Testing#####
#We use the Spearman (non-parametric) test since data is not normal

cor.test(Gannets$Growth, Gannets$Food, method="spearman")#Spearman test converts data to ranks, results in ties sometimes. (Equal data in same rank)
#rho is the equivalent of the parametric Pearsons R (-1<R<1) showing the strength of correlation. >0.7 suggests strong correlation
