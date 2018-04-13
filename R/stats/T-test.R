#####FAKE DATA SET 31.10.14##########

###Creating a Data Set###

Lux<-c(rep("A",5),rep("B",5),rep("C",5),rep("D",5))    #Generates labels for treatment levels

GroupA<-rnorm(5,15,1)#Creates 5 random numbers with a mean of 15 and a s.d. of 1#
GroupB<-rnorm(5,12,1)
GroupC<-rnorm(5,8,1)
GroupD<-rnorm(5,6,1)

Depth<-c(GroupA,GroupB,GroupC,GroupD)#Collates the groups as the dependent variable vector

Data<-data.frame(Lux,Depth)#Generates a dataframe (matrix) combining both ind. and dep. vectors

Data$Depth<-as.integer(Data$Depth)#The $ sign identifies the vector you want to work with from the dataframe

####################
#####T-Test######

#T-test can only compare two groups directly against each other. Therefore we will just take two.

GroupA<-Data[Data$Lux =="A",2]#Subsets only the Lux group "A". The 2 after the comma indicates we want the number column of the "Data" data frame only.
GroupB<-Data[Data$Lux =="B",2]#See above

boxplot(GroupA, GroupB)

hist(GroupA)#Creates a histogram of GroupA
shapiro.test(GroupA)#Tests the data for normality. Only use for data sets >10 or <70. Null hypothesis assumes normality, therefore a significant p-value suggests non-normality.
var(GroupA)#Gives the variance of GroupA

hist(GroupB)
shapiro.test(GroupB)
var(GroupB)

 #Rule of thumb for homogeneity of variances: Divide the variance of one group by the variance of the other and it should be lower than 2.

t.test(GroupA,GroupB)#This test includes the Welch adjustment which compensates for inhomogeneity of variances. Consumes some test power. Nonetheless homogeneity of variance is more important than normality.

t.test(GroupA,GroupB,var=TRUE)#This is just the normal t-test. Dont apply the Welch adjustment

###############
#####Power analysis######

#Power analysis shows the test power either pro- or retrospectively. First of all we need the components within the t-test formula.

mean(GroupA)-mean(GroupB)#The difference in means is required.

sd(GroupA)#We need the s.d. Use the s.d. that is larger always!! In this case GroupA.

power.t.test(delta=2, sd=1.483, n=5)#This is the power analysis syntax for the earlier t-test.

#The above power analysis gives a test power of 46.6%. Not good! This means that only in 46% of cases will you find the effect if there is one!
#Usually run a power analysis if an effect has not been found to see if the experiment was could enough to find it.

power.t.test(delta=2, sd=1.483, power=0.95)#This test gives us the required number of n required for a 95% test power experiment.
#The above test gives an n=15.3 showing that at least 16 replicates are required for a robust test, prob more.

###Non-parametric Wilcox ranked t-test###

wilcox.test()