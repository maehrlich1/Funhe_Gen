#####MANIPULATING DATA##########

###Creating a Data Set###

Lux<-c(rep("A",5),rep("B",5),rep("C",5),rep("D",5))    #Generates labels for treatment levels

GroupA<-rnorm(5,15,1)#Creates 5 random numbers with a mean of 15 and a s.d. of 1#
GroupB<-rnorm(5,12,1)
GroupC<-rnorm(5,8,1)
GroupD<-rnorm(5,6,1)

Depth<-c(GroupA,GroupB,GroupC,GroupD)#Collates the groups as the dependent variable vector

Data<-data.frame(Lux,Depth)#Generates a dataframe (matrix) combining both ind. and dep. vectors

Data$Depth<-as.integer(Data$Depth)#The $ sign identifies the vector you want to work with from the dataframe, turns numbers into integers.

################################

###Extracting data###

str(Data)#What is the structure of the dataframe? Is the Independent variable a factor?

length(Data$Depth)#Amount of obs in the vector

class(Data$Depth)#Type of data in the vector

levels(Data$Lux)#Type of data in the vector

table(Data$Lux)#Shows number of relicates per category. Can be used to check for a balanced design

max(Data$Depth)#Max value in the vector

min(Data$Depth)#Min value in the vector

tapply(Data$Depth, Data$Lux, min)#Finds the minimum response level (depth) across all factor levels (Lux)

tapply(Data$Depth, Data$Lux, max)#Finds the max response level (depth) across all factor levels (Lux)

tapply(Data$Depth, Data$Lux, mean)#Finds the mean response level (depth) across all factor levels (Lux)

tapply(Data$Depth, Data$Lux, hist)#Plots a histogram

all(Data$Depth>10)#Are all obs in Depth >10?

any(Data$Depth>10)#Are any of the obs in Depth >10?

sum(Data$Depth>10)#How many obs in Depth are larger than 10?

Data$Depth>10#Shows the data that is >10 (Boolean statement)?

########################
####Sorting####

Sorted0<-order(Data$Depth)#Sorts data in the vector in ascending order

Sorted1<-order(-Data$Depth)#Sorts the data in the vector in descendent order

Sorted3<-data[order(Data$Depth,)] #Sorts the entire dataframe according to a variable in a vector! See below

Data[Sorted0,]#Square brackets have [rows,columns] syntax, leaving it blank indicates "all". Forms subsets essentially,see below.

Data[Sorted1,]

#######################
#####Forming subsets#####

#Square brackets create data frames! Can be used for creating subsets

New1<-Data[Data$Depth>5,]#Take all elements of Depth that are larger than 5. DONT FORGET THE COMMA AT THE END!!!!

New2<-Data[Data$Depth>5 & Data$Depth<10,]#All depth data that is bigger than 5 and smaller than 10.

New3<-Data[Data$Depth==1 | Data$Depth==5 | Data$Depth==8,]#The | is the OR Operator. Takes values that are equal to 1,5 or 8.

New4<-Data[Data$Depth !=1 & Data$Depth !=5 & Data$Depth !=8,]#The ! excludes all data points that are 1,5 and 8.

New5<-Data[Data$Lux =="A" & Data$Depth >4,]#Onlytake elements that belong to Lux Group A (only work on one group). The second criterion determines the subgroup of A that I want.

######################
########Means & SD###############

my_means<-tapply(FRSS$NoFG, FRSS$Station, mean) #Gives the mean for every group
my_sd<-tapply(FRSS$NoFG, FRSS$Station, sd) #Gives the sd for every group

mean_sd<-data.frame(Mean = my_means, SD = my_sd) #Creates dataframe for means and associated sd
