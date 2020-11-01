#NYC Squirrels
library(readr)
library(ggplot2)
library(viridis) #cool color palette to use
library(cowplot) #minimalistic theme for ggplot 
setwd("~/Desktop") #set working directory

#Import
nyc_squirrels <- readr::read_csv("https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2019/2019-10-29/nyc_squirrels.csv")

#Heatmap of squirrel density
ggplot(nyc_squirrels)+
  stat_summary_hex(aes(x=long, y=lat), bins=50)+                      #Hexagons avoid graphing artefacts due to boundary conditions
  coord_fixed()+                                                      #this scale the x and y axes in the same way, makes sense for spatial data
  scale_fill_viridis(option="inferno", name = "Squirrel Density")+    #apply the awesome color palette
  #facet_wrap(~shift)+                                                #create subplot using other variables! (or not...)
  theme_cowplot()                                                     #apply minimalistic theme to remove background

############################################################

vars <- c("shift", "age", "primary_fur_color", "running", "eating")   #create a vector containing all of the variables you would like to plot spatially

data[vars] <- lapply(data[vars], as.factor) #turn those variables from characters into factors

#Loop over all variables making a plot for each
for (z in vars){
  
  plot <- ggplot(nyc_squirrels)+
            geom_point(aes_string(x="long", y="lat", col=z), alpha=0.6)+      #the variable determines the color here. alpha makes points slightly transparent
            labs(x="Longitude", y="Latitude")+
            coord_fixed()+
            theme_cowplot()

  pdf(file=paste(z, ".pdf", sep=""), height=8, width=8)                       #starts writing to a pdf file
  plot(plot)
  dev.off()                                                                   #closes the pdf file
  
}
