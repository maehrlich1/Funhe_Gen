######PLOTTING MAPS#########

#######STANDARD MAPS#######

library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency

pcontorta <- readShapePoly("pinucont.shp")   #layer of data for species range
samps <- read.csv("FieldSamples.csv")   #my data for sampling sites, contains a column of "lat" and a column of "lon" with GPS points in decimal degrees
map("worldHires","Canada", xlim=c(-140,-110),ylim=c(48,64), col="gray90", fill=TRUE)  #plot the region of Canada I want
map("worldHires","usa", xlim=c(-140,-110),ylim=c(48,64), col="gray95", fill=TRUE, add=TRUE)  #add the adjacent parts of the US; can't forget my homeland
plot(pcontorta, add=TRUE, xlim=c(-140,-110),ylim=c(48,64), col=alpha("darkgreen", 0.6), border=FALSE)  #plot the species range
points(samps$lon, samps$lat, pch=19, col="red", cex=0.5)  #plot my sample sites


####COASTLINES######

library(rworldxtra)

worldmap<-getMap(resolution="high") #Imports the world map from the package in high resolution

plot(worldmap) #Check the loaded map

atlantic <- plot(worldmap, xlim = c(10,10.2), ylim = c(53.9,54.8), col = "gray") #Coordinates of map window required

#####NOAA QUALITY MAPS#########

####BATHYMETRIC MAPS#######

library(marmap)

KBbathy <- getNOAA.bathy(lon1 = 9.975, lon2 = 10.25, lat1 = 54.3, lat2 = 54.6, resolution = 1) #Queries the NOAA database for the bathymetric map. A matrix is created and cropped to specified size. Resolution is given in minutes.

summary(KBbathy) #Summary data including min and max depths

plot(KBbathy, image = "TRUE") #The image command adds a bluescale to the depths

scaleBathy(KBbathy, deg = 0.05, x = "bottomleft") #Adds a scale to the map. The deg command sets the length of the scale in lon degrees.

blues <- colorRampPalette( c( "purple", "darkblue", "blue", "cadetblue1","white")) #Creates a colour palette that ranges from purple via blue to white.

plot(KBbathy, image = "TRUE", bpal = blues(100), n = 5) #The bpal command allows the use of a custom colour scale. The argument gives the resolution of the colour scale. The n command gives the number of isobaths.

plot(KBbathy, image = "TRUE", bpal = blues(100), n = 5,
     col = c("black","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey"),
     lwd = c(1.5,1,1,1,1,1))

#col command dictates the colours of isobaths, starting with the coastline. The lwd gives the linewidths in the same way.

beigepal <- colorRampPalette("beige") #Create another colour palette

plot(KBbathy, n = 10, image = TRUE, land = TRUE, 
     bpal = list(c(0,max(KBbathy),"beige"),c(min(KBbathy),0,blues(100))), 
     col = c("black","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey"), 
     lwd = c(1.5,1,1,1,1,1,1,1,1,1,1,1,1))

#The list argument in the bpal allows for the differentiated colouration of different regimes e.g. land.

plot(KBbathy, n=1, lwd= 1.5, add = TRUE) #The add command overlays a new map on the old! Useful for coastlines!

points(Stations$Longitude, Stations$Latitude, pch = 21 , col = "black", bg = "red") #Points can be added to the map and taken from a dataset.
