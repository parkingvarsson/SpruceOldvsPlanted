library("devtools")
library("data.table")
library("rgdal")
library("rgeos")
library("rworldmap")
library("RColorBrewer")
library("scales")
library("maptools")
library("fossil")
library("vegan")
library("envirem")
library("geosphere")

##Functions
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

##Read WorldClim data
#midH<-stack(list.files(path="~/Dropbox/RTools/envirem/Europe_holo_ccsm4_2.5arcmin",pattern=".tif",full.names=T))
##Envirem data
current<-stack(list.files(path="~/Dropbox/RTools/envirem/Eurasia_current_2.5arcmin",pattern=".tif",full.names=T))
rcp45<-stack(list.files(path="~/Dropbox/RTools/envirem/Scandinavia_ccsm4_rcp45_2070_2.5arcmin",pattern=".tif",full.names=T))
rcp85<-stack(list.files(path="~/Dropbox/RTools/envirem/Scandinavia_ccsm4_rcp85_2070_2.5arcmin",pattern=".tif",full.names=T))

##WorldClim data
#current<-stack(list.files(path="~/Dropbox/RTools/worldclim/current_2.5arcmin",pattern=".bil",full.names=T))


##Populations

##Current climate
spruce.coords.sp<-SpatialPoints(coords[,c("Long","Lat")],proj4string = current@crs)
currspruce.env<-extract(current,spruce.coords.sp,method='bilinear',df=T)
rcp45.env<-extract(rcp45,spruce.coords.sp,method='bilinear',df=T)
rcp85.env<-extract(rcp85,spruce.coords.sp,method='bilinear',df=T)


samples<-cbind(coords,currspruce.env[,-1])
samples4.5<-cbind(coords,rcp45.env[,-1])
samples8.5<-cbind(coords,rcp85.env[,-1])

##Create raster data Spruce native range
##Europe
e<-extent(c(-15,60,40,75))
##Sweden
s<-extent(c(6,25,60,69))


data(wrld_simpl)

##P. abies distribution
pabies.native <- readOGR("~/Dropbox/Spruce/spruce_3.0/distribution_data/Picea_abies/shapefiles/Picea_abies_plg_clip.shp")
pabies.nonnative <- readOGR("~/Dropbox/Spruce/spruce_3.0/distribution_data/Picea_abies/shapefiles/Picea_abies_syn_plg_clip.shp")


##Extract current climate data
spruce.current<-crop(current,e)
spruce.current<-crop(current,s)
spruce.rcp45<-crop(rcp45,s)

#Create raster data for Sweden
#Using shape file
#Sweden shapefile
e<-extent(c(5.5,30,54,71))
data(wrld_simpl)
SWE<-subset(wrld_simpl, NAME %in% c("Sweden"))

r.spruce<-crop(current,extent(SWE))
r.spruce<-mask(r.spruce,SWE)
r.spruce<-crop(r.spruce,extent(pabies.native))
r.spruce<-mask(r.spruce,pabies.native)
names(r.spruce)<-substr(names(r.spruce),19,50)
r.spruce.table<-as.data.frame(rasterToPoints(r.spruce))











##Plot current climate data
col1<-addalpha("indianred3",alpha=0.9)
col2<-addalpha("indianred1",alpha=0.9)
col3<-addalpha(rgb(1,0.85,0.85),alpha=0.9)
col4<-addalpha(rgb(0.85,0.95,1.0),alpha=0.9)
col5<-addalpha("lightsteelblue2",alpha=0.9)
col6<-addalpha("lightsteelblue4",alpha=0.9)
col.palette<-rev((colorRampPalette(c(col1,col2,col3,col4,col5,col6),alpha=TRUE)(100)))
col.palette<-((colorRampPalette(c("snow1","snow2","snow3","seagreen","orange","firebrick"),alpha=TRUE)(100)))
col7<-addalpha("grey85",alpha=0.15)
col8<-addalpha("grey75",alpha=0.15)

plot(spruce.current[[4]],col=col.palette,las=1)
plot(spruce.rcp45[[4]],col=col.palette,las=1)
par(mfrow=c(1,2))
plot(spruce.current[[1]],col=col.palette,las=1,legend=T,xaxt="n",yaxt="n",bty="n",xlab="Longitude",ylab="Latitude",xlim=c(-15,60))
plot(pabies.native,add=T,col=col5)
plot(pabies.nonnative,add=T,col=col6)
plot(spruce.coords.sp,add=T,pch=21,cex=0.3,bg=col.vector)
xplot(spruce.coords.sp,add=T,pch=21,cex=1.2,bg=col.vector)
axis(1)
axis(2,las=1)

