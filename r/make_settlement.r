library(maptools)
library(sp)
library(rgdal)
library(neotoma)
library(reshape)
library(ggplot2)
library(grid)

source('r/utils/helpers.r')

use_date = '2014-05-01'
mainDir  = getwd()

all.sites <- read.csv(paste('pollen_meta_', use_date, '.csv', sep=''), stringsAsFactors=FALSE)

clh.sites <- all.sites[substr(all.sites$datasetID,1,3) == 'CLH',]

neo.sites <- all.sites[substr(all.sites$datasetID,1,3) != 'CLH',]

centers_alb = ll2albers(neo.sites$lat, neo.sites$long)

neo.sites$x = centers_alb[,1]
neo.sites$y = centers_alb[,2]

for (i in 1:nrow(neo.sites)){
  
  print(i)
  
  # get the data from neotoma
  if (neo.sites$datasetID[i] == 2290){
    print('Rossburg')
    next
  }
  x = get_download(as.numeric(neo.sites$datasetID[i]))
  
  coords = c(centers_alb[i,1], centers_alb[i,2])
  
  # make the pollen diagram and
  # put the site on a map
  plot.pd(x, coords, map=us.shp)
  plot.map(x, coords, map=us.shp, type='neo')
}

clh <- read.csv('data/pollen_counts_calcote.csv', stringsAsFactors=FALSE)

centers_alb = ll2albers(clh.sites$lat, clh.sites$long)

clh.sites$x = centers_alb[,1]
clh.sites$y = centers_alb[,2]

for (i in 1:nrow(clh.sites)){
  print(i)
  
  idx = which(clh$site == clh.sites$site[i])
  
  x = clh[idx,1:ncol(clh)]
  coords = c(centers_alb[i,1], centers_alb[i,2])
  plot.pd.clh(x, coords, map=us.shp)
  plot.map(x, coords, map=us.shp, type='clh')
  
}


# simon thinks one big pdf will save work!
system("gs -sDEVICE=pdfwrite -o pollen_diagrams.pdf pollen_diagrams/*.pdf")

system("gs -sDEVICE=pdfwrite -o maps.pdf maps/*.pdf")

# why is my ggmap totally busted?
# myLocation <- c(lon = -95.3632715, lat = 29.7632836) 
# myMap <- get_map(location=myLocation, 
#                  source="stamen", maptype="watercolor", crop=FALSE) 
# ggmap(myMap) 