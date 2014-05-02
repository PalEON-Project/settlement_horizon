library(maps)

get_neo <- function(){
 
  gpids <- get_table(table.name='GeoPoliticalUnits')

  gp_rows <- vector(length=3)
  gp_rows[1] <- which(gpids$GeoPoliticalName == 'Wisconsin')
  gp_rows[2] <- which(gpids$GeoPoliticalName == 'Minnesota')
  gp_rows[3] <- which(gpids$GeoPoliticalName == 'Michigan')
  
  ID1   <- gpids[gp_rows,1]
  
  
  ID1   <- gpids[which(gpids$GeoPoliticalName == 'Wisconsin'),1]
  ID2   <- gpids[which(gpids$GeoPoliticalName == 'Minnesota'),1]
  ID3   <- gpids[which(gpids$GeoPoliticalName == 'Michigan'),1]
  
  meta1  <- get_datasets(datasettype='pollen', gpid=ID1, ageold=(2000))
  meta2  <- get_datasets(datasettype='pollen', gpid=ID2, ageold=(2000))
  meta3  <- get_datasets(datasettype='pollen', gpid=ID3, ageold=(2000))
  
  meta   <- c(meta1, meta2, meta3)
  #meta   <- meta1
  
  # create of vector of the DatasetIDs 
  umw = c('michigan:north', 'wisconsin', 'minnesota', 'michigan:south')
  
  n      <- length(meta)
  ids    <- vector(mode="numeric", length=0)
  states <- vector(mode="numeric", length=0)
  lat    <- vector(mode="numeric", length=0)
  long   <- vector(mode="numeric", length=0)
  site   <- vector(mode="character", length=0)
  handle <- vector(mode="character", length=0)
  PI     <- vector(mode="character", length=0)
  descriptor <- vector(mode="character", length=0)
  
  for (i in 1:n){
    if(meta[[i]]$DatasetType == 'pollen'){
      
      x = unlist(meta[[i]]$Site$LongitudeWest)
      y = unlist(meta[[i]]$Site$LatitudeNorth)
      
      state = map.where(database = "state", x, y)
      print(state)
      
      if (is.element(state, umw)){
        print(i)
        
        states = c(states, state)
        ids    = c(ids, unlist(meta[[i]]$DatasetID))# unlist(meta[[i]]$DatasetID)
        long   = c(long, x)
        lat    = c(lat, y)
        site   = c(site, as.character(unlist(meta[[i]]$Site$SiteName)))
        handle = c(handle, unlist(meta[[i]]$CollUnitHandle))
        PI     = c(PI, as.character(unlist(meta[[i]]$DatasetPIs[1,2])))
        about = as.character(meta[[i]]$Site$SiteDescription)
        descriptor = c(descriptor, strsplit(about, '\\.')[[1]][1])
      }
    }
  }
  
  empty = rep(NA, length(ids))
  
  pollen <- data.frame(datasetID=ids, handle=handle, site=site,  long=long, lat=lat, state=states, 
                       pi=PI, pre=empty, settlement=empty, notes=empty,
                       stringsAsFactors=FALSE)#, description=descriptor)
  
  return(pollen)
}

get_calcote <- function(){
  
  clh.sites  <- read.csv('data/hotchkiss_lynch_calcote_meta.csv', stringsAsFactors = FALSE)
  clh.counts <- read.csv('data/hotchkiss_lynch_calcote_counts.csv', stringsAsFactors = FALSE)
  
  clh.sites$name <- gsub(" ","", clh.sites$name, fixed=TRUE)
  clh.counts$name <- gsub(" ","", clh.counts$name, fixed=TRUE)
  
  n <- nrow(clh.sites)
  
  ids    <- vector(mode="numeric", length=0)
  states <- vector(mode="numeric", length=0)
  lat    <- vector(mode="numeric", length=0)
  long   <- vector(mode="numeric", length=0)
  site   <- vector(mode="character", length=0)
  PI     <- vector(mode="character", length=0)
  descriptor <- vector(mode="character", length=0)
  
  site.count = 0
  
  for (i in 1:n){
    
    site.i <- as.character(clh.sites$name[i])
    idx    <- which(clh.counts$name == site.i)
    type   <- clh.counts[idx,1]
    
    if (length(type) > 2){
      
      site.count = site.count + 1
      
      ids    = c(ids, paste('CLH', site.count, sep=''))
      states = c(states, 'wisconsin')
      site   = c(site, site.i)
      lat    = c(lat, clh.sites$lat[i])
      long   = c(long, clh.sites$long[i])
      PI     = c(PI, clh.counts$analyst[idx[1]])
      
    }
    
  }
  
  empty = rep(NA, length(ids))
  
  pollen_meta <- data.frame(datasetID=ids, handle=toupper(site), site=site, long=long, lat=lat, state=states, 
                            pi=PI, pre=empty, settlement=empty, notes=empty, 
                            stringsAsFactors=FALSE)
  
  # combine the meta and counts into a single table
  pollen_all = NULL
  for (i in 1:nrow(pollen_meta)){
    idx <- which(clh.counts$name == pollen_meta$site[i])
    
    counts.i <- clh.counts[idx,!(colnames(clh.counts) %in% c('name', 'analyst'))]
    
    meta.i <- pollen_meta[rep(i, length(idx)),]
    
    pollen_all = rbind(pollen_all, data.frame(meta.i, counts.i, row.names=NULL))
  }
  
  write.table(pollen_all, 'data/pollen_counts_calcote.csv', row.names=FALSE, sep=',')
  
  return(pollen_meta)
  
}

pollen_meta_neo     <- get_datasetids()
pollen_meta_calcote <- get_calcote()

pollen_meta <- rbind(pollen_meta_neo, pollen_meta_calcote)
pollen_meta <- pollen_meta[with(pollen_meta, order(handle)),]

write.table(pollen_meta, paste('pollen_meta_', Sys.Date(), '.csv', sep=''), row.names=FALSE, sep=',', na='')
