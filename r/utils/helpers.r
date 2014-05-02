us.shp <- readShapeLines('data/us_alb.shp',
                         proj4string=CRS('+init=epsg:3175'))

ll2albers <- function(lat, long){
  
  centers_ll = data.frame(cbind(long, lat))
  colnames(centers_ll) = c('x', 'y')

  coordinates(centers_ll) <- ~ x + y
  proj4string(centers_ll) <- CRS('+proj=longlat +ellps=WGS84')

  centers_alb <- spTransform(centers_ll, CRS('+init=epsg:3175'))
  centers_alb <- as.matrix(data.frame(centers_alb))
  
  return(centers_alb)
}

plot.map <- function(x, coords_albers, map=us.shp, type){
 
  if (type=='neo'){
    site = paste(x$metadata$dataset$collection.handle, '_', x$metadata$dataset$dataset.id, sep='')
  } else if (type=='clh'){
    site = paste(toupper(x$site[1]), '_', x$datasetID[1], sep='')
  }
  
  path = 'maps'#file.path(mainDir, 'maps1')
  
  if (!file.exists(path)){
    dir.create(path)
  }
  
  while(length(dev.list()) > 1)dev.off()
  try(dev.off())
  
  pdf(file=paste(path, '/', site, '_MAP.pdf', sep=''), 
      width=6, height = 4)
  plot(coords_albers[1], coords_albers[2], type='n', xlim = c(5000, 1001000), ylim=c(600000, 1430000), 
       asp=1, xaxt='n', yaxt='n', ann=FALSE)
  plot(us.shp, add=T, lwd=2)
  points(coords_albers[1], coords_albers[2],col='orange', pch=19, cex=1.5)
  title(main=site)
  dev.off()
  
}


plot.pd <- function(x, coords, map=us.shp){
  
  site = paste(x$metadata$dataset$collection.handle, '_', x$metadata$dataset$dataset.id, sep='')
#   path = file.path(mainDir, site)
  
  path_pd =  'pollen_diagrams'#file.path(mainDir, 'pollen_diagrams')
  if (!file.exists(path_pd)){
    dir.create(path_pd)
  }
  
  path_count = 'counts'#file.path(mainDir, 'counts')
  if (!file.exists(path_count)){
    dir.create(path_count)
  }
  
  no_age = any(is.na(x$sample.meta$Age))
  
  if (no_age){
    age = x$chronologies[1][[1]]$Age
  } else{
    age = x$sample.meta$Age
  }
  
  good.samps <- which(age < 2000)
  
  if (length(good.samps) == 0){
    print(paste('WARNING: ', site, ' has no samples younger than 2000 YBP!', sep=''))
    return(NULL)
  }
  
  pol <- x$counts[good.samps,]
  pol.counts <- compile_list(pol, 'WhitmoreSmall')
  pol <- pol.counts / rowSums(pol.counts)
  pol <- pol[ ,!colnames(pol) %in% 'Other']
  
  markers = c('AMBROSIA', 'RUMEOXYR', 'POACEAE')
  
  most.common <- colnames(pol)[rank(colMeans(pol)) > (ncol(pol) - 10)]
  most.common <- unique(c(most.common, markers[markers  %in% colnames(pol)]))
  
  while(length(dev.list()) > 1)dev.off()
  
  try(dev.off())
  
  core.dat <- data.frame(depths = x$sample.meta$depths[good.samps], pol.counts[,most.common])
  write.table(core.dat, file=paste(path_count, '/', site, '.csv', sep=''),
              row.names=FALSE, sep=',')
  
  id.mod5   <- seq(1,length(good.samps))%%5
  sample.id <- rep(0, length(good.samps))
  sample.id[which(id.mod5 == 0)] = 1
  sample.id <- as.factor(sample.id)
  
  props.dat <- data.frame(depths = x$sample.meta$depths[good.samps],  sample.id=sample.id, pol[,most.common] * 100)
  
  m2 <- props.dat
  m2[,3:ncol(m2)] <- m2[,3:ncol(m2)] * 5
  
  m.melt  <- melt(props.dat, id = c('depths', 'sample.id'))
  m2.melt <- melt(m2, id = c('depths', 'sample.id'))
  
  m <- data.frame(rbind(m.melt, m2.melt),
                  exagr = factor(rep(c('no', 'yes'), each = nrow(m.melt)), 
                                 levels=c('no', 'yes')))
  
  
  levels(m$variable)[levels(m$variable) == 'RUMEOXYR'] <- 'RUMEX'
  
  
  m$value[m$value > 15 & m$exagr == 'yes'] <- 15
  m$value[m$value > 15 & m$exagr == 'yes'] <- NA
  
  taxa   = unique(m$variable)
  for (i in 1:length(taxa)){
    taxon = taxa[i]
    if (any(m$value[m$variable == taxon]>15,na.rm=TRUE)){
      m$value[m$exagr == 'yes' & m$variable == taxon] <- m$value[m$exagr == 'no' & m$variable== taxon]
    } 
  }
  
  hline_top = min(m$depths) - (max(m$depths) - min(m$depths))/100
  hline_bottom = max(m$depths) + (max(m$depths) - min(m$depths))/100
  
  p <- ggplot(m, aes(y=depths, x=value))+
    geom_path(colour='steelblue4', aes(x=value, y=depths, linetype = exagr)) +
    geom_vline(x=15, alpha=0) + geom_hline(aes(yintercept=hline, alpha=0), data=data.frame(hline=hline_top)) +
    geom_hline(aes(yintercept=hline, alpha=0), data=data.frame(hline=hline_bottom)) +
    geom_segment(m[m$exagr == 'no',], mapping=aes(y=depths, x=0, yend=depths, xend=value, colour=sample.id)) +
    scale_color_manual(values=c('steelblue4', 'darkorange1')) + 
    facet_grid(.~variable, space='free_x', scale='free_x') + 
    theme_bw(base_size=10) + ggtitle(site) +
    theme(legend.position='none',axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=rel(1)), panel.margin = unit(0.5, "lines")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_reverse(expand=c(0,0)) + xlab('percent')
#    print(p)
  ggsave(file=paste(path_pd, '/', site, '.pdf', sep=''), width=12, height=8)
}

plot.pd.clh <- function(x, coords_albers, map=us.shp){
  
  site = paste(toupper(x$site[1]), '_', x$datasetID[1], sep='')
  
  path_pd = 'pollen_diagrams'#file.path(mainDir, 'pollen_diagrams')
  if (!file.exists(path_pd)){
    dir.create(path_pd)
  }
  
  path_count = 'counts'#file.path(mainDir, 'counts')
  if (!file.exists(path_count)){
    dir.create(path_count)
  }
  
  no_age = any(is.na(x$age))
  
  good.samps <- which(x$age < 2000)
  
  if (length(good.samps) == 0){
    print(paste('WARNING: ', site, ' has no samples younger than 2000 YBP!', sep=''))
    return(NULL)
  }
  
  counts = x[,12:ncol(x)]
  
  pol <- counts[good.samps,]
  pol.counts <- compile_list(pol, 'WhitmoreSmall')
  pol  <- pol.counts / rowSums(pol.counts)
  
  pol <- pol[ ,!colnames(pol) %in% 'Other']
  
  most.common <- colnames(pol)[rank(colMeans(pol)) > (ncol(pol) - 10)]
  
  markers = c('AMBROSIA', 'RUMEOXYR', 'POACEAE')
  
  most.common <- unique(c(most.common, markers[markers  %in% colnames(pol)]))
  
  if (any(colSums(pol.counts[,most.common]) == 0)){
    most.common = most.common[which(colSums(pol.counts[,most.common]) != 0)]
  }
  
  while(length(dev.list()) > 1)dev.off()
  
  try(dev.off())
  
  core.dat <- data.frame(depths = x$depth_mid[good.samps], pol.counts[,most.common])
  write.table(core.dat, file=paste(path_count, '/', site, '.csv', sep=''),
              row.names=FALSE, sep=',')
  id.mod5   <- seq(1,length(good.samps))%%5
  sample.id <- rep(0, length(good.samps))
  sample.id[which(id.mod5 == 0)] = 1
  sample.id <- as.factor(sample.id)
  
  props.dat <- data.frame(depths = x$depth_mid[good.samps],  sample.id=sample.id, pol[,most.common] * 100)

  m2 <- props.dat
  m2[,3:ncol(m2)] <- m2[,3:ncol(m2)] * 5
  
  m.melt  <- melt(props.dat, id = c('depths', 'sample.id'))
  m2.melt <- melt(m2, id = c('depths', 'sample.id'))
  
  m <- data.frame(rbind(m.melt, m2.melt),
                  exagr = factor(rep(c('no', 'yes'), each = nrow(m.melt)), 
                                 levels=c('no', 'yes')))
  
  
  levels(m$variable)[levels(m$variable) == 'RUMEOXYR'] <- 'RUMEX'
  
  
  m$value[m$value > 15 & m$exagr == 'yes'] <- 15
  m$value[m$value > 15 & m$exagr == 'yes'] <- NA
  
  taxa   = unique(m$variable)
  for (i in 1:length(taxa)){
    taxon = taxa[i]
    if (any(m$value[m$variable == taxon]>15,na.rm=TRUE)){
      m$value[m$exagr == 'yes' & m$variable == taxon] <- m$value[m$exagr == 'no' & m$variable== taxon]
    } 
  }
  
  hline_top = min(m$depths) - (max(m$depths) - min(m$depths))/100
  hline_bottom = max(m$depths) + (max(m$depths) - min(m$depths))/100
  
  p <- ggplot(m, aes(y=depths, x=value))+
    geom_path(colour='steelblue4', aes(x=value, y=depths, linetype = exagr)) +
    geom_vline(x=15, alpha=0) + geom_hline(aes(yintercept=hline, alpha=0), data=data.frame(hline=hline_top)) +
    geom_hline(aes(yintercept=hline, alpha=0), data=data.frame(hline=hline_bottom)) +
    geom_segment(m[m$exagr == 'no',], mapping=aes(y=depths, x=0, yend=depths, xend=value, colour=sample.id)) +
    scale_color_manual(values=c('steelblue4', 'darkorange1')) + 
    facet_grid(.~variable, space='free_x', scale='free_x') + 
    theme_bw(base_size=10) + ggtitle(site) +
    theme(legend.position='none',axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=rel(1)), panel.margin = unit(0.5, "lines")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_reverse(expand=c(0,0)) + xlab('percent') + ylab('depth')
  
  ggsave(file=paste(path_pd, '/', site, '.pdf', sep=''), width=12, height=8)

}