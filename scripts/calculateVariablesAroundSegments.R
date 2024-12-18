if(generateTable == "yes"){
  #Species distribution model. Developed using code "Species occurrence submodel"
  
  spocras <- rast(paste0(EMdata_directoryName,"/",spoc))
  #plot(spocras)
  #Conver NA values of distribution model to 0
  spocras[is.na(spocras[])]<-0
  
  #Land cover categories as explained in table 1 of the final report
  vegetationtypes<- rast(paste0(EMdata_directoryName,"/",lccat))
  
  #Impervious surface
  impsur<- rast(paste0(EMdata_directoryName,"/",impervsur))
  
  #Vegetation raster
  vegetation<-rast(paste0(EMdata_directoryName,"/",veg))
  #plot(vegetation)
  
  #add all rasters to a list
  rasters<-list(spocras,vegetation,impsur,vegetationtypes)
  
  #--------------------#
  #Process rasters######
  #--------------------#
  
  # Check if rasters is already a list
  if (!is.list(rasters)) {
    # If not a list, create a list from the input
    rasters <- list(rasters)
  }
  
  # Ensure all elements in the list are SpatRaster objects
  if (!all(sapply(rasters, inherits, "SpatRaster"))) {
    stop("All inputs must be SpatRaster objects or a list of SpatRaster objects.")
  }
  
  # Check that all rasters are have the same extent and same crs as roads
  # Get the extent and CRS of the first raster
  base_extent <- ext(rasters[[1]])
  base_crs <- crs(roadshp)
  
  # Check each raster against the base extent and CRS
  
  checkcrs<-sapply(rasters, function(r) {
    crs(r) == base_crs 
  })
  
  if(all(checkcrs) == FALSE){
    print("Projecting all rasters to roads coordinate system")
    reference_raster<-project(spocras,crs(roadshp))
    rasters <- lapply(rasters, function(r) project(r, reference_raster))
  }
  
  #check that all rasters are in the same extent
  checkextent <- sapply(rasters, function(r) {
    identical(ext(r), base_extent) 
  })
  
  if(all(checkextent) == FALSE){
    print("Converting all rasters to the same extent and spatial resolution ")
    rasters <- lapply(rasters, function(r) extend(r, reference_raster))
  }
  
  
  # Combchecks# Combine rasters into a SpatRaster using `c()` and ensure alignment if needed
  allrast <- do.call(c, rasters)
  
  
  distweightedmean<-function(roadshp,rasters,buffer){
    
    # Generate buffers along roads
    
    roadbuff <- st_buffer(roadshp,buffer)
    #plot(roadbuff[1])
    
    #--------------------------------------------------------------------------------#
    #calculate weighted distance from the road to raster values within a buffer####
    #--------------------------------------------------------------------------------#
    
    pb <- txtProgressBar(min = 0, max = nrow(roadshp), style = 3)  # Style 3 for a dynamic bar
    
    m<-list()
    for(i in 1:nrow(roadshp)){
      setTxtProgressBar(pb, i)
      r_crop<-crop(rasters[[1]],roadbuff[i,])#first raster must always be the species distribution model
      r_mask<-mask(r_crop, roadbuff[i,])
      cells <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
      colnames(cells)<-c("x","y","sp")
      cells<-subset(cells,!is.na(cells$sp))
      points_sf <- st_as_sf(cells, coords = c("x", "y"), crs = crs(roadshp)) 
      distfromroad <- st_distance(points_sf, roadshp[1,], byid=TRUE)
      m[[i]]<-data.frame(cells,distfromroad=as.numeric(distfromroad),osm_id=roadshp[i,"osm_id"][[1]])
    }
    
    #now let's calculate the weighted mean of species probability occurrence for each site
    
    weighted.meandist<-function(value,distance,lambda)
    {
      value[is.na(value)]<-0
      denom<-value*(exp(-distance*lambda))
      numer<-exp(-lambda*distance)
      a<-sum(denom)
      b<-sum(numer)
      return(a/b)
    }
    
    #Calculate weighted mean for each road segment
    weighted.meanspecies<-list()
    
    for (i in 1:length(m)){
      wmean<-weighted.meandist(value=m[[i]]$sp,distance=m[[i]]$distfromroad,lambdasp)
      osm_id<-unique(m[[i]]$osm_id)
      weighted.meanspecies[[i]]<-data.frame(wmean,osm_id)
    }  
    
    #create a dataframe
    weighted.meanspecies<-as.data.frame(do.call(rbind,weighted.meanspecies))
    
    return(weighted.meanspecies)
  }
  
  #Calculate proportion of vegetation types around roads
  calcvar<-function(roadshp,r1,idcolumn){
    pb <- txtProgressBar(min = 0, max = nrow(roadbuff), style = 3) #To check the progress
    lcc<-list()
    for(i in 1:nrow(roadshp)){
      setTxtProgressBar(pb,i)
      r1crop<-crop(r1,roadbuff[i,])
      roadr1<-terra::rasterize(roadbuff[i,],r1crop)
      roadr1<-mask(r1crop,roadr1)
      #plot(roadlcc)
      cells<-data.frame(roadr1)
      #add road ID
      cells[,idcolumn]<-unique(roadbuff[[i,idcolumn]])
      
      colnames(cells)<-c("var","road")
      
      
      df_counts <- cells %>% 
        group_by(road)%>%
        dplyr::count(var) 
      
      colnames(df_counts)<-c(idcolumn,"var","n")
      
      #multiply by resolution
      df_counts$area<-df_counts$n * res(r1crop)[1]^2
      
      total<-sum(df_counts$area, na.rm = T)
      
      df_counts$prop<-df_counts$area/total
      
      df_counts$n<-NULL
      df_counts$area<-NULL
      
      
      #now spread de df
      
      df_counts<- spread(df_counts, key=var, value=prop)
      
      lcc[[i]]<-df_counts
    }
    return(do.call(plyr::rbind.fill,lcc))
  }
  
  #proportion of land cover types around roads
  veg100m<-calcvar(roadbuff,rasters[[4]],"osm_id")
  
  varriskmodel<-merge(veg100m,weighted.meanspecies,"osm_id")
  
  varriskmodel$imp.prop<-varriskmodel %>% sum(varriskmodel$`1-Impervious|1-Low`,
                                              varriskmodel$`1-Impervious|2-Med`,
                                              varriskmodel$`1-Impervious|3-High`,
                                              varriskmodel$`1-Impervious|4-very High`,
                                              na.rm = T)
  varriskmodel$barren.prop<-varriskmodel %>% sum(varriskmodel$`2-Barren|1-Low`,
                                              varriskmodel$`2-Barren|2-Med`,
                                              varriskmodel$`2-Barren|3-High`,
                                              varriskmodel$`2-Barren|4-very High`,
                                              na.rm = T)
  
  varriskmodel$veghigh<-varriskmodel$`3-Vegetated|3-High`
  varriskmodel$vegmid<-varriskmodel$`3-Vegetated|2-Med`
  varriskmodel$veglow<-varriskmodel$`3-Vegetated|1-Low`
  varriskmodel$vegveryhigh<-varriskmodel$`4-Vegetated|4-very High`
  
  varriskmodel<-varriskmodel[c("osm_id", "wmean","imp.prop","barren.prop",            
                               "veghigh", "vegmid","veglow")]
  
}else{
  varriskmodel<-read.csv(paste0(EMdata_directoryName,"/",myvar))
}