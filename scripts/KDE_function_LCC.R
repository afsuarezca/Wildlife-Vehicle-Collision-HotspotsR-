#set working directory
setwd(path_master)

#Road network
roadshp<-st_read(paste0(KDEdata_directoryName,"/",roads))

# Ensure the geometry type is LINESTRING
if (!all(st_geometry_type(roadshp) == "LINESTRING")) {
  print("All geometries must be LINESTRING.")
  #convert to LINESTRING#
  print("Converting LINESTRING.")
  roadshp<-st_cast(roadshp, "LINESTRING")
}

#Roadkills
roadkills<-st_read(paste0(KDEdata_directoryName,"/",roadkills))

if(plotroad == "yes"){
  # Plot the data####
  tm_shape(roadshp) + 
    tm_lines("black") + 
    tm_shape(roadkills) + 
    tm_dots("red", size = 0.2)
}

#Split roads to segments of 300 m######

# Function to split a line into segments of a specified length
split_line <- function(line, max_length) {
  # Calculate the number of points needed along the line
  total_length <- st_length(line)
  num_points <- ceiling(as.numeric(total_length) / max_length)
  
  # Generate points along the line
  points <- st_line_sample(line, n = num_points, type = "regular")
  
  # Split the line at the generated points
  segments <- st_split(line, points)
  
  # Return the resulting segments
  return(st_collection_extract(segments, "LINESTRING"))
}

kde_LCC<-function(roadshp,roadkills,roadl,mindist,
                  bw,method,kernel){
  # Calculate lixels to use as sampling points#####
  #Cut the lines of a feature collection of linestrings into lixels with a specified minimal distance 
  lixels <- lixelize_lines(roadshp,#road network
                           roadl,#lixel length
                           mindist = mindist#Minimum length of a lixel
  )
  
  samples <- lines_center(lixels)
  
  densities <- nkde(roadshp, 
                    events = roadkills,
                    w = rep(1,nrow(roadkills)),
                    samples = samples,
                    kernel_name = kernel,
                    bw = bw, div= "bw", 
                    method = method, digits = 1, tol = 1,
                    grid_shape = c(1,1), max_depth = 8,
                    agg = 5, #we aggregate events within a 5m radius (faster calculation)
                    sparse = TRUE,
                    verbose = FALSE)
  samples$density <- densities
  
  #Map the density values estimated for each lixel centre:####
  
  # rescaling to help the mapping
  samples$density <- samples$density*1000
  samples2 <- samples[order(samples$density),]
  
  return(samples2)
  
}

samples2<- kde_LCC(roadshp,
                   roadkills, 
                   roadl,
                   mindist,
                   bw ,
                   method,
                   kernel)

#Generate plot#######

#Define color palette for plotting
colorRamp <- brewer.pal(n = 7, name = "Spectral")
colorRamp <- rev(colorRamp)

#Plot final results
tm_shape(roadshp) + 
  tm_lines("black") + 
  tm_shape(samples2) + 
  tm_dots("density", style = "kmeans", palette = colorRamp, n = 7, size = 0.1) + 
  tm_layout(legend.outside = TRUE, 
            main.title = paste0("vehicle-wildlife collisions density by kilometres in LCC,",
                                "\nwithin a radius of 300 metres"),
            main.title.size = 1)


