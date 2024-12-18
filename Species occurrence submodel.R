# Author: Andres Felipe Suarez-Castro
# Date: November 2024

-------------------------------------------------------------------------------------------------
  # Title: Species occurrence models
  # Component: Exposure model
  # Spatial occupancy models for each species using survey records from the Atlas of Living Australia (ALA)
  -------------------------------------------------------------------------------------------------

# Set working directory, clean workspace and load packages#######
rm(list = ls())

# Set working directory
setwd("//uq.edu.au/UQ-Research/LOGANROAD-A1321")

library(biomod2)
library(randomForest)
library(terra)
library(caret)
library(spatial.tools)
library(tidyselect)
library(dismo)
library(reshape)
library(sf)

#Load data#########

#load ALA data
data_sp<-read.csv("./Data_species/Macropus_rufogriseus/Macropus_rufogriseus.csv")#
#data_sp<-read.csv("sp_coordinate_presence.csv")#for presence
#data_sp<-read.csv("Wallabi_records.csv")

#remove data with missing coordinates
sp_data<-subset(data_sp, !is.na(data_sp$Longitude) & data_sp[35] == "Queensland")
# which records are duplicates?
dups <- duplicated(sp_data)
# remove duplicates
sp_data <- sp_data[!dups, ]


# Save graphics defaults
par.defaults <- par(no.readonly=TRUE)
save(par.defaults, file="R.default.par.RData")

# Tell biomod2 which parts of the database refer to which biomod2 object
myRespName <- "Macropus_rufogriseus"
myRespXY <- sp_data[,c("Longitude","Latitude")]
coordinates.spes<- sp_data[,c("Longitude","Latitude")]
#transform to appropiate coordinate system
coordinates(coordinates.spes) <- c("Longitude","Latitude")
proj4string(coordinates.spes) <- CRS("+init=epsg:4326") # WGS 84
CRS.new <- CRS("+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
xy<-spTransform(coordinates.spes, CRS.new)
SPDF <-SpatialPointsDataFrame(coords=xy, data=myRespXY)


#Load environmental covariates#########

myExpl<-stack("MyExpl_LCC2.grd")

# Load in environmental covariates

FPC<- raster("./Data_layers/FPC_SEQ")#Foliage projective cover
#system.time(FPC<-projectRaster(FPC,res = c(30,30), crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
DEM<-raster("./Data_layers/dem_seq_30")#Digital elevation model
#system.time(DEM<-projectRaster(DEM,res = c(30,30), crs = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
map.r <- raster("./Data_layers/map_30_seq")# Mean annual precipitation/rainfall; 
mat.r <- raster("./Data_layers/mat_30_seq2")# Mean annual temperature
bvg<-raster("./Data_layers/BVG/bgv17_seqf")
slope<-raster("./Data_layers/slope_seq30")


#reading built up areas from geodatabase
#require(rgdal)
#fgdb <- "./Data_layers/built_up/built_up/data.gdb"
#subset(ogrDrivers(), grepl("GDB", name))
#fc_list <- ogrListLayers(fgdb)
#print(fc_list)
# Read the feature class
#built_up <- readOGR(dsn=fgdb,layer="Built_up_areas")

-------------------------------------------------------------------
#it is necessary to convert the layers that are factors to numeric
#  to avoid Error in sample.int(length(x), size, replace, prob) : 
# invalid first arguments in BIOMOD_FormatingData
-------------------------------------------------------------------

is.factor(map.r)
#[1] TRUE
map.r<-deratify(map.r,"ID")
FPC<-deratify(FPC,"ID")
bvg<-deratify(bvg, "ID")

# Create stack of all environmental covariate layers
myExpl <- stack(list(precipitation=map.r,temperature=mat.r,FPC=FPC,
                     elevation=DEM, bvg=bvg, slope=slope))#,distance=dist.sync

writeRaster(myExpl, file="MyExpl_LCC2.grd", format="raster",overwrite = TRUE)

myExpl<-subset(myExpl,c(3,4,6,7))#exclude temperaturE and precipitation

#read shapefile Logan
LCC<-readOGR(dsn = "./Data_layers", layer = "Logan_lga_GDA")
#read shapefile SEQ
SEQ<-readOGR(dsn = "./Data_layers", layer = "SEQ_project")

#check layers extent and stack predictors
extent_Logan<-extent(LCC)
extent_SEQ<-extent(SEQ)

myExplLCC<-crop(myExpl,extent_Logan)#this generates a rasterbrick (does not work in biomod)
myExplLCC<-stack(myExplLCC)#convert again to stack

#myExpl<-subset(myExpl,c(3:6))#exclude temperaturE and precipitation

myExpl$bgv17_seqf<-as.factor(myExpl$bgv17_seqf)

# replacing NA's by zero
myExpl$FPC[is.na(myExpl$FPC[])] <- 100

#save(myExpl,file = "myExpl.RData")# if you saved a large raster like this, the values of large Raster* objects are written to file, 
#to avoid memory limitation problems. If you do not explicitly provide a filename, they are stored in the temporary data folder 
#that will be removed when the R session ends.o avoid that from happening, you can  write them to a permanent file with writeRaster

#Extract records that are inside the study area and sample points based on distance#####

DF <- as.data.frame(SPDF)
DF$values<-raster::extract(myExpl$bgv17_seqf,SPDF)#extract coordinates presence points

DF<-subset(DF,!is.na(DF$values))

DF<-DF[sample(nrow(DF),500),]  

#require(spatialEco)
#require(spatstat)
#require(sp) 

#n = round(length(SPDF) * 0.1, digits=0) # select five percent data 
#SPDF.wrs <- pp.subsample(SPDF, n=n, window='extent') 

#SPDF.wrs  <- SPDF[sample(1:length(SPDF),n),]

# And then convert it (back) to a data.frame  

#DF<-DF[sample(nrow(DF),100),]
myRespXY <- DF[,c("Longitude.1","Latitude.1")]

#define response (presence)

DF$Presence<- 1

myResp <- as.numeric(DF[,"Presence"])

########Create the ensemble models using biomod2########


# Create object (myBiomodData) to contain all the previous objects within it, 
# formatted correctly

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.absences = 100,
                                     PA.nb.rep = 3,
                                     PA.strategy = "random")#,
                                     #PA.sre.quant = 0.025
                                     #PA.dist.min = 350)#when selecting disk

#Error in .check.params.sre(Response, Explanatory, NewData, Quant) : 
#SRE algorithm does not handle factorial variables

# Set the options that you have chosen for the different model algorithms
# - here using defaults (empty arguments) 
myBiomodOption <- BIOMOD_ModelingOptions()


# Run the biomod2 models that you have chosen on the data provided.
myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('RF','GBM','GAM'),
  models.options = myBiomodOption,
  NbRunEval = 15,
  DataSplit = 70,
  Prevalence = 0.5,
  VarImport = 3,
  models.eval.meth = c('ROC',"TSS"),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName, "Species1", sep=""))

save(myBiomodModelOut, file = "./Occurrence_models/Macropus_rufogriseus/PC_presence_4mar.RData")

#load("run1_presencesp.RData")

# Get evaluation statistics for each model####
myBiomodModelEval <- get_evaluations(myBiomodModelOut)

save(myBiomodModelEval, file = "./Occurrence_models/Macropus_rufogriseus_myBiomodModelEval_4mar.RData")

#load("myBiomodModelEval1.RData")

# Display the model evaluation statistics
TSS.results<-myBiomodModelEval["TSS", "Testing.data" ,,,]

TSS.results<-melt(TSS.results, id=c("GLM","GBM","RF"))
TSS.results.top<-TSS.results[order(-TSS.results["value"]),][1:5,]
TSS.results.top<-apply(TSS.results.top,1,function(x) paste(c(x[1],x[2]),collapse="_"))


ROC.results<-myBiomodModelEval["ROC", "Testing.data",,,]
ROC.results<-melt(ROC.results, id=c("GLM","GAM","RF"))
ROC.results.top<-ROC.results[order(-ROC.results["value"]),][1:3,]
ROC.results.top<-apply(ROC.results.top,1,function(x) paste(c("Macropus.rufogriseus",x[3],x[2],x[1]),collapse="_"))

get_variables_importance(myBiomodModelOut)

# Inspect model response curves for implausible responses

#myGLMModels <- BIOMOD_LoadModels(myBiomodModelOut, models=c('GLM'))

myRespPlotGLM <- response.plot2(models  = myGLMModels,
                                Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                                show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                                do.bivariate = FALSE,
                                fixed.var.metric = 'median',
                                col = c("red", "blue", "green"),
                                legend = TRUE,
                                data_species = get_formal_data(myBiomodModelOut,'resp.var'))


#myRFModels <- BIOMOD_LoadModels(myBiomodModelOut, models=c('RF'))

myRespPlotRF <- response.plot2(models  = c(ROC.results.top),
                               Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                               show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                               do.bivariate = FALSE,
                               fixed.var.metric = 'median',
                               #col = c("red", "blue", "green"),
                               legend = TRUE,
                               data_species = get_formal_data(myBiomodModelOut,'resp.var'))

# Project selected models for later building ensemble model######
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplLCC,
  proj.name = 'current',
  selected.models = c(ROC.results.top),
  binary.meth = 'ROC',
  compress = 'xz',
  build.clamping.mask = FALSE,
  output.format = '.grd')

# Plot selected models, if you want
#plot(myBiomodProj, str.grep = "RF")

save(myBiomodProj, file = "./Occurrence_models/Macropus_rufogriseus/PC_myBiomodProj_4mar.RData")

# Build ensemble model from selected models########
myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                       chosen.models = c(ROC.results.top),
                                       em.by = 'all',
                                       eval.metric = c('ROC'),
                                       eval.metric.quality.threshold = c(0.6),
                                       prob.mean = TRUE,
                                       prob.cv = FALSE,
                                       prob.ci = FALSE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = FALSE,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional')

save(myBiomodEM,file = "./Occurrence_models/Macropus_rufogriseus/PC_myBiomodEM_4mar.RData")

# Project and plot predictions from the ensemble model 
EMplot <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                     EM.output = myBiomodEM)


plot(EMplot)

# save map as a raster######
Ensemble_raster <- raster("Macropus.rufogriseus/proj_current/proj_current_Macropus.rufogriseus_ensemble.grd")

writeRaster(Ensemble_raster, file="Occurrence_models/Final_models/Macropus_rufogriseus_ensemble1_4mar", format="GTiff", overwrite = TRUE)
