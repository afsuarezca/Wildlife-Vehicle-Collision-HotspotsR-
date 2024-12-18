#LIBRARIES
library("foreign")
library("runjags")
library("coda")
library("rjags")
library("snowfall")
library("parallel")


#SET UP DATA FOR JAGS

#get data into format for analysis

Data<- varriskmodel

#select relevant variables
Data <- Data[c("Road_Block","veg500", "veg100m", "sinuosity", "Road_AADT", "wmean.kangaroo","wmean.koala","wmean.wallaby","Road_Speed","Road_Width","Length", "Kangaroo_counts",
               "Koala_counts","Wallabies_counts")]

#convert road speed to numeric and Length to km
Data$Road_Speed <- as.numeric(Data$Road_Speed)
Data$Length <- Data$Length / 1000

#standardize data except counts and length
scale.data <- as.data.frame(scale(Data[c("veg500", "veg100m", "sinuosity", "Road_AADT", "wmean.kangaroo","wmean.koala","wmean.wallaby", "Road_Speed","Road_Width")]))

#bind count data and standardized variables
Data<-cbind(Data$Kangaroo_counts,Data$Koala_counts,Data$Wallabies_counts, scale.data, Data$Length)

#assign names
names(Data) <- c("Kangaroo_counts", "Koala_counts","Wallabies_counts", "veg500", "veg100m", "sinuosity", "Road_AADT", "wmean.kangaroo","wmean.koala","wmean.wallaby", "Road_Speed","Road_Width","Length")

#SET UP JAGS MODEL

JAGS.Data <- list(U = length(Data[,"Wallabies_counts"]), K = Data[,"Wallabies_counts"], LENGTH = Data[,"Length"],
                  X = as.matrix(Data[c("veg100m", "sinuosity", "Road_AADT", "Road_Speed","Road_Width")]),
                  Y = as.matrix(Data[,c("veg500", "wmean.wallaby")]),
                  Nx = ncol(Data[,c("veg100m", "sinuosity", "Road_AADT", "Road_Speed","Road_Width")]),
                  Ny = length(Data[,c("veg500", "wmean.wallaby")]))

#get functions
source(paste0(riskmodel__directoryName,"/", "functions_jags.R"))

#run jags
#sfInit( parallel = TRUE, cpus = 3)

#export data, functions and libraries to workers
#sfExportAll()
#sfClusterEval(library(runjags))
#sfClusterEval(library(coda))
#sfClusterEval(library(rjags))
#sfClusterEval(library(parallel))
#sfClusterEval(library(modeest))

#Jags.Fits <- sfLapply(JAGS.Data,get.jags)

#set JAGS location - MODIFY FOR YOUR CASE
JagsLoc <- "C:/ProgramData/Microsoft/Windows/Start Menu/Programs/JAGS"

system.time(Jags.Fit_wallaby <- get.jags(JAGS.Data, JagsLoc))

#sfStop()

save(Jags.Fit_wallaby, file = "jags_fit_wallaby.RData")
