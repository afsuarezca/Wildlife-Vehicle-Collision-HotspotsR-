

#------------------------------------------------------------------------------#

# Author: Andres Felipe Suarez-Castro
# Date: November 2024
#Control file to run the workflow as described in the report

#IDENTIFY AND PRIORITISE WILDLIFE-VEHICLE COLLISION HOTSPOTS FOR PRIORITY SPECIES IN THE CITY OF LOGAN 
#(REFERENCE NUMBER VP120463)

#------------------------------------------------------------------------------#


#---------------------------------------------------------------------------#
### MAIN DIRECTORY AND PACKAGES ####
#---------------------------------------------------------------------------#

# Add a relative or absolute path name (No spaces are allowed)
path_master <- "C:/Users/s2993062/documents/Logan_project/LCC_wildlife_vehicle_collision"

#Add a folder path for the data######
#I recommend that you save your data in the following folders:

KDEdata_directoryName<-"Data/KDEdata" # Data for Kernel density estimates
EMdata_directoryName<-"Data/ExposureModeldata" #Data for exposure model
riskmodel__directoryName<-"RiskModel"
results_directoryName<-"results"

# Load packages#####

#devtools::install_github("JeremyGelb/spNetwork")

x <- c("sf","terra","tmap","dbscan",
            "lwgeom","RColorBrewer","spNetwork","dplyr",
       "tidyr","plyr","rjags","runjags","blavaan","psych","tidyverse")


lapply(x, require, character.only = TRUE)

#install packages:
# if you don't have these packages installed yet, please use the install.packages("package_name") command.

#----------------------------------------------------#
### KERNEL DENSITY DATA AND ESTIMATES ####
#----------------------------------------------------#

#read road and roadkill data####
#Road network
roads<-"LCCroadsSample.shp"

#Roadkills
roadkills<-"roadkillsSample.shp"

#Plot data?

plotroad<-"yes" #yes or no

#Function parameters#####

#Road segment length in meters
roadl<-300

#Minimum length at which a road will be cut (check ?lixelize_lines())
mindist<-50

# The kernel bandwidth (using the scale of the lines), 
# can be a single float or a numeric vector 
# if a different bandwidth must be used for each event.
bw<-300

#Method
# The method to use when calculating the NKDE, 
# must be one of simple / discontinuous / continuous 
# (see nkde details for more information)

method<-"discontinuous"

#kernel_name
# The name of the kernel to use. Must be one of triangle, gaussian, tricube, cosine, 
# triweight, quartic, epanechnikov or uniform.

kernel<-"gaussian"

#---------------------------------------------------#
# Apply the Network Kernel Density Estimation#####
#--------------------------------------------------#

source(paste0(path_master,"/", "KDE_function_LCC.R"), 
       echo=TRUE)

#-------------------------------------------#
### RISK MODEL VARIABLES####
#-------------------------------------------#

#specify if it is necessary to generate the table with variables used in the exposure model
generateTable<-"yes" #yes or no. If no, you must specify file name in myvar 

#_____________________________________________________#
#If you have table with variables, add file name here:
#_____________________________________________________#

#Your table MUST have the following variables with the exact names:

c("Road_Block","veg500", "veg100m", "sinuosity", "Road_AADT", 
  "wmean.kangaroo","wmean.koala","wmean.wallaby",
  "Road_Speed","Road_Width","Length", "Kangaroo_counts",
  "Koala_counts","Wallabies_counts")


myvar<-"variables_road_logan.csv"

#_________________________________________________#
#Otherwise, add raster names to calculate variables
#_________________________________________________#


#add path, file names for exposure model

#Species distribution model. Developed using code "Species occurrence submodel"
spoc <- "Occurrence_models/Final_models/Macropus_rufogriseus_noimp.tif"

#Land cover categories as explained in table 1 of the final report
lccat<- "LandCoverCategories/BroadLandCoverCategories_fromSAVIandHeightReclassified_20170626.tif"

#Impervious surface
impervsur<- "LandCoverCategories/imp_sur"

#Vegetation raster
veg<-"localveglcc"

#value (in meters) to calculate variables at buffers around roads
bufferoads<- 100

#lambda value to calculate weighted distance
lambdasp<-0.005

#table with road variables and roadkill data
#It must have the following columns
#"Road_Block"        "Road_Speed"        "Road_Width"        "Road_AADT"        
#"Length"            "wmean.koala"       "wmean.kangaroo"    "wmean.wallaby"     "veg500"            "veg100m"          
#"sinuosity"         "Wallabies_counts"  "Koala_counts"      "Kangaroo_counts"

roadvar<-"roadvar.csv"

source(paste0(path_master,"/", "calculateVariablesAroundSegments.R"), 
       echo=TRUE)


#__________________________________________________#
#-----------Run risk model---------------------####
#_________________________________________________#

#Specify location of your text file for your model to run in JAGS
model = "jags_model.txt"
#set JAGS location - MODIFY FOR YOUR CASE
JagsLoc <- "C:/ProgramData/Microsoft/Windows/Start Menu/Programs/JAGS"

source(paste0(riskmodel__directoryName,"/","jags_code.R"))

