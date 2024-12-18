rm(list = ls())
#------------------------------


# Set working directory
setwd("//uq.edu.au/UQ-Research/LOGANROAD-A1321")
#install packages:
# if you dont have these packages installed yet, please use the install.packages("package_name") command.

data<- read.csv("variables_roads_logan.csv")

data<-subset(data,data$counts!=0)

data<-data[c("counts","veg500", "sinuosity", "Road_AADT", "wmean.kangaroo")]

scaled.dat <- scale(data)

hayesrisk <- ' # direct effect

counts ~ veg500 + sinuosity + Road_AADT + wmean.kangaroo

#show that dependent variable has variance
                    counts ~~ counts'

data<-data.frame(counts = rnorm(100,1,1),veg500 =rnorm(100,1,1) , sinuosity=rnorm(100,1,1), Road_AADT = rnorm(100,1,1), wmean.kangaroo = rnorm(100,1,1))


fit.bayes <- blavaan(hayesrisk, data=data, convergence= "auto",  test="none", seed=c(19,02,2018))

summary(fit.bayes)

model.informative.priors1 <- 
  '#the regression model with priors
counts ~ prior("dnorm(0,0.001)")*sinuosity + prior("dnorm(0,0.001)")*veg500 + prior("dnorm(0,0.001)")*Road_AADT + prior("dnorm(0,0.001)")*Road_AADT

counts ~~ counts

#we want to have an intercept (with normal prior)
counts ~ 1'

fit.bayes <- blavaan(model.informative.priors1, data=data[sample(nrow(data),1000),], convergence= "auto",  test="none", seed=c(19,02,2018))

# the test="none" input stops the calculations of some posterior checks, we do not need at this moment and speeds up the process. 
# the seed command is simply to guarantee the same exact result when running the sampler multiple times. You do not have to set this.