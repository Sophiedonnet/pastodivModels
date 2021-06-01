
source('functions/modules_simulation.R')
source('functions/complete_simulation.R')
library(igraph,ggplot2)
library(tidyverse)
######  General parameters ################# 



n.generations <- 50 # Nb de générations


# Parameters 
param.default <- list(n.ram = 5,
                      n.ewe = 40,
                      career.ram = 8,
                      career.ewe = 8,
                      age.min.ram = 0,
                      age.min.ewe = 0,
                      age.repro.ewe = 3,
                      age.repro.ram = 1)
param.default$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
names(param.default$rate.repro) = c('nb.lambs','probability')



####################################################################""
#################" SIMULATION 1 : 1 herd with accident
#######################################################################"

# Network
n.herds  <- 1 #Nb de troupeaux
param= lapply(1:n.herds,function(i) param.default) # here same parameters for all the Herds. 

# Simulation 
res <- Simulate.herds(n.herds,n.generations = 50,param.allHerds=param,herds.Network = herds.Network,LHerds=NULL,computeInbreeding  = TRUE)
LHerds <- res$LHerds

# InBreeding  along time
inBreedingTime <- res$inBreeding
inBreedingTime$herd <- as.factor(inBreedingTime$herd)
inBreedingTime$gen <- as.factor(inBreedingTime$gen)
ggplot(inBreedingTime,aes(col=gen,y=inBreed,x=gen)) + geom_boxplot() + ggtitle('One herd') +  theme(legend.position='none')

ggsave("slides/alongtime_oneherd.png")

