rm(list=ls())
source('functions/modules_simulation.R')
source('functions/complete_simulation.R')
library(igraph,ggplot2)
library(tidyverse)
######  General parameters ################# 

 
n.herds  <- 10#Nb de troupeaux
n.generations <- 50 # Nb de générations


# Parameters 
param.default <- list(n.ram = 2,
                      n.ewe = 40,
                      age.max.repro.ram = 8,
                      age.max.repro.ewe = 8,
                      age.min.repro.ewe = 3,
                      age.min.repro.ram = 1,
                      career.ram  = 3)
param.default$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
names(param.default$rate.repro) = c('nb.lambs','probability')
param = lapply(1:n.herds,function(i) param.default) # here same parameters for all the Herds. 

####################################################################""
#################" SIMULATION 1 : independent herds : no exchange
#######################################################################"

# Network

ram.Network <- diag(1,n.herds)
plot(graph_from_adjacency_matrix(ram.Network, mode = c("directed")))
herds.Network = list(ewe.replace= NULL,ram.replace = ram.Network,ram.circulation = NULL)


# Simulation 
myparam <- param
myparam[[1]]$career.ram = 3
param.allHerds <- myparam
res <- Simulate.herds(n.herds,n.generations,param.allHerds,herds.Network = herds.Network,LHerds=NULL,computeInbreeding  = FALSE)
LHerds <- res$LHerds

inBreeding <- compute.inbreeding(LHerds)  
ggplot(inBreeding,aes(x=inBreed)) + geom_histogram()
ggplot(inBreeding,aes(col=herd,y=inBreed,x=herd)) + geom_boxplot() + ggtitle('independant ')


#######################################################################"
#################" SIMULATION 2 : Hub : one people gives to all the others
#######################################################################"
ram.Network <- diag(0,n.herds)
ram.Network[,1] <- 1
ram.Network[n.herds,1] <- 0
ram.Network[n.herds,n.herds] <- 1
herds.Network = list(ewe.replace= NULL,ram.replace = ram.Network)
plot(graph_from_adjacency_matrix(t(ram.Network), mode = c("directed")))


# Simulation 
res <- Simulate.herds(n.herds,n.generations,param.allHerds = param,herds.Network = herds.Network,LHerds = NULL)
LHerds <- res$LHerds
# InBreeding 
inBreeding <- compute.inbreeding(LHerds)  


ggplot(inBreeding,aes(x=inBreed)) + geom_histogram()
ggplot(inBreeding,aes(col=herd,y=inBreed,x=herd)) + geom_boxplot() +   ggtitle('Hub network')

#######################################################################"
#################" SIMULATION 3 : random network
#######################################################################"
set.seed(1)
test  = TRUE
while(test){
  ram.Network <- matrix(sample(c(0,1),(n.herds)^2,replace=TRUE,prob = c(9/10,1/10)),(n.herds),(n.herds))
  diag(ram.Network) <- 0
  ram.Network[1,] <- ram.Network[,1] <- 0
  ram.Network[1,1] <- 1
  test  = sum((rowSums(ram.Network)==0)) >0
}
plot(graph_from_adjacency_matrix(t(ram.Network), mode = c("directed")),main = "Ram")
herds.Network = list(ram.replace = ram.Network, ewe.replace =NULL)

res  <- Simulate.herds(n.herds ,n.generations,param.allHerds = param,herds.Network)
LHerds <- res$LHerds
inBreeding <- compute.inbreeding(LHerds)  
ggplot(inBreeding,aes(col=herd,y=inBreed,x=herd)) + geom_boxplot()  + ggtitle('Random network')



#################" SIMULATION 4 :  chain

param.default <- list(n.ram = 2,
                      n.ewe = 40,
                      age.max.repro.ram = 8,
                      age.max.repro.ewe = 8,
                      age.min.repro.ewe = 3,
                      age.min.repro.ram = 1,
                      career.ram  = 3)
param.default$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
names(param.default$rate.repro) = c('nb.lambs','probability')
param = lapply(1:n.herds,function(i) param.default) # here same parameters for all the Herds. 

ram.Network <- diag(0,n.herds)
for (i in 2:n.herds){
  ram.Network[i,i-1] <- 1
}
ram.Network[1,n.herds-1] <- 1
ram.Network[n.herds,n.herds-1] <- 0
ram.Network[n.herds,n.herds] <- 1
herds.Network = list(ram.circulation = ram.Network, ewe.replace =NULL)
plot(graph_from_adjacency_matrix(t(ram.Network), mode = c("directed")))


res <- Simulate.herds(n.herds ,n.generations,param.allHerds = param,herds.Network)
LHerds <- res$LHerds
inBreeding <- compute.inbreeding(LHerds)  
ggplot(inBreeding,aes(col=herd,y=inBreed,x=herd)) + geom_boxplot()+ ggtitle('Chain network')



#################" SIMULATION 5 :  complete network

n.herds = 5
ram.Network <- diag(0,2*n.herds)
for (i in 2:n.herds){
  ram.Network[i,i-1] <- 1
}
ram.Network[1,n.herds] <- 1
diag(ram.Network)[n.herds + (1:n.herds)] = 1


ewe.Network <- diag(0,2*n.herds)
for (i in n.herds+(2:n.herds)){
  ewe.Network[i,i-1] <- 1
}
ewe.Network[n.herds+1,2*n.herds] <- 1
diag(ewe.Network)[(1:n.herds)] = 1



plot(graph_from_adjacency_matrix(t(ram.Network), mode = c("directed")))
plot(graph_from_adjacency_matrix(t(ewe.Network), mode = c("directed")))


herds.Network = list(ram.for.replace = ram.Network,ram.for.repro = NULL, ewe.for.replace =ewe.Network)

n.herds  =10
res <- Simulate.herds(n.herds ,n.generations,param.allHerds = param,herds.Network)
LHerds <- res$LHerds
inBreeding <- computeInbreedingFunction(LHerds)  
inBreeding$herd <- as.factor(inBreeding$herd)
ggplot(inBreeding,aes(col=herd,y=inBreed,x=herd)) + geom_boxplot()+ ggtitle('Chain network')


#################" SIMULATION 6 :  change size of one herd
param[[1]]$n.ram = 6
param[[1]]$n.ewe = 120

ram.Network <- diag(0,n.herds)
ram.Network[,1] <- 1
ram.Network[n.herds,1] <- 0
ram.Network[n.herds,n.herds] <- 1
plot(graph_from_adjacency_matrix(t(ram.Network), mode = c("directed")))


# Simulation 
LHerds <- Simulate.herds(n.herds,n.generations,param.allHerds = param,herds.Network)

# InBreeding 
inBreeding <- computeInbreedingFunction(LHerds)  
inBreeding$herd <- as.factor(inBreeding$herd)

ggplot(inBreeding,aes(x=inBreed)) + geom_histogram()
ggplot(inBreeding,aes(col=herd,y=inBreed,x=herd)) + geom_boxplot() +   ggtitle('Hub network')


#################" SIMULATION 7 :  compare exchange on ram and exchanges on ewe

n.herds = 5
ram.Network <- diag(0,2*n.herds)
for (i in 2:n.herds){
  ram.Network[i,i-1] <- 1
}
ram.Network[1,n.herds] <- 1
diag(ram.Network)[n.herds + (1:n.herds)] = 1


ewe.Network <- diag(0,2*n.herds)
for (i in n.herds+(2:n.herds)){
  ewe.Network[i,i-1] <- 1
}
ewe.Network[n.herds+1,2*n.herds] <- 1
diag(ewe.Network)[(1:n.herds)] = 1



plot(graph_from_adjacency_matrix(t(ram.Network), mode = c("directed")))
plot(graph_from_adjacency_matrix(t(ewe.Network), mode = c("directed")))


herds.Network = list(ram.for.replace = ram.Network,ram.for.repro = NULL, ewe.for.replace =ewe.Network)

n.herds  =10
res <- Simulate.herds(n.herds ,n.generations,param.allHerds = param,herds.Network)
LHerds <- res$LHerds
inBreeding <- computeInbreedingFunction(LHerds)  
inBreeding$herd <- as.factor(inBreeding$herd)
ggplot(inBreeding,aes(col=herd,y=inBreed,x=herd)) + geom_boxplot()+ ggtitle('Compare exchange on ewe/ram')



