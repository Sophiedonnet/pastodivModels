
source('functions/modules_simulation.R')
source('functions/complete_simulation.R')
library(igraph,ggplot2)
library(tidyverse)
######  General parameters ################# 

 

n.generations <- 50 # Nb de générations


# Parameters 
param.default <- list(n.ram = 5,
                      n.ewe = 80,
                      career.ram = 8,
                      career.ewe = 8,
                      age.min.ram = 0,
                      age.min.ewe = 0,
                      age.repro.ewe = 3,
                      age.repro.ram = 1)
param.default$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
names(param.default$rate.repro) = c('nb.lambs','probability')



####################################################################""
#################" SIMULATION 1 : independent herds : no exchange
#######################################################################"

# Network
n.herds  <- 1 #Nb de troupeaux
param= lapply(1:n.herds,function(i) param.default) # here same parameters for all the Herds. 

ram.Network <- diag(1,n.herds)
plot(graph_from_adjacency_matrix(ram.Network, mode = c("directed")))
herds.Network = list(ewe.for.replace= NULL,ram.for.replace = ram.Network,ram.for.repro = ram.Network)


# Simulation 
seed = sample(1:100,1)
set.seed(seed)
res <- Simulate.herds(n.herds,n.generations,param.allHerds=param,herds.Network = herds.Network,LHerds=NULL,computeInbreeding  = TRUE)
LHerds <- res$LHerds

# InBreeding  along time
inBreedingTime <- res$inBreeding
inBreedingTime$herd <- as.factor(inBreedingTime$herd)
inBreedingTime$gen <- as.factor(inBreedingTime$gen)
ggplot(inBreedingTime,aes(col=gen,y=inBreed,x=gen)) + geom_boxplot() + ggtitle('One herd') +  theme(legend.position='none')

ggsave("slides/alongtime_oneherd.png")

 


#######################################################################"
#################" SIMULATION 2 : Hub : one people gives to all the others
#######################################################################"
n.herds = 10;
param= lapply(1:n.herds,function(i) param.default)

ram.Network <- diag(0,n.herds)
ram.Network[,1] <- 1
ram.Network[n.herds,1] <- 0
ram.Network[n.herds,n.herds] <- 1
herds.Network = list(ewe.for.replace= NULL,ram.for.replace =ram.Network,ram.for.repro = ram.Network+diag(n.herds))
plot(graph_from_adjacency_matrix(t(ram.Network), mode = c("directed")))
plot(graph_from_adjacency_matrix(t(herds.Network$ram.for.repro), mode = c("directed")))

# Simulation 
res <- Simulate.herds(n.herds,n.generations,param.allHerds = param,herds.Network = herds.Network,LHerds = NULL,computeInbreeding  = TRUE)
LHerds <- res$LHerds
inBreedingTime <- res$inBreeding
inBreedingTime <- inBreedingTime %>% group_by(herd,gen) %>%summarise(mean_inbreed = mean(inBreed)) %>%mutate(herd = as.factor(herd))
ggplot(inBreedingTime,aes(y=mean_inbreed,x=gen,group=herd)) + geom_line(aes(linetype=herd, color=herd)) + ggtitle('Star + one isolated ') 
ggsave("slides/alongtime_meantenherds.png")



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
plot(graph_from_adjacency_matrix(t(ram.Network), mode = c("directed")))
herds.Network = list(ram.for.replace = ram.Network,ram.for.repro = ram.Network, ewe.for.replace =NULL)

res  <- Simulate.herds(n.herds ,n.generations,param.allHerds = param,herds.Network)
LHerds <- res$LHerds
inBreeding <- computeInbreedingFunction(LHerds)  
inBreeding$herd <- as.factor(inBreeding$herd)
ggplot(inBreeding,aes(col=herd,y=inBreed,x=herd)) + geom_boxplot()  + ggtitle('Random network')



#################" SIMULATION 4 :  chain


ram.Network <- diag(0,n.herds)
for (i in 2:n.herds){
  ram.Network[i,i-1] <- 1
}
ram.Network[1,n.herds-1] <- 1
ram.Network[n.herds,n.herds-1] <- 0
ram.Network[n.herds,n.herds] <- 1
herds.Network = list(ram.for.replace = ram.Network,ram.for.repro = diag(n.herds)+ram.Network, ewe.for.replace =NULL)
plot(graph_from_adjacency_matrix(t(ram.Network), mode = c("directed")))


res <- Simulate.herds(n.herds ,n.generations,param.allHerds = param,herds.Network)
LHerds <- res$LHerds
inBreeding <- computeInbreedingFunction(LHerds)  
inBreeding$herd <- as.factor(inBreeding$herd)
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



