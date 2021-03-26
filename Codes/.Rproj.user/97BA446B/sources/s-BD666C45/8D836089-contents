
source('functions/modules_simulation.R')
source('functions/complete_simulation.R')
# fonction generale
n.herds  <- 10

herds.Network <- matrix(sample(c(0,1),n.herds^2,replace=TRUE,prob = c(2/3,1/3)),n.herds,n.herds)
diag(herds.Network) <- 0
test  = sum((rowSums(herds.Network)==0))>0
while(test){
  herds.Network <- matrix(sample(c(0,1),n.herds^2,replace=TRUE,prob = c(2/3,1/3)),n.herds,n.herds)
  diag(herds.Network) <- 0
  test  = sum((rowSums(herds.Network)==0)) >0
}

#herds.Network <- diag(1,n.herds)

# essai differents reseaux
n.generations <- 50

LHerds <- Simulate.herds(n.herds,n.generations,param.allHerds=NULL,herds.Network,computeInbreeding = FALSE)

inBreeding <- computeInbreedingFunction(LHerds)  
hist(inBreeding)

Herds.hub = matrix(0,n.herds,n.herds)
Herds.hub[-1,1] =  1
Herds.hub[1,-1] =  1

LHerds <- Simulate.herds(n.herds,n.generations,param.allHerds=NULL,Herds.hub,computeInbreeding = FALSE)
inBreeding <- computeInbreedingFunction(LHerds)  
hist(inBreeding)

LHerds <- Simulate.herds(n.herds,n.generations,param.allHerds=NULL,diag(n.herds),computeInbreeding = FALSE)
inBreeding <- computeInbreedingFunction(LHerds)  
hist(inBreeding)





# essai debug -------------------------------------------------------------


param.default <- list(n.ram = 1,
                      n.ewe = 10,
                      career.ram = 8,
                      career.ewe = 8,
                      age.min.ram = 0,
                      age.min.ewe = 0,
                      age.repro.ewe = 3,
                      age.repro.ram = 1)
param.default$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
names(param.default$rate.repro) = c('nb.lambs','probability')
param= lapply(1:n.herds,function(i) param.default)
LHerds <- Simulate.herds(n.herds,10,param.allHerds=param,diag(n.herds),computeInbreeding = FALSE)
