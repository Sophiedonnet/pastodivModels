
source('functions/modules_simulation.R')

# fonction generale
n.herds  <- 10
herds.Network <- matrix(sample(c(0,1),n.herds^2,replace=TRUE,prob = c(2/3,1/3)),n.herds,n.herds)
diag(herds.Network) <- 0

n.generations <- 1

myHerds <- Simulate.herds(n.herds,n.generations,param.allHerds=NULL,herds.Network)
  
