rm(list=(ls))
source('functions/modules_simulation.R')
source('functions/complete_simulation.R')
library(igraph,ggplot2)

######  General parameters ################# 



n.generations <- 30 # Nb de générations


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
n.herds  <- 2 #Nb de troupeaux
param= lapply(1:n.herds,function(i) param.default) # here same parameters for all the Herds. 

# Simulation 
res <- Simulate.herds(n.herds,n.generations,param.allHerds=param,herds.Network = NULL,LHerds=NULL,computeInbreeding  = TRUE)
LHerds <- res$LHerds

# InBreeding  along time
inBreedingTime <- res$inBreeding
inBreedingTime$gen <- as.factor(inBreedingTime$gen)
ggplot(inBreedingTime,aes(col=gen,y=inBreed,x=gen)) + geom_boxplot() + ggtitle('One herd') +  theme(legend.position='none')
ggplot(inBreedingTime,aes(y=mean_inbreed,x=gen,group=herd)) + geom_line(aes(linetype=herd, color=herd)) + ggtitle('Star + one isolated ') 


############ accident
LHerds[[1]] <- module.lose(LHerds[[1]],loseProportion = 0.5)
compute.herds.size(LHerds)
res2 <- Simulate.herds(n.herds,10,param.allHerds=param,herds.Network = NULL,LHerds=LHerds,computeInbreeding  = TRUE)


inBreedingTime <- res2$inBreeding
inBreedingTime <- inBreedingTime %>% group_by(herd,gen) %>%summarise(mean_inbreed = mean(inBreed)) %>%mutate(herd = as.factor(herd))
ggplot(inBreedingTime,aes(y=mean_inbreed,x=gen,group=herd)) + geom_line(aes(linetype=herd, color=herd)) + ggtitle('Star + one isolated ') 
ggsave("slides/alongtime_meantenherds.png")




size.herds <- res2$herds_size 
ggplot(res2$herds_size,aes(col=herd,y=size,x=gen)) + geom_line() + ggtitle('Two herds after accident') +  theme(legend.position='none')


