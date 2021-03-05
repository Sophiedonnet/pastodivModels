
####################################################################
module.initialize.oneHerd <- function(num.herd,param=list(),seed=NULL){
####################################################################
# num.herd num to identify the herd
# param list of useful parameters for initialising the herd
#----------------------------------------------
  if(!is.null(seed)){set.seed(seed)}
  param.default <- list(n.ram = 2,
                        n.ewe = 40,
                        career.ram = 8,
                        age.min.ram=0,
                        career.ewe = 8,
                        age.min.ewe = 0,
                        age.repro.ewe = 3
  )
  #---------------------------------------------
  param.default[names(param)] <- param
  
  pop.table = data.frame(id=paste("0",num.herd,1:(param.default$n.ewe+param.default$n.ram),sep="-"),
                         father = "0",mother = "0",
                         herd=num.herd,sex=c(rep("M",param.default$n.ram),rep("F",param.default$n.ewe)),
                         age = c(sample(param.default$age.min.ram:(param.default$career.ram-1),replace = T,param.default$n.ram),
                                 sample(param.default$age.min.ewe:(param.default$career.ewe-1),replace=T,param.default$n.ewe)))
# - 1 car on les fait vieillir d'un an, que pour ram ? pour etre sur d'avoir des reproducteurs en age
  return(pop.table)
}


####################################################################
module.aging.oneHerd = function(pop.table,param=list()){
####################################################################
  
  if(!is.data.frame(pop.table)){stop()}
  # age  +1 for any animal. 
  # if age > career.ram or  age > career.ewe, we set herd  = -1, meaning that they are leaving the herds forever. 
  
  param.default <- list(n.ram = 2,
                        n.ewe = 40,
                        career.ram = 8,
                        age.min.ram=0,
                        career.ewe = 8,
                        age.min.ewe = 0,
                        age.repro.ewe = 3
  )
  param.default[names(param)] <- param
  #---------------------------------------------
 
  

  pop.table$age <- pop.table$age  + 1
  wM <- which((pop.table$sex=='M') & (pop.table$age > param.default$career.ram))
  wF <- which((pop.table$sex=='F') & (pop.table$age > param.default$career.ewe))
  
  ### For rams and ewes too old : removed from to the Herds. pop.table$herd set to -1 
  pop.table$herd[wM] = - 1
  pop.table$herd[wF] = - 1
  
  return(pop.table)
}


####################################################################
# reproduction intra troupeau
####################################################################
module.reproduction.intraHerd = function(pop.table,num.herd,num.gen,param=list())
{
  ######################################################
  param.default <- list(n.ram = 2,
                        n.ewe = 40,
                        career.ram = 8,
                        age.min.ram=0,
                        career.ewe = 8,
                        age.min.ewe = 0,
                        age.repro.ewe = 3,
                        age.repro.ram = 1
  )
  param.default$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
  names(param.default$rate.repro) = c('nb.lambs','probability')
  param.default[names(param)] <- param
  ########################################################
  
  
  repro.ewe.table <- pop.table[pop.table$sex == "F" & (pop.table$age > param.default$age.repro.ewe),] # ewe able to do babies
  repro.ram.table <- pop.table[pop.table$sex == "M" & (pop.table$age > param.default$age.repro.ram),] # ram able to do babies
  
  n.repro.ewe <- nrow(repro.ewe.table)             # nbre de ewe able to do babies
  n.newborns.per.ewe <- sample(param.default$rate.repro[,1],n.repro.ewe,replace = TRUE,prob  = param.default$rate.repro[,2]) # nb of newborns per ewe.
  n.newborns <- sum(n.newborns.per.ewe)
  # creation table newborn
  newborn.table <- as.data.frame(matrix(data = NA, nrow = n.newborns , ncol=ncol(pop.table)))
  colnames(newborn.table) <- colnames(pop.table)
  
  # creation et attribution new id:
  newborn.table$id=paste(num.gen,num.herd,1:n.newborns,sep="-")
  newborn.table$age <- 0
  newborn.table$sex <- sample(c("F","M"), size = n.newborns, replace = T)
  newborn.table$herd <- num.herd; 
  newborn.table$mother <- rep(repro.ewe.table$id,n.newborns.per.ewe)
  possible.father <- repro.ram.table$id # extraction id father :
  newborn.table$father <- rep(sample(possible.father,size = n.repro.ewe,replace = T),n.newborns.per.ewe)
  return(newborn.table)
}

 
  
 
  
n.herds <- 3
param.allHerds<- lapply(1:n.herds, function(i){list()})
pop.table.allHerds <-  lapply(1:n.herds,function(i){module.initialize.oneHerd(num.herd=i,param=param.allHerds[[i]],seed=i*10)})
pop.table.allHerds <- lapply(1:n.herds, function(i){module.aging.oneHerd(pop.table.allHerds[[i]],param= param.allHerds[[i]])})
pop.table <- pop.table.allHerds[[1]]
newborn.table <- module.reproduction.intraHerd(pop.table.allHerds[[2]],2,1)




newborn.table.ewe <- newborn.table[newborn.table$sexe=="F",]
n.newborn.ewe <- nrow(newborn.table.ewe)
# control ewe ages and remove old ones
newpop.table.ewe <- pop.table[!pop.table$age == 8 & pop.table$sexe=="F",]
n.oldewe <- nrow(pop.table[pop.table$age == 8 & pop.table$sexe=="F",])



############################################# ESSAI






# fonction generale
n.generations <- 1
Simulate.herds = function(n.herds,n.generations,param.allHerds,herds.Network)
{
  
  ################### initialisation   = génération 0 
  LHerds = lapply(1:n.herds,function(i){module.initialize.oneHerd(i,param = param.allHerds[[i]])})
  
  for (gen in 1:nb.generations){ 
    
    LHerds <- lapply(1:n.herds,function(i){module.aging.oneHerd(pop.table = LHerds[[i]],param = param.allHerds[[i]])})
    newBorns <- lapply(1:n.herds,function(i){module.reproduction.intraHerd(pop.table = LHerds[[i]],i,gen,param = param.allHerds[[i]])})
    
    }
    
  
  for (t in 1:nherds)
  {
    
    module.aging()
    module.reproduction()
    module.repartition()
    save(LHerds)
    
  }
  
}
