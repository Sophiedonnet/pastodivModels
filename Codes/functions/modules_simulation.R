
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
module.reproduction.intraHerd = function(pop.table,num.gen,param=list())
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
  num.herd <- pop.table$herd[which(pop.table$herd != -1)[1]]
  ######################################################## 
  
  repro.ewe.table <- pop.table[(pop.table$herd != -1) & (pop.table$sex == "F") & (pop.table$age > param.default$age.repro.ewe),] # ewe able to do babies
  repro.ram.table <- pop.table[(pop.table$herd != -1) & (pop.table$sex == "M") & (pop.table$age > param.default$age.repro.ram),] # ram able to do babies
  
  n.repro.ewe <- nrow(repro.ewe.table)             # nbre de ewe able to do babies
  n.newborns.per.ewe <- sample(param.default$rate.repro[,1],n.repro.ewe,replace = TRUE,prob  = param.default$rate.repro[,2]) # nb of newborns per ewe.
  n.newborns <- sum(n.newborns.per.ewe)
  # creation table newborn
  newborn.table <- as.data.frame(matrix(data = NA, nrow = n.newborns , ncol=ncol(pop.table)))
  colnames(newborn.table) <- colnames(pop.table)
  # fullfilling newborn.table
  newborn.table$id=paste(num.gen,num.herd,1:n.newborns,sep="-")
  newborn.table$age <- 0
  newborn.table$sex <- sample(c("F","M"), size = n.newborns, replace = T)
  newborn.table$herd <- num.herd; 
  newborn.table$mother <- rep(repro.ewe.table$id,n.newborns.per.ewe)
  possible.father <- repro.ram.table$id # extraction id father :
  newborn.table$father <- rep(sample(possible.father,size = n.repro.ewe,replace = T),n.newborns.per.ewe)
  return(newborn.table)
}


####################################################################
# reproduction inter troupeaux
####################################################################
module.reproduction.interHerds = function(pop.table.ewe,pop.table.ram,num.gen,param=list())
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
  num.herd <- pop.table.ewe$herd[which(pop.table.ewe$herd != -1)[1]]
  ######################################################## 
  
  repro.ewe.table <- pop.table.ewe[(pop.table.ewe$herd != -1) & (pop.table.ewe$sex == "F") & (pop.table.ewe$age > param.default$age.repro.ewe),] # ewe able to do babies
  repro.ram.table <- pop.table.ram[(pop.table.ram$herd != -1) & (pop.table.ram$sex == "M") & (pop.table.ram$age > param.default$age.repro.ram),] # ram able to do babies
  
  n.repro.ewe <- nrow(repro.ewe.table)             # nbre de ewe able to do babies
  n.newborns.per.ewe <- sample(param.default$rate.repro[,1],n.repro.ewe,replace = TRUE,prob  = param.default$rate.repro[,2]) # nb of newborns per ewe.
  n.newborns <- sum(n.newborns.per.ewe)
  # creation table newborn
  newborn.table <- as.data.frame(matrix(data = NA, nrow = n.newborns , ncol=ncol(pop.table.ewe)))
  colnames(newborn.table) <- colnames(pop.table.ewe)
  # fullfilling newborn.table
  newborn.table$id=paste(num.gen,num.herd,1:n.newborns,sep="-")
  newborn.table$age <- 0
  newborn.table$sex <- sample(c("F","M"), size = n.newborns, replace = T)
  newborn.table$herd <- num.herd; 
  newborn.table$mother <- rep(repro.ewe.table$id,n.newborns.per.ewe)
  possible.father <- repro.ram.table$id # extraction id father :
  newborn.table$father <- rep(sample(possible.father,size = n.repro.ewe,replace = T),n.newborns.per.ewe)
  return(newborn.table)
}



####################################################################
# replacement intra troupeau
####################################################################
module.replaceEwe.intraHerd = function(pop.table,newborn.table,param=list()){
  
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
  
  
  ################"" REPLACE the ewe that are too old. 
  w.F <- which(newborn.table$sex == 'F')
  w.TooOld <- which(pop.table$herd != -1 & pop.table$sex == 'F' & pop.table$age >= param.default$career.ewe)
  n.TooOld <- length(w.TooOld)
  if (n.TooOld > 0){
    pop.table$herd[-w.TooOld]<- -1
    u <- sample(w.F,n.TooOld,replace=FALSE)
    pop.table <- rbind(pop.table,newborn.table[u,])
    newborn.togive <- newborn.table[-u,]
  }
  res <- list(pop.table = pop.table,newborn.togive  = newborn.togive)
  return(res)
}



####################################################################
# choose rams for each herd
####################################################################  
choose.Ram <- function(network,n,mode = "intra") 
{
  switch(mode,
         intra = {return(as.list(1:n))},
         melange = {}
         )
  
}

# 
#   
# n.herds <- 3
# param.allHerds<- lapply(1:n.herds, function(i){list()})
# pop.table.allHerds <-  lapply(1:n.herds,function(i){module.initialize.oneHerd(num.herd=i,param=param.allHerds[[i]],seed=i*10)})
# pop.table.allHerds <- lapply(1:n.herds, function(i){module.aging.oneHerd(pop.table.allHerds[[i]],param= param.allHerds[[i]])})
# pop.table <- pop.table.allHerds[[1]]
# newborn.table <- module.reproduction.intraHerd(pop.table.allHerds[[1]],1,1)
# 
# 
# 
# 
# newborn.table.ewe <- newborn.table[newborn.table$sexe=="F",]
# n.newborn.ewe <- nrow(newborn.table.ewe)
# # control ewe ages and remove old ones
# newpop.table.ewe <- pop.table[!pop.table$age == 8 & pop.table$sexe=="F",]
# n.oldewe <- nrow(pop.table[pop.table$age == 8 & pop.table$sexe=="F",])
# 


############################################# ESSAI

herds.Network <- matrix(sample(c(0,1),n.herds^2,replace=TRUE,prob = c(2/3,1/3)),n.herds,n.herds)
diag(herds.Network) <- 0


# fonction generale
n.generations <- 1
Simulate.herds = function(n.herds,n.generations,param.allHerds=NULL,herds.Network,LHerds=NULL)
{
  
  if(is.null(param.allHerds)){
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
    param.allHerds <- lapply(1:n.herds,function(i){param.default})
  }
  
  ################### initialisation   = génération 0 
  if (is.null(LHerds))
    LHerds = lapply(1:n.herds,function(i){module.initialize.oneHerd(i,param = param.allHerds[[i]])})
  
  for (gen in 1:n.generations){ 
    
    LHerds <- lapply(1:n.herds,function(i){module.aging.oneHerd(pop.table = LHerds[[i]],param = param.allHerds[[i]])})
    
    
    pat <- choose.Ram(herds.Network,n.herds)
    newBorns <- lapply(1:n.herds,function(i){module.reproduction.interHerds(pop.table.ewe = LHerds[[i]],pop.table.ram = do.call("rbind",LHerds[pat[[i]]]),gen,param = param.allHerds[[i]])})
    
    
    resultsReplace <- lapply(1:n.herds,function(i){module.replaceEwe.intraHerd(pop.table = LHerds[[i]],newBorns[[i]],param = param.allHerds[[i]])})
    
    LHerds <- lapply(1:n.herds,function(i){resultsReplace[[i]]$pop.table})
    Lnewborns.togive <- lapply(1:n.herds,function(i){resultsReplace[[i]]$newborn.togive})
    
    for (i in 1:n.herds){
      
      pop.table.i <- LHerds[[i]]
      param.i <- param.allHerds[[i]]
      w.R.TooOld <- which((pop.table.i$herd != -1) & (pop.table.i$sex=='M') & pop.table.i$age > param.i$age.repro.ram)  
      n.R.TooOld <- length(w.R.TooOld)
      if (n.R.TooOld > 0){
        donnor.i <-sample(which(herds.Network[i,]==1),1)
        L.i <- Lnewborns.togive[[donnor.i]]
        L.i <- L.i[L.i$sex=='M',][sample(1:nrow(L.i),n.R.TooOld,replace=FALSE),]
        L.i$herd <- i
        pop.table.i$herd[w.R.TooOld]<- -1
        pop.table.i <- rbind(pop.table.i,L.i)
        LHerds[[i]] <- pop.table.i
      }
    }
  }
  return(LHerds)
}
  