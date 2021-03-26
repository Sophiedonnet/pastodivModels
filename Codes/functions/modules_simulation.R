
####################################################################
module.initialize.oneHerd <- function(num.herd,param,seed=NULL){
  ####################################################################
  # num.herd num to identify the herd
  # param list of useful parameters for initialising the herd
  #----------------------------------------------
  if(!is.null(seed)){set.seed(seed)}
  
  
  pop.table = data.frame(ind=paste("0",num.herd,1:(param$n.ewe+param$n.ram),sep="-"),
                         father = "0",mother = "0",
                         herd=num.herd,
                         sex=c(rep("M",param$n.ram),rep("F",param$n.ewe)),
                         age = c(sample(param$age.min.ram:(param$career.ram-1),replace = T,param$n.ram),
                                 sample(param$age.min.ewe:(param$career.ewe-1),replace=T,param$n.ewe)))
  return(pop.table)
}


####################################################################
module.aging.oneHerd = function(pop.table,param=list()){
  ####################################################################
  
  if(!is.data.frame(pop.table)){stop()}
  # age  +1 for any animal. 
  # if age > career.ram or  age > career.ewe, we set herd  = -1, meaning that they are leaving the herds forever. 
  
  # param <- list(n.ram = 2,
  #                       n.ewe = 40,
  #                       career.ram = 8,
  #                       age.min.ram=0,
  #                       career.ewe = 8,
  #                       age.min.ewe = 0,
  #                       age.repro.ewe = 3
  # )
  # param[names(param)] <- param
  #---------------------------------------------
  pop.table$age <- pop.table$age  + 1
  #wM <- which((pop.table$sex=='M') & (pop.table$age > param$career.ram))
  #wF <- which((pop.table$sex=='F') & (pop.table$age > param$career.ewe))
  
  ### For rams and ewes too old : removed from to the Herds. pop.table$herd set to -1 
  #pop.table$herd[wM] = - 1
  #pop.table$herd[wF] = - 1
  
  return(pop.table)
}


###########reproduction intra troupeau  ################################
module.reproduction.intraHerd = function(pop.table,num.gen,param=list())
  ####################################################################
{
  
  ########################################################
  num.herd <- pop.table$herd[which(pop.table$herd != -1)[1]]
  ######################################################## 
  
  repro.ewe.table <- pop.table[(pop.table$herd != -1) & (pop.table$sex == "F") & (pop.table$age >= param$age.repro.ewe),] # ewe able to do babies
  repro.ram.table <- pop.table[(pop.table$herd != -1) & (pop.table$sex == "M") & (pop.table$age >= param$age.repro.ram),] # ram able to do babies
  
  n.repro.ewe <- nrow(repro.ewe.table)             # nbre de ewe able to do babies
  n.newborns.per.ewe <- sample(param$rate.repro[,1],n.repro.ewe,replace = TRUE,prob  = param$rate.repro[,2]) # nb of newborns per ewe.
  n.newborns <- sum(n.newborns.per.ewe)
  # creation table newborn
  newborn.table <- as.data.frame(matrix(data = NA, nrow = n.newborns , ncol=ncol(pop.table)))
  colnames(newborn.table) <- colnames(pop.table)
  # fullfilling newborn.table
  newborn.table$ind=paste(num.gen,num.herd,1:n.newborns,sep="-")
  newborn.table$age <- 0
  newborn.table$sex <- sample(c("F","M"), size = n.newborns, replace = T)
  newborn.table$herd <- num.herd; 
  newborn.table$mother <- rep(repro.ewe.table$ind,n.newborns.per.ewe)
  possible.father <- repro.ram.table$ind # extraction id father :
  newborn.table$father <- rep(sample(possible.father,size = n.repro.ewe,replace = T),n.newborns.per.ewe)
  return(newborn.table)
}


####################################################################
# reproduction inter troupeaux
####################################################################
module.reproduction = function(mothers,fathers,num.gen,param=list())
{
  
  ########################################################
  num.herd <- mothers$herd[which(mothers$herd != -1)[1]]
  ######################################################## 
  
  n.mothers <- nrow(mothers)             # nbre de ewe able to do babies
  possible.father <- fathers$ind
  
  if ((length(possible.father)>0) & (n.mothers > 0)){
    
    n.newborns.per.mother <- sample(param$rate.repro[,1],n.mothers,replace = TRUE,prob  = param$rate.repro[,2]) # nb of newborns per ewe.
    n.newborns <- sum(n.newborns.per.mother)
    # creation table newborn
    newborn.table <- as.data.frame(matrix(data = NA, nrow = n.newborns , ncol=ncol(mothers)))
    colnames(newborn.table) <- colnames(mothers)
    # fullfilling newborn.table
    newborn.table$ind=paste(num.gen,num.herd,1:n.newborns,sep="-")
    newborn.table$age <- 0
    newborn.table$sex <- sample(c("F","M"), size = n.newborns, replace = T)
    newborn.table$herd <- num.herd; 
    newborn.table$mother <- rep(mothers$ind,n.newborns.per.mother)
    newborn.table$father <- rep(sample(possible.father,size = n.mothers,replace = T),n.newborns.per.mother)
  }else{
    newborn.table = NULL
  }
  return(newborn.table)
}



####################################################################
# replacement intra troupeau
####################################################################
module.replaceEwe.intraHerd = function(pop.table,newborn.table,param=list()){
  
  # param <- list(n.ram = 2,
  #                       n.ewe = 40,
  #                       career.ram = 8,
  #                       age.min.ram=0,
  #                       career.ewe = 8,
  #                       age.min.ewe = 0,
  #                       age.repro.ewe = 3,
  #                       age.repro.ram = 1
  # )
  # param$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
  # names(param$rate.repro) = c('nb.lambs','probability')
  # param[names(param)] <- param
  
  
  ################"" REPLACE the ewe that are too old. 
  
  w.F <- which(newborn.table$sex == 'F')
  w.TooOld <- which(pop.table$herd != -1 & pop.table$sex == 'F' & pop.table$age >= param$career.ewe)
  n.TooOld <- length(w.TooOld)
  
  n.Lacking <- param$n.ewe - sum((pop.table$herd) != -1 & (pop.table$sex == 'F'))
  if ((n.TooOld > 0) & (length(w.F)>0)){
    pop.table$herd[w.TooOld]<- -1
    u <- sample(w.F,min(n.TooOld + n.Lacking,length(w.F)),replace=FALSE)
    pop.table <- rbind(pop.table,newborn.table[u,])
    newborn.togive <- newborn.table[-u,]
  } else {
    newborn.togive <- newborn.table
  }
  
  
    res <- list(pop.table = pop.table,newborn.togive  = newborn.togive)
  return(res)
}



####################################################################
# choose rams for each herd
####################################################################  
choose.Ram.reproduction <- function(network,n,mode = "intra") 
{
  switch(mode,
         intra = {return(as.list(1:n))},
         melange = {}
  )
}


choose.Ram.replace <- function(network,n,mode = "intra") 
{
  switch(mode,
         intra = {return(as.list(1:n))},
         melange = {}
  )
}


####################################################################
# choose rams for each herd
####################################################################
computeInbreedingFunction = function(LHerds){
  
  n.herds <- length(LHerds)
  size.herds <- sapply(1:n.herds,function(i){nrow(LHerds[[i]])})
  allHerds <- do.call("rbind",LHerds)
  ped <- allHerds
  ped=ped[!duplicated(ped$ind),]
  
  LEV = unique(c(( allHerds$ind),( allHerds$father),( allHerds$mother)))
  ped$ind=as.numeric(factor(ped$ind,levels=LEV))
  ped$father=as.numeric(factor(ped$father,levels=LEV))
  ped$mother=as.numeric(factor(ped$mother,levels=LEV))
  
  
  w.without.knwon.parents  = which(ped$father==which(LEV == "0"))
  ped$father[w.without.knwon.parents] = 0
  ped$mother[w.without.knwon.parents] = 0
  
  w.M <- which(ped$sex=='M')
  w.F <- which(ped$sex=='F')
  ped$sex[w.M] <- 1
  ped$sex[w.F] <- 2
  geneal <- gen.genealogy(ped)
  u <- which(ped$herd!=-1)
  inBreed <- gen.f(geneal,u)
}




