library(plyr)
library(tidyverse)

######### Initialize one herd  ########################
module.initialize.oneHerd <- function(num.herd,param,seed=NULL){
####################################################################
# num.herd num to identify the herd
# param list of useful parameters for initialising the herd
#----------------------------------------------
  if(!is.null(seed)){set.seed(seed)}
  
  if(is.null(param$age.min.ram)){param$age.min.ram = 0}
  if(is.null(param$age.min.ewe)){param$age.min.ewe = 0}
  
  
  pop.table = data.frame(ind=paste("0",num.herd,1:(param$n.ewe+param$n.ram),sep="-"),
                         father = "0",mother = "0",
                         herd=num.herd,
                         sex=c(rep("M",param$n.ram),rep("F",param$n.ewe)),
                         age = c(sample(param$age.min.ram:(param$age.max.repro.ram-1),replace = T,param$n.ram),
                                 sample(param$age.min.ewe:(param$age.max.repro.ewe-1),replace=T,param$n.ewe)))
  
  pop.table$time.in.current.herd <- 0
  return(pop.table)
}
######### Initialize all herds ########################
module.initialize <- function(param.allHerds, seed = NULL){
####################################################################
  n.herds <- length(param.allHerds)
  LHerds <- lapply(1:n.herds,function(i){module.initialize.oneHerd(i,param = param.allHerds[[i]])})
  return(LHerds)
}

############ Vieillissement intra herds ###########################################
module.aging.oneHerd = function(pop.table,param){
####################################################################
  
  if(!is.data.frame(pop.table)){stop()}
  w <- which(pop.table$herd != -1)
  pop.table$age[w] <- pop.table$age[w]  + 1
  pop.table$time.in.current.herd[w] <- pop.table$time.in.current.herd[w] + 1
  return(pop.table)
}
############ Vieillissement ###########################################
module.aging = function(LHerds,param.allHerds){
####################################################################
  n.herds <- length(LHerds)
  LHerds <- lapply(1:n.herds,function(i){module.aging.oneHerd(pop.table = LHerds[[i]],param = param.allHerds[[i]])})
  return(LHerds)
}

#################### select parents in each Herd #########
module.select.parents <- function(LHerds, param.allHerds,sex){
###########################################################
  n.herds <- length(param.allHerds)
  Lparent <- lapply(1:n.herds,function(i){
    L.i <- LHerds[[i]]
    age.min.repro.i = switch(sex,
                             "M" =   param.allHerds[[i]]$age.min.repro.ram,
                             "F" = param.allHerds[[i]]$age.min.repro.ewe)
    age.max.repro.i = switch(sex,
                             "M" =   param.allHerds[[i]]$age.max.repro.ram,
                             "F" = param.allHerds[[i]]$age.max.repro.ewe)
    
    w.i <- which((L.i$herd != -1) & (L.i$sex == sex) & (L.i$age >=  age.min.repro.i) & (L.i$age <=  age.max.repro.i))
    return(L.i[w.i,])}
  )
  return(Lparent)
}


####### Reproduction one herd #############################
module.reproduction.oneHerd = function(mothers,fathers,num.gen,param)
########################################################
  {
  num.herd <- mothers$herd[which(mothers$herd != -1)[1]]
  ######################################################## 
  n.mothers <- nrow(mothers)             # nbre de ewe able to do babies
  possible.father <- fathers$ind
  n.fathers <- nrow(fathers)
  if ((n.fathers >0) & (n.mothers > 0)){
    n.newborns.per.mother <- sample(param$rate.repro[,1],n.mothers,replace = TRUE,prob  = param$rate.repro[,2]) # nb of newborns per ewe.
    n.newborns <- sum(n.newborns.per.mother)
    # creation table newborn
    newborn.table <- as.data.frame(matrix(data = NA, nrow = n.newborns , ncol=ncol(mothers)))
    colnames(newborn.table) <- colnames(mothers)
    # fullfilling newborn.table
    newborn.table$ind <- paste(num.gen,num.herd,1:n.newborns,sep="-")
    newborn.table$age <- 0
    newborn.table$time.in.current.herd <- 0
    newborn.table$sex <- sample(c("F","M"), size = n.newborns, replace = T)
    newborn.table$herd <- num.herd; 
    newborn.table$mother <- rep(mothers$ind,n.newborns.per.mother)
    u <- sample(possible.father,size = n.mothers,replace = T)
    newborn.table$father <- rep(u,n.newborns.per.mother)
  }else{
    newborn.table = NULL
  }
  return(newborn.table)
}

####### Reproduction all herd #############################
module.reproduction <- function(Lmothers, Lfathers,num.gen,param.allHerds){
########################################################
  n.herds <- length(param.allHerds)
  LnewBorns <- lapply(1:n.herds,function(i){module.reproduction.oneHerd(mothers = Lmothers[[i]],fathers = Lfathers[[i]],num.gen,param = param.allHerds[[i]])})
  return(LnewBorns)
}




#################### Exchange fathers ##############
module.exchange.ram.for.one.year <- function(Lfathers,ram.for.repro.Network){
########################################################  
    n.herds <- length(Lfathers)
    newFathers <- Lfathers
    for (i in 1:n.herds){
      prob.i <- ram.for.repro.Network[i, ]
      prob.i <- prob.i/sum(prob.i)
      where.i <-sample(1:n.herds,1,prob = prob.i)
      newFathers[[i]] <- Lfathers[[where.i]]
    }
    Lfathers <- newFathers
    return(Lfathers)
}

#################### select.rams.togive  #################" 
module.select.rams.tomove <- function(LHerds,param.allHerds){
######################################################"

  n.herds <- length(LHerds) 
  L.tomove <- vector(mode = "list", length = n.herds) 
  
  for (i in 1:n.herds){
    pop.table.i <- LHerds[[i]]
    param.i  <- param.allHerds[[i]]
    career.ram.i <- param.i$career.ram
    age.lim.i <- param.i$age.max.repro.ram
    # the ones we give are in the herd,male, still able to reproduct but spent too much time in the current herd
    
    w.i <- which((pop.table.i$herd != -1) 
                 & (pop.table.i$sex == 'M') 
                 & (pop.table.i$age < age.lim.i) 
                 & (pop.table.i$time.in.current.herd >= career.ram.i)
    )
    if(length(w.i) >0){
      L.tomove[[i]] <- pop.table.i[w.i, ]
      LHerds[[i]]$herd[w.i] <- -1 
    }
  }
  res <- list(L.tomove = L.tomove,LHerds = LHerds)
  return(res)
  
}


#################### module.replace.interHerd  #################" 
module.replace.interHerd = function(LHerds,L.togive,ExchangeNetwork,param.allHerds,sex){
########################################################  
  
  # LHerds : composition of all the herds
  # L.togive : list of all animals that are available to be given in all the herds
  # ram.Network : network of exchanges of ram
  # param.allHerds : parameter inside all the herds
  
  
  
  n.herds <- length(LHerds)
  
  order_ex <- sample(1:n.herds,n.herds,replace=FALSE)
  for (i in order_ex){
    
    pop.table.i <- LHerds[[i]]
    param.i <- param.allHerds[[i]]
    #########################"" test
    age.lim.i <- switch(sex,
                       "M" = param.i$age.max.repro.ram,
                       "F" = param.i$age.max.repro.ewe)
    size.i <- switch(sex,
                       "M" = param.i$n.ram,
                       "F" = param.i$n.ewe)
    career.ram.i <- switch(sex,
                           "M" = param.i$career.ram,
                           "F" = Inf)
    
    ### on va chercher les trop vieux
    w.TooOld.i <- which((pop.table.i$herd != -1) 
                        & (pop.table.i$sex == sex) 
                        & (pop.table.i$age >= age.lim.i)
                        )

    n.TooOld.i <- length(w.TooOld.i)
    ### prendre en compte aussi ceux qu'on a enlevé car passé trop de temps (ils sont déjà -1)
    n.Lacking.i <-  size.i - sum((pop.table.i$herd != -1) & (pop.table.i$sex == sex)) + n.TooOld.i
    
    if ( n.Lacking.i > 0){
      
      if (n.TooOld.i > 0) {pop.table.i$herd[w.TooOld.i]<- -1}
      #### chose donnor using ExchangeNetwork. Exchange.network may be weigthed 
      prob.i <- ExchangeNetwork[i, ]
      if (sum(prob.i) > 0){
        prob.i <- prob.i/sum(prob.i)
        donnor.i <-sample(1:n.herds,1,prob = prob.i)
      
        #### select young animals
        L.i <- L.togive[[donnor.i]]
        w.i <- which(L.i$sex == sex)
      
        ### replace too old rams and update LHerds and newborns to give
        if (length(w.i) > 0){
          
          u.i <- sample(w.i,min(length(w.i),n.Lacking.i),replace=FALSE)
          L.i.given <- L.i[u.i,]
          L.togive[[donnor.i]] <- L.i[-u.i,] 
          L.i.given$herd <- i
          L.i.given$time.in.current.herd <- 0
          pop.table.i <- rbind(pop.table.i,L.i.given)
          LHerds[[i]] <- pop.table.i
        }
      }
    }
  }
  return(res = list(LHerds = LHerds, L.togive = L.togive))
}



############# module loose part of herds
module.lose <- function(pop.table,loseProportion = 0.5){
####################################################################
  
  if ((loseProportion > 1) |  (loseProportion < 0)){stop("loseProportion must be between 0 and 1")}
  w_inHerd <- which(pop.table$herd != -1)
  n_Dead <- floor(loseProportion * length(w_inHerd))
  w_Dead <- sample(w_inHerd,n_Dead,replace=FALSE)
  pop.table$herd[w_Dead] = -1
  return(pop.table)
}


################ compute.Inbreeding.Function ###################
compute.inbreeding = function(LHerds){
####################################################################
  
  n.herds <- length(LHerds)
  size.herds <- sapply(1:n.herds,function(i){nrow(LHerds[[i]])})
  
  ### construction de la généalogie
  ped <- do.call("rbind",LHerds)
  ped <- ped[!duplicated(ped$ind),]
  if(!is.factor(ped$ind)){ped$ind <- as.factor(ped$ind)}
  if(!is.factor(ped$father)){ped$father <- as.factor(ped$father)}
  if(!is.factor(ped$mother)){ped$mother <- as.factor(ped$mother)}
  
  wf0 <- which(ped$father=='0')
  LEV = unique(c(levels(ped$ind),levels(ped$father),levels(ped$mother)))
  
  
  ped2 <- ped
  ped2$ind <- as.numeric(factor(ped$ind,levels=LEV))
  ped2$father <- as.numeric(factor(ped$father,levels=LEV))
  ped2$mother <- as.numeric(factor(ped$mother,levels=LEV))
  
  ped2$father[wf0] <- ped2$mother[wf0]  <- 0 
  ped2$sex  <- 2*(ped$sex=='F') + 1*(ped$sex=='M')
  geneal <- gen.genealogy(ped2)
  u <- which(ped2$herd!=-1)
  inBreed <- gen.f(geneal,ped2$ind[u])
  
  U <- as.data.frame(cbind(ped2$herd[u],inBreed))
  names(U) <- c('herd','inBreed')
  U$herd <- as.factor(U$herd)
  return(U)
}

################ compute.Inbreeding.Function ###################
compute.herds.size = function(LHerds){
####################################################################
  return(vapply(LHerds,function(u){sum(u$herd!= - 1)},1))
}