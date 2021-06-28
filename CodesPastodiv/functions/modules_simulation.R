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


####### Reproduction one herd #############################
module.reproduction.oneHerd = function(mothers,fathers,num.gen,param)
########################################################
  {
  
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
    newborn.table$ind <- paste(num.gen,num.herd,1:n.newborns,sep="-")
    newborn.table$age <- 0
    newborn.table$sex <- sample(c("F","M"), size = n.newborns, replace = T)
    newborn.table$herd <- num.herd; 
    newborn.table$mother <- rep(mothers$ind,n.newborns.per.mother)
    newborn.table$father <- rep(sample(possible.father,size = n.mothers,replace = T),n.newborns.per.mother)
    newborn.table$time.in.current.herd <- 0
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


#################### select parents in each Herd #########
module.select.parents <- function(LHerds, param.allHerds,sex){
########################################################
  n.herds <- length(param.allHerds)
  Lparents <- lapply(1:n.herds,function(i){
    L.i <- LHerds[[i]]
    age.min.repro.i = switch(sex,
        "M" =   param.allHerds[[i]]$age.min.repro.ram,
        "F" = param.allHerds[[i]]$age.min.repro.ewe)
    max.time.in.herd.i = switch(sex,
          "M" =   param.allHerds[[i]]$max.time.in.herd,
          "F" = Inf)
    
    w.i <- which((L.i$herd != -1) & (L.i$sex == sex) & (L.i$age >=  age.min.repro.i) & (L.i$time.in.herd <= max.time.in.herd.i))
    return(L.i[w.i,])}
  )
  return(Lparents)
}

#################### Exchange fathers ##############
module.exchange.ram <- function(Lfathers,ram.for.repro.Network){
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



#################### module.replace.interHerd  #################" 
module.replace.interHerd = function(LHerds,Lnewborns.togive,ExchangeNetwork,param.allHerds,sex){
########################################################  
  
  # LHerds : composition of all the herds
  # Lnewborns.togive : list of all newborn that are available to be given in all the herds
  # ram.Network : network of exchanges of ram
  # param.allHerds : parameter inside all the herds
  
  
  #if(sex=='F'){browser()}
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
    max.time.in.herd.i <- switch(sex,
                                 "M" = param.i$max.time.in.herd,
                                 "F" = Inf)
    w.TooOld.i <- which((pop.table.i$herd != -1) 
                        & (pop.table.i$sex == sex) 
                        & (pop.table.i$age >= age.lim.i)
                        & (pop.table.i$time.in.herd >= max.time.in.herd.i ))
    
    n.TooOld.i <- length(w.TooOld.i)
    n.Lacking.i <-  size.i - sum((pop.table.i$herd != -1) & (pop.table.i$sex == sex)) + n.TooOld.i
    #if(n.Lacking.i < 0){browser()}
    
    if ( n.Lacking.i > 0){
      
      if (n.TooOld.i > 0) {pop.table.i$herd[w.TooOld.i]<- -1}
      
      #### chose donnor using ExchangeNetwork. Exchange.network may be weigthed 
      prob.i <- ExchangeNetwork[i, ]
      prob.i <- prob.i/sum(prob.i)
      donnor.i <-sample(1:n.herds,1,prob = prob.i)
      
      #### select young rams
      L.i <- Lnewborns.togive[[donnor.i]]
      w.i <- which(L.i$sex == sex)
      
      ### replace too old rams and update LHerds and newborns to give
      if (length(w.i) > 0){
        u <- sample(1:length(w.i),min(length(w.i),n.Lacking.i),replace=FALSE)
        L.i.given <- L.i[w.i[u],]
        Lnewborns.togive[[donnor.i]] <- L.i[-w.i[u],] 
        L.i.given$herd <- i
        L.i.given$time.in.herd <- 0
        pop.table.i <- rbind(pop.table.i,L.i.given)
        LHerds[[i]] <- pop.table.i
      }
    }
  }
  return(res = list(LHerds = LHerds, Lnewborns.togive = Lnewborns.togive))
}

#################### select.rams.togive  #################" 
module.select.rams.togive <- function(LHerds,param.allHerds){
######################################################"
  n.herds <- length(LHerds) 
  L.ramtogive <- vector(mode = "list", length = n.herds) 
 
 
  for (i in 1:n.herds){
    L.i <- LHerds[[i]]
    param.i  <- param.allHerds[[i]]
    w.i <- which((L.i$herd != -1) 
                 & (L.i$sex == 'M') 
                 & (L.i$time.in.herd > param.i$career.ram)
                 & (L.i$age < param.i$age.max.repro.ram)
    )
    if(length(w.i) >0){
      L.ramtogive[[i]] <- L.i[w.i, ]
      LHerds[[i]]$herd[w.i] <- -1 
    }
  }
  res <- list(L.ramtogive = L.ramtogive,LHerds = LHerds)
  return(res)
    
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
  

  code.no.parents <- as.character(ped$father[which(ped$father=='0')][1])
  
  LEV = unique(c(levels(ped$ind),levels(ped$father),levels(ped$mother)))
  ped$ind=as.numeric(factor(ped$ind,levels=LEV))
  ped$father=as.numeric(factor(ped$father,levels=LEV))
  ped$mother=as.numeric(factor(ped$mother,levels=LEV))
  ped$father[ped$father == ped$father[1]]=0
  ped$mother[ped$mother == ped$mother[1]]=0
  
  
  num_sex<- rep(2,length(ped$sex))
  num_sex[ped$sex=='M'] <- 1
  ped$sex = num_sex
  
  
 
  geneal <- gen.genealogy(ped)
  u <- which(ped$herd!=-1)
  inBreed <- gen.f(geneal,u)
  
  U <- as.data.frame(cbind(ped$herd[u],inBreed))
  names(U) <- c('herd','inBreed')
  
  return(U)
}

################ compute.Inbreeding.Function ###################
compute.herds.size = function(LHerds){
####################################################################
  return(vapply(LHerds,function(u){sum(u$herd!= - 1)},1))
}