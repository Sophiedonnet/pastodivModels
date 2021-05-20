library(plyr)
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
  pop.table$age <- pop.table$age  + 1
  return(pop.table)
}
############ Vieillissement ###########################################
module.aging = function(LHerds,param.allHerds){
####################################################################
  n.herds <- length(LHerds)
  LHerds <- lapply(1:n.herds,function(i){module.aging.oneHerd(pop.table = LHerds[[i]],param = param.allHerds[[i]])})
  return(LHerds)
}




####################################################################
# reproduction from mothers and fathers
####################################################################
module.reproduction.oneHerd = function(mothers,fathers,num.gen,param)
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
    newborn.table$ind <- paste(num.gen,num.herd,1:n.newborns,sep="-")
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

module.reproduction <- function(Lmothers, Lfathers,num.gen,param.allHerds){
  
  n.herds <- length(param.allHerds)
  LnewBorns <- lapply(1:n.herds,function(i){module.reproduction.oneHerd(mothers = Lmothers[[i]],fathers = Lfathers[[i]],num.gen,param = param.allHerds[[i]])})
  return(LnewBorns)
}

########################################################
#################### select fathers in each Herd
########################################################
module.select.parents <- function(LHerds, param.allHerds,sex){
  
  n.herds <- length(param.allHerds)
  Lparents <- lapply(1:n.herds,function(i){
    L.i <- LHerds[[i]]
    age.lim.i = switch(sex,
        "M" =   param.allHerds[[i]]$age.repro.ram,
        "F" = param.allHerds[[i]]$age.repro.ewe)
    w.i <- which((L.i$herd != -1) & (L.i$sex == sex)&(L.i$age >=  age.lim.i))
    return(L.i[w.i,])}
  )
  return(Lparents)
}

########################################################
#################### Exhange fathers
########################################################
module.exchange.ram <- function(Lfathers,ram.for.repro.Network){
  
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
  
    #where.fathers <- module.exchange.ram(ram.for.repro.Network,n.herds)
    #Lfathers<- lapply(1:n.herds,function(i){do.call("rbind",Lfathers[where.fathers[[i]]])})
}



module.replace.interHerd = function(LHerds,Lnewborns.togive,ExchangeNetwork,param.allHerds,sex){
  
  
  # LHerds : composition of all the herds
  # Lnewborns.togive : list of all newborn that are available to be given in all the herds
  # ram.Network : network of exchanges of ram
  # param.allHerds : parameter inside all the herds
  
  n.herds <- length(LHerds)
  for (i in 1:n.herds){
    
    pop.table.i <- LHerds[[i]]
    param.i <- param.allHerds[[i]]
    #########################"" test
    age.lim.i <- switch(sex,
                       "M" = param.i$career.ram,
                       "F" = param.i$career.ewe)
    size.i <- switch(sex,
                       "M" = param.i$n.ram,
                       "F" = param.i$n.ewe)
    
    w.TooOld.i <- which((pop.table.i$herd != -1) & (pop.table.i$sex == sex) & (pop.table.i$age >= param.i$career.ram))
    n.TooOld.i <- length(w.TooOld.i)
    n.Lacking.i <-  size.i- sum((pop.table.i$herd != -1) & (pop.table.i$sex == sex))
    
    if (n.TooOld.i > 0){
      pop.table.i$herd[w.TooOld.i]<- -1
      
      #### chose donnor using ram.network. ram.network may be weigthed 
      prob.i <- ExchangeNetwork[i, ]
      prob.i <- prob.i/sum(prob.i)
      donnor.i <-sample(1:n.herds,1,prob = prob.i)
      
      #### select young rams
      L.i <- Lnewborns.togive[[donnor.i]]
      w.i <- which(L.i$sex == sex)
      
      ### replace too old rams and update LHerds and newborns to give
      if (length(w.i) > 0){
        u <- sample(1:length(w.i),min(length(w.i),n.TooOld.i + n.Lacking.i),replace=FALSE)
        L.i.given <- L.i[w.i[u],]
        Lnewborns.togive[[donnor.i]] <- L.i[-w.i[u],] 
        L.i.given$herd <- i
        pop.table.i <- rbind(pop.table.i,L.i.given)
        LHerds[[i]] <- pop.table.i
      }
    }
  }
  return(res = list(LHerds = LHerds, Lnewborns.togive = Lnewborns.togive))
}


# ####################################################################
# # replacement intra troupeau
# ####################################################################
# module.replaceEwe.intraHerd = function(pop.table,newborn.table,param=list()){
#   
# 
#   w.F <- which(newborn.table$sex == 'F')
#   w.TooOld <- which(pop.table$herd != -1 & pop.table$sex == 'F' & pop.table$age >= param$career.ewe)
#   n.TooOld <- length(w.TooOld)
#   n.Lacking <- param$n.ewe - sum((pop.table$herd != -1) & (pop.table$sex == 'F'))
#   newborn.togive <- newborn.table
#   if (n.TooOld > 0){
#     pop.table$herd[w.TooOld]<- -1
#     if (length(w.F)>0){
#           u <- sample(1:length(w.F),min(n.TooOld + n.Lacking,length(w.F)),replace=FALSE)
#           pop.table <- rbind(pop.table,newborn.table[w.F[u],])
#           newborn.togive <- newborn.togive[-w.F[u],]
#     }
#   }
#    
#   
#   res <- list(pop.table = pop.table,newborn.togive  = newborn.togive)
#   return(res)
# }
####################################################################
# replacement of ram inter herds
####################################################################
# module.replaceRam.interHerd = function(LHerds,Lnewborns.togive ,ram.Network,param.allHerds){
# 
#   
#   # LHerds : composition of all the herds
#   # Lnewborns.togive : list of all newborn that are available to be given in all the herds
#   # ram.Network : network of exchanges of ram
#   # param.allHerds : parameter inside all the herds
# 
#   n.herds <- length(LHerds)
#     for (i in 1:n.herds){
#   
#     pop.table.i <- LHerds[[i]]
#     param.i <- param.allHerds[[i]]
#     #########################"" test
#     w.R.TooOld <- which((pop.table.i$herd != -1) & (pop.table.i$sex=='M') & (pop.table.i$age >= param.i$career.ram))
#     n.R.TooOld <- length(w.R.TooOld)
#     if (n.R.TooOld > 0){
#       pop.table.i$herd[w.R.TooOld]<- -1
#     
#       #### chose donnor using ram.network. ram.network may be weigthed 
#       p.i.r <- ram.Network[i, ]
#       p.i.r <- p.i.r/sum(p.i.r)
#       donnor.i <-sample(1:n.herds,1,prob = p.i.r)
#       
#       #### select young rams
#       L.i <- Lnewborns.togive[[donnor.i]]
#       w.i <- which(L.i$sex == 'M')
#       
#       ### replace too old rams and update LHerds and newborns to give
#       if (length(w.i) > 0){
#         u <- sample(1:length(w.i),min(length(w.i),n.R.TooOld),replace=FALSE)
#         L.i.given <- L.i[w.i[u],]
#         Lnewborns.togive[[donnor.i]] <- L.i[-w.i[u],] 
#         L.i.given$herd <- i
#         pop.table.i <- rbind(pop.table.i,L.i.given)
#         LHerds[[i]] <- pop.table.i
#       }
#     }
#     }
#     return(res = list(LHerds = LHerds, Lnewborns.togive = Lnewborns.togive))
# }


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
  ped <- ped[!duplicated(ped$ind),]
  
  LEV <- unique(c(levels(allHerds$ind), levels (allHerds$father), levels(allHerds$mother)))
  indNew <- mapvalues(ped$ind, from = LEV, to = 1:length(LEV))
    
  #LEV <- unique(c(( as.character(allHerds$ind)),as.character( allHerds$father),as.character( allHerds$mother)))
  ped$ind<- as.numeric(factor(ped$ind,levels=LEV))
  ped$father <- as.numeric(factor(ped$father,levels=LEV))
  ped$mother <- as.numeric(factor(ped$mother,levels=LEV))
  
  
  w.without.knwon.parents  <-  which(ped$father==which(LEV == "0"))
  ped$father[w.without.knwon.parents] <- 0
  ped$mother[w.without.knwon.parents] <- 0
  
  
  ped$sex <- as.numeric(mapvalues(ped$sex, from = c("M", "F"), to = c(1, 2)))
  
 
  geneal <- gen.genealogy(ped)
  u <- which(ped$herd!=-1)
  inBreed <- gen.f(geneal,u)
  
  U <- as.data.frame(cbind(ped$herd[u],inBreed))
  names(U) <- c('herd','inBreed')
  
  return(U)
}




