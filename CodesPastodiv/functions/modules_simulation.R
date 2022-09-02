library(plyr)
library(tidyverse)

######### Initialize one herd  ########################
module.initialize.oneHerd <- function(num.herd,param,seed=NULL,init.genes = NULL){
####################################################################
# num.herd num to identify the herd
# param list of useful parameters for initialising the herd
#----------------------------------------------
  if(!is.null(seed)){set.seed(seed)}
  
  if(is.null(param$age.min.ram)){param$age.min.ram = 0}
  if(is.null(param$age.min.ewe)){param$age.min.ewe = 0}
  if(is.null(param$nbgenes.coding)){param$nbgenes.coding  = 10}
  if(is.null(param$nbgenes.noncoding)){param$nbgenes.noncoding  = 10}
  
  
  pop.table = data.frame(ind=paste("0",num.herd,1:(param$n.ewe+param$n.ram),sep="-"),
                         father = "0",mother = "0",
                         herd=num.herd,
                         sex=c(rep("M",param$n.ram),rep("F",param$n.ewe)),
                         age = c(sample(param$age.min.ram:(param$age.max.repro.ram-1),replace = T,param$n.ram),
                                 sample(param$age.min.ewe:(param$age.max.repro.ewe-1),replace=T,param$n.ewe)))
  
  pop.table$time.in.current.herd <- 0; 

  
  ### initialisaiton des genes 
  nb.genes <- param$nbgenes.coding + param$nbgenes.noncoding
  if (nb.genes>0){
    nb.ind <-   param$n.ewe+param$n.ram
    Genes  <- matrix(sample(c(0,1,2), nb.ind * nb.genes,replace= TRUE, prob = c(1,1,1)), nb.ind, nb.genes)
    namesCols <- c()
    if (param$nbgenes.coding>0){
      namesCols <- c(namesCols, paste('codingGene', 1: param$nbgenes.coding, sep='.'))
      }
    if (param$nbgenes.noncoding>0){
      namesCols <- c(namesCols, paste('noncodingGene', 1: param$nbgenes.noncoding, sep='.'))
    }
    colnames(Genes)<- namesCols
    pop.table <- cbind(pop.table, Genes)
  }
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
    u <- sample(possible.father,size = n.mothers, prob = compute.val.repro(fathers),replace = T)
    newborn.table$father <- rep(u,n.newborns.per.mother)
    
    
    genes.places.1 <- which(substr(names( newborn.table),1,6)=='coding')[1] # indice de la premiere colonne des genes 
    genes.places <- genes.places.1:ncol(mothers)
    browser()
    for(i in 1:n.newborns){
     
      newborn.i <- newborn.table[i, ]
      genes.mother <- mothers[mothers$ind==newborn.i$mother, genes.places]
      genes.father <- fathers[fathers$ind==newborn.i$father, genes.places]
      genes.newborn.i <- simulation.genes.one.newborn(genes.mother,genes.father)
      newborn.table[i, genes.places] <- genes.newborn.i
    }
    
    ### attribute genes to newborns; 
    
  }else{
    newborn.table = NULL
  }
  return(newborn.table)
}
############ Calcul de la valeur reproductive des males à partir de leurs genes codant
compute.val.repro = function(pop.father,selectionWeights=NULL){
  
  
  fatherColNames <- names(pop.father)
  coding <- which(substr(fatherColNames,1,6)=='coding')
  if(is.null(selectionWeights)){
    selectionWeights = seq(length(coding),1)
  }
  codingGenes <- pop.father[,coding]
  #browser()
  
  return(c(as.matrix(codingGenes) %*%selectionWeights))
}

############ Calcul des genes des enfants à partir des genes des peres et meres. 
simulation.genes.one.newborn = function(genes.mother,genes.father){
  
  numeric.genes <- 3*genes.mother + genes.father + 1; 
  
  myHereditaryMatrix <- matrix(0,9,3);
  myHereditaryMatrix[1,1] <- 1;     # mother  0  father 0 
  myHereditaryMatrix[2,1:2] <- 1/2; # mother 0  father 1
  myHereditaryMatrix[3,2] <- 1;     # mother 0  father 2
  myHereditaryMatrix[4,1:2] <- c(1/2,1/2);   # mother 1  father 0
  myHereditaryMatrix[5,] <- c(1/4,1/2,1/4);   # mother 1  father 1
  myHereditaryMatrix[6,2:3] <- c(1/2,1/2);   # mother 1  father 2
  myHereditaryMatrix[7,2] <- 1 ;   # mother 2  father 0
  myHereditaryMatrix[8,2:3] <- c(1/2,1/2) ;   # mother 2  father 1
  myHereditaryMatrix[9,3] <- 1 ;   # mother 2  father 2
  
  return(sapply(1:length(numeric.genes),function(g){sample(c(0,1,2),1, myHereditaryMatrix[numeric.genes[g],],replace = TRUE)}))
  
  
  
  
  
  
  
  
  
  
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
compute.geneal = function(LHerds){
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
  res <- list(geneal=geneal, ped = ped2)
  return(res)
}
################ compute.Inbreeding.Function ###################
compute.inbreeding = function(geneal, ped){
####################################################################
  
# geneal : result of "gen.genealogy"
# ped : pedegree of all the individuals
  
  u <- which(ped$herd!=-1)
  inBreed <- gen.f(geneal,ped$ind[u])
  U <- as.data.frame(cbind(ped$herd[u],inBreed))
  names(U) <- c('herd','inBreed')
  U$herd <- as.factor(U$herd)
  return(U)
}

################ compute.Herd Sized .Function ###################
compute.herds.size = function(LHerds){
####################################################################
  return(vapply(LHerds,function(u){sum(u$herd!= - 1)},1))
}



################ compute.Kinship.Function ###################
compute.kinship = function(geneal, ped){
####################################################################
  
  #------------------
  u <- which(ped$herd!=-1)
  kinship<-gen.phi(geneal,ped$ind[u], MT = TRUE)
  res <- list(kinshipMatrix  = kinship, herd = ped$herd[u])
  return(res)
}

