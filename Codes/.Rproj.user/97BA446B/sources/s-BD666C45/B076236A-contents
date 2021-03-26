library(GENLIB)
############################################# ESSAI

Simulate.herds = function(n.herds,n.generations,param.allHerds=NULL,herds.Network,LHerds=NULL,computeInbreeding = FALSE)
{
  
  if(is.null(param.allHerds)){
    param.default <- list(n.ram = 2,
                          n.ewe = 40,
                          career.ram = 8,
                          career.ewe = 8,
                          age.min.ram = 0,
                          age.min.ewe = 0,
                          age.repro.ewe = 3,
                          age.repro.ram = 1
    )
    param.default$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
    names(param.default$rate.repro) = c('nb.lambs','probability')
    param.allHerds <- lapply(1:n.herds,function(i){param.default})
  }
  
  ################### initialisation   = génération 0 
  if (is.null(LHerds)){
    LHerds <- lapply(1:n.herds,function(i){module.initialize.oneHerd(i,param = param.allHerds[[i]])})
    gen <-  0
  }else{
    gen <- max(as.numeric(unique(str_split_fixed(LHerds[[1]]$ind,"-",3)[,1]))) 
  }
  ###################  GENERATIONS
  while(gen <= n.generations){ 
    gen <- gen + 1
    
    #-------------------  tout le monde prend 1 an
    LHerds <- lapply(1:n.herds,function(i){module.aging.oneHerd(pop.table = LHerds[[i]],param = param.allHerds[[i]])})
    
    #------------------ reproduction with possibly Ram exchanges 
    #--- select fathers
    fathers <- lapply(1:n.herds,function(i){
      L.i <- LHerds[[i]]
      w.i <- which((L.i$herd != -1) & (L.i$sex == 'M')&(L.i$age >= param.allHerds[[i]]$age.repro.ram))
      return(L.i[w.i,])}
    )
    #--- exchange fathers
    where.fathers <- choose.Ram.reproduction(herds.Network,n.herds)
    fathers<- lapply(1:n.herds,function(i){do.call("rbind",fathers[where.fathers[[i]]])})
    
    #---- select mothers
    mothers <- lapply(1:n.herds,function(i){
      L.i <- LHerds[[i]]
      w.i <- which((L.i$herd != -1) & (L.i$sex == 'F') & (L.i$age >= param.allHerds[[i]]$age.repro.ewe))
      return(L.i[w.i,])}
    )
    
    #---------- reproduction
    #browser()
    if (inherits(try(lapply(1:n.herds,function(i){module.reproduction(mothers = mothers[[i]],fathers = fathers[[i]],gen,param = param.allHerds[[i]])})),"try-error")) browser()
    newBorns <- lapply(1:n.herds,function(i){module.reproduction(mothers = mothers[[i]],fathers = fathers[[i]],gen,param = param.allHerds[[i]])})
    
    
    #---------- replace old ewe by newborns
    if (inherits(try(lapply(1:n.herds,function(i){module.replaceEwe.intraHerd(pop.table = LHerds[[i]],newborn.table = newBorns[[i]],param = param.allHerds[[i]])})),"try-error")) browser()
    resultsReplace <- lapply(1:n.herds,function(i){module.replaceEwe.intraHerd(pop.table = LHerds[[i]],newborn.table = newBorns[[i]],param = param.allHerds[[i]])})
    
    LHerds <- lapply(1:n.herds,function(i){resultsReplace[[i]]$pop.table})
    
    #-------- newborns to give to other herds
    Lnewborns.togive <- lapply(1:n.herds,function(i){resultsReplace[[i]]$newborn.togive})

    for (i in 1:n.herds){
      pop.table.i <- LHerds[[i]]
      param.i <- param.allHerds[[i]]
      w.R.TooOld <- which((pop.table.i$herd != -1) & (pop.table.i$sex=='M') & (pop.table.i$age >= param.i$career.ram))
      n.R.TooOld <- length(w.R.TooOld)
      if (n.R.TooOld > 0){
      
        donnor.i <-sample(which(herds.Network[i,]==1),1)
        L.i <- Lnewborns.togive[[donnor.i]]

        if(!is.null(L.i)){
          w.i <- which(L.i$sex == 'M')
          u <- sample(w.i,min(length(w.i),n.R.TooOld),replace=FALSE)
          L.i.given <- L.i[u,]
          Lnewborns.togive[[donnor.i]] <- L.i[-u,] 
          L.i.given$ herd <- i
          pop.table.i$herd[w.R.TooOld]<- -1
          pop.table.i <- rbind(pop.table.i,L.i.given)
          LHerds[[i]] <- pop.table.i
          nbRam.i <-sum((LHerds[[i]]$herd != -1) & (LHerds[[i]]$sex== 'M'))
          print(c(i,nbRam.i))
          if(nbRam.i> param.allHerds[[i]]$n.ram ){browser()}
      
        }
      }
    }
    
    nbRam <- sapply(1:n.herds,function(i){sum(LHerds[[i]]$herd != -1 & LHerds[[i]]$sex== 'M')})
    nbEwe <-sapply(1:n.herds,function(i){sum(LHerds[[i]]$herd != -1 & LHerds[[i]]$sex== 'F')})
    print(param.allHerds[[1]]$n.ram)
   
    if(any(nbRam>param.allHerds[[1]]$n.ram)){browser()}
    print(nbRam)
    print(nbEwe)

    #if(computeInbreeding){inBreeding <- computeInbreedingFunction(LHerds)}
    
    
  }
  return(LHerds)
}
