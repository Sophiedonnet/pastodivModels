library(GENLIB)
library(tidyverse)
############################################# ESSAI

Simulate.herds = function(n.herds,n.generations,param.allHerds=NULL,herds.Network = list(),LHerds=NULL,computeInbreeding = FALSE)
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
  
  ######################
  if(is.null(herds.Network$ram.for.repro)){herds.Network$ram.for.repro = diag(n.herds)}
  if(is.null(herds.Network$ram.for.replace)){herds.Network$ram.for.replace = diag(n.herds)}
  if(is.null(herds.Network$ewe.for.replace)){herds.Network$ewe.for.replace = diag(n.herds)}
  
  ram.for.repro.Network <- herds.Network$ram.for.repro
  ram.for.replace.Network <- herds.Network$ram.for.replace
  ewe.for.replace.Network <- herds.Network$ewe.for.replace
  
  if (sum(sapply(herds.Network,function(E){dim(E)}) != n.herds) > 0) {
    stop('At least one network is not of the correct dimension')
  }
  
  #####################
  
  ################### initialisation   = génération 0 or continue from LHerds given
  if (is.null(LHerds)){
    LHerds <- module.initialize(param.allHerds)
    gen <- gen0 <- 0
    }else{
    gen <- gen0 <-  max(as.numeric(unique(str_split_fixed(LHerds[[1]]$ind,"-",3)[,1]))) 
    n.generations <- gen + n.generations
  }
  
  herds_size <- matrix(compute.herds.size(LHerds),nrow = 1)
  print(herds_size)
  if(computeInbreeding){
    inBreeding <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(inBreeding) <- c("herd","inBreed", "gen")
  }
  ###################  GENERATIONS
  while(gen < n.generations){ 
    gen <- gen + 1
    print(gen)
    #-------------------  tout le monde prend 1 an
    LHerds <- module.aging(LHerds,param.allHerds)
    
    #------------------ reproduction with possibly Ram exchanges 
    
    #---- select mothers
    Lmothers <- module.select.parents(LHerds, param.allHerds,sex = 'F')
    #--- select fathers
    Lfathers <- module.select.parents(LHerds, param.allHerds,sex = 'M')
    
    #--- exchange fathers
    
    Lfathers <- module.exchange.ram(Lfathers,ram.for.repro.Network)
   
    #---------- reproduction
    
    Lnewborns <- module.reproduction(Lmothers,Lfathers,num.gen = gen,param.allHerds)
    
    #---------- replace old ewe by newborns
    resultsReplaceEwe  <-  module.replace.interHerd(LHerds,Lnewborns ,ExchangeNetwork = ewe.for.replace.Network,param.allHerds,sex = 'F')
    LHerds <- resultsReplaceEwe$LHerds
    Lnewborns.togive <- resultsReplaceEwe$Lnewborns.togive
    
   
    #--------- replace ram inter Herds
    resultsReplaceRam  <-  module.replace.interHerd(LHerds,Lnewborns.togive ,ExchangeNetwork = ram.for.replace.Network,param.allHerds,sex = 'M')
    LHerds <- resultsReplaceRam$LHerds
    
    #------------- size of herds 
    herds_size <- rbind(herds_size,compute.herds.size(LHerds))
    # ----- inBreeding                    
                        
   if(computeInbreeding){
       inBreeding_gen <- compute.inbreeding(LHerds)
       inBreeding_gen$gen <- gen
       inBreeding <- rbind(inBreeding, inBreeding_gen)
      }
  }
  if(computeInbreeding){
    inBreeding$herd <- as.factor(inBreeding$herd)
  }
  
  
  herds_size <- as.data.frame(herds_size)
  colnames(herds_size) <- 1:n.herds
  herds_size$gen <- gen0:gen
  herds_size <- herds_size %>% gather(herd,size,  -c(gen)) 
  herds_size$herd <- as.factor(herds_size$herd)
  res = list(LHerds = LHerds,herds_size = herds_size)
  if (computeInbreeding){res$inBreeding = inBreeding}
  return(res)
}