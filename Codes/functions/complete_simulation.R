library(GENLIB)
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
  
  #####################
  
  
  ################### initialisation   = génération 0 or continue from LHerds given
  if (is.null(LHerds)){
    LHerds <- module.initialize(param.allHerds)
    gen <-  0
  }else{
    gen <- max(as.numeric(unique(str_split_fixed(LHerds[[1]]$ind,"-",3)[,1]))) 
  }
  
  if(computeInbreeding){
    inBreeding <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(inBreeding) <- c("herd","inBreed", "gen")
  }
  ###################  GENERATIONS
  while(gen <= n.generations){ 
    gen <- gen + 1
    
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
    
    
   if(computeInbreeding){
       inBreeding_gen <- computeInbreedingFunction(LHerds)
       inBreeding_gen$gen <- gen
       inBreeding <- rbind(inBreeding, inBreeding_gen)
     }
  }
  res = list(LHerds = LHerds)
  if (computeInbreeding){res$inBreeding = inBreeding}
  return(res)
}
