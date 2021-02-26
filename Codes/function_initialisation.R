initialize_herd <- function(num_herd,param=list()){
  # num_herd num to identify the herd
  # param list of useful parameters for initialising the herd
  #----------------------------------------------
  param_default <- list(n.ram = 2,
                        n.ewe = 40,
                        career.ram = 8,
                        age.min.ram=0,
                        career.ewe = 8,
                        age.min.ewe = 0
  )
  
  param_default$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
  names(param_default$rate.repro) = c('nb.lambds','probability')
  
  #---------------------------------------------
  param_default[names(param)] <- param
  
  pop.table = data.frame(id=paste("0",num_herd,1:(param_default$n.ewe+param_default$n.ram),sep="."),
                         father = "0",mother = "0",
                         herd=num_herd,sex=c(rep("M",param_default$n.ram),rep("F",param_default$n.ewe)),
                         age = c(sample(param_default$age.min.ram:(param_default$career.ram-1),replace = T,param_default$n.ram),
                                 sample(param_default$age.min.ewe:(param_default$career.ewe-1),replace=T,param_default$n.ewe)))
# - 1 car on les fait vieillir d'un an, que pour ram ? pour etre sur d'avoir des reproducteurs en age
  return(pop.table)
}

#initialize_herd(num_herd=2)

# reproduction intra troupeau
module_reproduction_intra(Lherds,param=list())
{
  param_default = list(age.repro.ewe=1,age.repro.ram=0) 
  param_default$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
  names(param_default$rate.repro) = c('nb.lambs','probability')
  
  Lherds = lapply(Lherds,function(herd)
    {
    # reproduction par herd au sein des troupeaux, ajouter les enfant au herd actuel
    offspring 
    return(rbind(herd,offspring))
  })
  
  
}



# fonction generale
Simulate_herds = function()
{
  
  LHerds = lapply(1:nherds,initialize_herd())
  
  for (t in 1:nherds)
  {
    
    module_aging()
    module_reproduction()
    module_repartition()
    save(LHerds)
    
  }
  
}