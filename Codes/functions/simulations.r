simulInitGeneration <- function(param=list()){
  
  
  
  #----------------------------------------------
  param_default <- list(n.ram = 2,
                        n.ewe = 40, 
                        life.expec.ewe = 8, 
                        age.repro.ewe = 1, 
                        career.ram = 8 
                        )
  
  param_default$rate.repro = as.data.frame(cbind(c(0,1,2),c(0,1,0)))
  names(param_default$rate.repro) = c('nb.lambds','probability')
  
  #---------------------------------------------
  param_default[names(param)] <- param
  #--------------------------------------------- 
  n.ram <- param_default$n.ram
  n.ewe <- param_default$n.ewe
  life.expec.ewe <- param_default$life.expec.ewe
  age.repro.ewe <- param_default$age.repro.ewe
  career.ram <- param_default$career.ram
  successrate.repro <- param_default$successrate.repro
  n.pop <- n.ram + n.ewe
  
  #--- NAMES --------------------------------
  u <- c(t(outer(letters,letters,paste0)))
  listNames <- c(t(outer(u,0:9,paste0)))  
  while (length(listNames)<n.pop){
    listNames = c(outer(letters,listNames,paste0))
  }  
  
  
  #----------------------------- 
  a <- 0 # generation id
  
  
  colnames(pop.table) <- c("id","father","mother","age","sexe", "career")
  # set sexe values
  pop.table$sexe[1:n.ram] <- "M"
  pop.table$sexe[(n.ram+1):n.pop] <-"F"
  # set id 
  pop.table$id <- paste0(a,listNames[1:n.pop])
  
  # set age classes
  pop.table$age[(n.ram+1):n.pop] <- age.seed.ewe   # set age for ewes
  pop.table$age[1:n.ram] <- age.seed.ram
  # set initial career time for rams & ewes (default = 0)
  pop.table[,6] <- 0
  # set father and mother id (random for first run):
  pop.table[,2] <- '0.0'
  pop.table[,3] <-

# check the Gen 0 dataframe :
#View(pop.table)
# check ewe age distribution in gen 0:
table(pop.table$age[pop.table$sexe=="F"])
# check ram age distribution in gen 0:
table(pop.table$age[pop.table$sexe=="M"])


}
