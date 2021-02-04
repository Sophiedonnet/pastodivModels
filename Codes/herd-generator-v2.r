#### avant de commencer : sélectionner répertoire de travail

setwd("~/Documents/.../r-consanguinity")

######################
# variables à déclarer :
n.ram <- 2                           # ram number
n.ewe <- 40                         # ewe number
n.gen <- 8                         # number of generations to be simulated
#life.expec.ewe <-8
age.repro.ewe <- 1                   # age at which an ewe will give birth / lower limit (only > will reproduce)
career.ram <-8                     # time a ram will be kept in the herd
#successrate.repro <-  1            # % of chance of giving birth for an ewe : to be added
  
# variables internes au modèle :
n.pop <- n.ram + n.ewe                                              # size of the total flock (including rams)
age.seed.ewe <-  sample(c(1:8), size = n.ewe, replace = T)          # random distribution of ages of ewes
age.seed.ram <- sample(c(2:8), size = n.ram, replace = T)           # set ages of rams / can be a sample

#########################################################################
############# Step 1 : create Generation 0 population  ##################
#########################################################################

a <- 0 # generation id
b <- sample(letters, size = n.pop, replace = T)
c <- sample(letters, size = n.pop, replace = T)
d <- sample(c(0:9), size = n.pop, replace = T)
id <- paste(a,b,c,d, sep="")

pop.table <- as.data.frame(matrix(data = NA, nrow= n.pop, ncol=6))
colnames(pop.table) <- c("id","father","mother","age","sexe", "career")
# set sexe values
pop.table[c(1:n.ram),5] <- "M"
pop.table[c((n.ram+1):n.pop),5] <-"F"
# set id 
pop.table[,1] <- id
# set age classes
pop.table[c((n.ram+1):n.pop),4] <- age.seed.ewe   # set age for ewes
pop.table[c(1:n.ram),4] <- age.seed.ram
# set initial career time for rams & ewes (default = 0)
pop.table[,6] <- 0
# set father and mother id (random for first run):
a <- "F" 
b <- sample(letters, size = n.pop, replace = T)
c <- sample(letters, size = n.pop, replace = T)
d <- sample(c(0:9), size = n.pop, replace = T)
idf <- paste(a,b,c,d, sep="")
a <- "M" 
b <- sample(letters, size = n.pop, replace = T)
c <- sample(letters, size = n.pop, replace = T)
d <- sample(c(0:9), size = n.pop, replace = T)
idm <- paste(a,b,c,d, sep="")
pop.table[,2] <- idf
pop.table[,3] <- idm

# check the Gen 0 dataframe :
#View(pop.table)
# check ewe age distribution in gen 0:
table(pop.table$age[pop.table$sexe=="F"])
# check ram age distribution in gen 0:
table(pop.table$age[pop.table$sexe=="M"])

###################################################################
############ Step 2 : Loop to simulate flock reproduction #########
###################################################################

for (i in 1:n.gen){
  # extraction brebis se reproduisant
  repro.ewe.table <- pop.table[pop.table$sexe=="F" & pop.table$age > age.repro.ewe,]
  n.repro.ewe <- nrow(repro.ewe.table)             # nbre de brebis en repro
  # creation table newborn
  newborn.table <- as.data.frame(matrix(data = NA, nrow= n.repro.ewe, ncol=6))
  colnames(newborn.table) <- c("id","father","mother","age","sexe", "career")
  n.newborn <- nrow(newborn.table)
  # creation et attribution new id:
  a <- i   # cf n.gen
  b <- sample(letters, size = n.newborn, replace = T)
  c <- sample(letters, size = n.newborn, replace = T)
  d <- sample(c(0:9), size = n.newborn, replace = T)
  newborn.id <- paste(a,b,c,d, sep="")
  newborn.table[,1] <- newborn.id
  newborn.table[, 4] <- 0
  newborn.table[, 5] <- sample(c("F","M"), size = n.repro.ewe, replace = T)
  newborn.table[, 6] <- 0
  # extraction id mother :
  newborn.table[, 3] <- repro.ewe.table$id
  # extraction id father :
  ram.id.table <- pop.table[pop.table$sexe=="M",]
  id.ram <- ram.id.table[,1]
  newborn.table[, 2] <- sample(id.ram, size = n.newborn, replace = T)
  # extract only females in newborn 
  newborn.table.ewe <- newborn.table[newborn.table$sexe=="F",]
  n.newborn.ewe <- nrow(newborn.table.ewe)
  # control ewe ages and remove old ones
  newpop.table.ewe <- pop.table[!pop.table$age == 8 & pop.table$sexe=="F",]
  n.oldewe <- nrow(pop.table[pop.table$age == 8 & pop.table$sexe=="F",])
  
  # control ram career and replace old ones
  ############################
  table.ram <- pop.table[pop.table$sexe=="M" & pop.table$age > 0,]   # select only males with age > 0
  n.oldram <- as.integer(nrow(table.ram[table.ram$career > career.ram,]) )        # extract old rams
  
  if (n.oldram > 0){
    
    table.oldram <- table.ram[table.ram$career > career.ram,]
    a <- "NR"   # to be adapted
    c <- sample(letters, size = n.oldram, replace = T)
    d <- sample(c(0:9), size = n.oldram, replace = T)
    newram.id <- paste(a,c,d, sep="")
    a <- "RF" 
    c <- sample(letters, size = n.oldram, replace = T)
    d <- sample(c(0:9), size = n.oldram, replace = T)
    newramf.id <- paste(a,c,d, sep="")
    a <- "RM" 
    c <- sample(letters, size = n.oldram, replace = T)
    d <- sample(c(0:9), size = n.oldram, replace = T)
    newramm.id <- paste(a,c,d, sep="")
    age.new.ram <- sample(c(2:8), size = n.oldram, replace = T)
    table.oldram[,1] <- newram.id
    table.oldram[,2] <- newramf.id
    table.oldram[,3] <- newramm.id
    table.oldram[,4] <- age.new.ram
    table.oldram[,5] <- "M"
    table.oldram[,6] <- 0
    
    table.keptram <- table.ram[!table.ram$career > career.ram,]
    newpop.table.ram <- rbind(table.oldram, table.keptram)
  } 
  else {
    newpop.table.ram <-table.ram
  }
  # bind the two tables after 
  newpop.table <- rbind(newpop.table.ewe, newpop.table.ram)
  # check number of animals in herd and add only enough new ewes
  n.oldpop <- nrow(newpop.table)
  n.newpop.total <- n.oldpop + n.newborn.ewe
  surplus <- n.newpop.total-n.pop
  n.newewe.kept <- n.newborn.ewe-surplus
  final.ewe.table <- newborn.table.ewe[sample(nrow(newborn.table.ewe), n.newewe.kept), ]
  
  # stack all df to create new pop for next step:
  pop.table.1 <- rbind(final.ewe.table, newpop.table)
  
  print(paste("gen",i,"completed", sep = " "))
  print(paste("total number birth for gen",i,":", n.newborn, sep =" "))
  print(paste("total ewe birth for gen",i,":", n.newborn.ewe, sep =" "))
  print(paste("total number of old ewe removed for gen",i,":", n.oldewe, sep =" "))
  print(paste("total number of old rams removed for gen",i,":", n.oldram, sep =" "))
  
  # save data
  #save(pop.table.1, file = paste("gen", i, ".Rdata", sep = ""))
  write.csv(pop.table.1, file = paste("gen", i, ".csv", sep = ""), row.names=F)
  
  # increment age at the end of the cycle
  pop.table.1$age <- pop.table.1$age +1
  # increment career of all animals
  pop.table.1$career <- pop.table.1$career+1
  
  pop.table <- pop.table.1
  print(nrow(pop.table))
}

