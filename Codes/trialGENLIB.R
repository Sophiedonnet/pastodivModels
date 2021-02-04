my_family <- as.data.frame(matrix(NA,1,4))
names(my_family)<- c('ind','father','mother','sex')
my_family[1,] <- c(1,0,0,1)
my_family[2,] <- c(2,0,0,2)
my_family[3,] <- c(3,1,2,2)
my_family[4,] <- c(4,1,3,2)
my_family[5,] <- c(5,1,4,1)

#### transform my data frame into a genealogy object
gen_my_family <- gen.genealogy(my_family)

#### plot
gen.graph(gen_my_family)

#### compute inbreeding of individual k

inb<- gen.f(gen_my_family,c(1:5))
inb
