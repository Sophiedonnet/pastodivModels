LHall = Reduce(rbind,LHerds)
LHall$genbirth = sapply(LHall$ind,function(identif){
  temp=substr(identif,1,3)
  if (substr(temp,2,2)=="-") temp=substr(temp,1,1)
  if (substr(temp,3,3)=="-") temp=substr(temp,1,2)
  as.numeric(temp)
  })
summary(as.numeric(LHall$genbirth))
LHall$valrepro = compute.val.repro(LHall,selectionWeights=NULL)


library(ggplot2)
LHallsans0 = LHall |> filter(genbirth>0)
ggplot(LHallsans0,aes(x=genbirth,fill=factor(codingGene.5)))+geom_bar(aes(y = (..count..)/sum(..count..)))

ggplot(LHallsans0,aes(x=genbirth,y=valrepro))+geom_point()+geom_smooth()

