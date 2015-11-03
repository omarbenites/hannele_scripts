##Data Sample for alpha design

library(agricolae)
library(nlme)
 
Genotype<-paste("geno",1:30,sep="")
ntr<-length(Genotype)
r<-2
k<-3
s<-10
obs<-ntr*r
b <- s*r
book<-design.alpha(Genotype,k,r,seed=5)
book$book[,3]<- gl(20,3)
# dataset
y<-c(5,2,7,6,4,9,7,6,7,9,6,2,1,1,3,2,4,6,7,9,8,7,6,4,3,2,2,1,1,2,
     1,1,2,4,5,6,7,8,6,5,4,3,1,1,2,5,4,2,7,6,6,5,6,4,5,7,6,5,5,4)
dbook<-data.frame(book$book,yield=y)

