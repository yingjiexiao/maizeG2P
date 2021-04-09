#!/bin/sh

# This script is to permute the correlations between genotype and phenotype, resulting in a bunch of fake phenotypes for testing the GWAS threshold under different type-I error.

rm(list=ls())
set.seed(20210125)

sz=c(207,400,600,800,1000,1428)
acc=c(0.1,0.2,0.4,0.6,0.8,1)

for(i in 1:length(sz)){
  cat(' ====',sz[i],'====\n')
  id=sample(1:sz[i],sz[i],replace=FALSE)  

  for(j in 1:length(acc)){
  cat(' ======',acc[j],'======\n')
  x=read.delim(paste('./',sz[i],'/simCUBIC',sz[i],'_h',acc[j],'.gemma',sep=''),head=FALSE)
  
   xp=list()
   for(k in 1:1000){
   if(k%%200==0) cat('=======',k,'of 1000 ===== \n')
   id=sample(1:sz[i],sz[i],replace=FALSE)
   xp[[k]]=x[id,k] 
   }
   
  xps=do.call(cbind,xp)
  stopifnot(nrow(xps)==sz[i])
  stopifnot(ncol(xps)==1000) 
  write.table(xps,paste('./',sz[i],'/permCUBIC',sz[i],'_h',acc[j],'.gemma',sep=''),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
}

cat('++++++++++ permutated phenotype make done! +++++++++++++++\n')

  



