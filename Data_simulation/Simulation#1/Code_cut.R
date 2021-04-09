#!/bin/sh

# This script is used to generate the GWAS p-value threshold for different type-I error based on permutation tests.

sz=c(207,400,600,800,1000,1428)
acc=c(0.1,0.2,0.4,0.6,0.8,1)
alpha=seq(0.05,0.95,0.05)

cutGWAS2=list()
for(i in 1:length(sz)){
  cat('====',sz[i],'====\n')  
  
  cutGWAS1=list()
  for(j in 1:length(acc)){
    cat('======',acc[j],'======\n')

    x1=vector()
    x2=vector()
    for(k in 1:500){
    x1[k]=min(read.delim(paste('./',sz[i],'/output/',acc[j],'/',sz[i],'_',acc[j],'_',k,'.assoc.txt',sep=''),stringsAs=F)$p_wald,na.rm=T)
    x2[k]=quantile(read.delim(paste('./',sz[i],'/output/',acc[j],'/',sz[i],'_',acc[j],'_',k,'.assoc.txt',sep=''),stringsAs=F)$p_wald,0.01,na.rm=T,type=1)
    }
    
    x11=-log10(x1)
    x22=-log10(x2)
#    alpha=seq(0.05,0.95,0.05)
    cuts1=quantile(x11,1-alpha,na.rm=T,type=1)
    cuts2=quantile(x22,1-alpha,na.rm=T,type=1)
    cutGWAS1[[j]]=data.frame(acc=acc[j],typeI=alpha,cuts1=cuts1,cuts0.99=cuts2)
  }

  cutGWAS2[[i]]=data.frame(size=sz[i],do.call(rbind,cutGWAS1))
}

cutGWASs=do.call(rbind,cutGWAS2)

write.table(cutGWASs,'./cutGWAS.perm.txt',sep='\t',row.names=F,quote=F)

cat('++++++++++++++++++  output done!  ++++++++++++++++++++\n')




