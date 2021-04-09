#!/bin/sh

# simulate 20 QTN effect following an exponential distribution of a^n (a=0.95),deliver the effects to 20 random SNP positions. The weighted summation of 20 SNP values is set to be genotypic value for a specific line. This process repeats for 1000 times.


rm(list=ls())

set.seed(20210108)
a=0.96
add=a^(1:20)
h2=c(0.1,0.2,0.4,0.6,0.8)

x=read.delim('cubicQTNsgeno.txt',stringsAs=F,row.names=1,check.names=F)

qtn=list()
y=list()
r2=list()
y0.1=list()
y0.2=list()
y0.4=list()
y0.6=list()
y0.8=list()

for(i in 1:1000){
if(i%%10==0) cat('=====',i,' of 1000 simulations ====\n')
id=sample(1:nrow(x),20,replace=FALSE)
qtn[[i]]=rownames(x)[id]
qtn_mat=t(x[id,])

blup=qtn_mat%*%as.matrix(add)
add_mat=matrix(rep(add,nrow(qtn_mat)),nrow(qtn_mat),byrow=T)
qtn_va=apply(qtn_mat*add_mat,2,var)

va=var(blup[,1])
y[[i]]=blup
r2[[i]]=qtn_va/va

# produce the simulated phenotype by adding residuals to blup, with different h2.
y0.1[[i]]=blup+rnorm(nrow(qtn_mat),0,sqrt(va*(1-0.1)/0.1))
y0.2[[i]]=blup+rnorm(nrow(qtn_mat),0,sqrt(va*(1-0.2)/0.2))
y0.4[[i]]=blup+rnorm(nrow(qtn_mat),0,sqrt(va*(1-0.4)/0.4))
y0.6[[i]]=blup+rnorm(nrow(qtn_mat),0,sqrt(va*(1-0.6)/0.6))
y0.8[[i]]=blup+rnorm(nrow(qtn_mat),0,sqrt(va*(1-0.8)/0.8))
}

qtns=data.frame(id=1:20,eff=add,do.call(cbind,qtn))
names(qtns)[-1:-2]=paste('s',1:1000,sep='')
r2s=data.frame(id=1:20,eff=add,do.call(cbind,r2))
names(r2s)[-1:-2]=paste('s',1:1000,sep='')

write.table(qtns,'Sim20QTN_ID.txt',sep='\t',row.names=F,quote=F)
write.table(r2s,'Sim20QTN_R2.txt',sep='\t',row.names=F,quote=F)

y1s=data.frame(line=colnames(x),do.call(cbind,y))
y0.1s=data.frame(line=colnames(x),do.call(cbind,y0.1))
y0.2s=data.frame(line=colnames(x),do.call(cbind,y0.2))
y0.4s=data.frame(line=colnames(x),do.call(cbind,y0.4))
y0.6s=data.frame(line=colnames(x),do.call(cbind,y0.6))
y0.8s=data.frame(line=colnames(x),do.call(cbind,y0.8))
names(y1s)[-1]=paste('s',1:1000,sep='')
names(y0.1s)[-1]=paste('s',1:1000,sep='')
names(y0.2s)[-1]=paste('s',1:1000,sep='')
names(y0.4s)[-1]=paste('s',1:1000,sep='')
names(y0.6s)[-1]=paste('s',1:1000,sep='')
names(y0.8s)[-1]=paste('s',1:1000,sep='')

# randomly select a subset samples with different sizes.
id207=sample(1:nrow(y1s),207,replace=FALSE)
id400=sample(1:nrow(y1s),400,replace=FALSE)
id600=sample(1:nrow(y1s),600,replace=FALSE)
id800=sample(1:nrow(y1s),800,replace=FALSE)
id1000=sample(1:nrow(y1s),1000,replace=FALSE)

#output simulated phenos with series setups.
write.table(y1s,'./1428/simCUBIC1428_h1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.1s,'./1428/simCUBIC1428_h0.1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.2s,'./1428/simCUBIC1428_h0.2_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.4s,'./1428/simCUBIC1428_h0.4_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.6s,'./1428/simCUBIC1428_h0.6_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.8s,'./1428/simCUBIC1428_h0.8_1000phe.txt',sep='\t',row.names=F,quote=F)

write.table(y1s[id207,],'./207/simCUBIC207_h1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.1s[id207,],'./207/simCUBIC207_h0.1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.2s[id207,],'./207/simCUBIC207_h0.2_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.4s[id207,],'./207/simCUBIC207_h0.4_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.6s[id207,],'./207/simCUBIC207_h0.6_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.8s[id207,],'./207/simCUBIC207_h0.8_1000phe.txt',sep='\t',row.names=F,quote=F)

write.table(y1s[id400,],'./400/simCUBIC400_h1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.1s[id400,],'./400/simCUBIC400_h0.1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.2s[id400,],'./400/simCUBIC400_h0.2_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.4s[id400,],'./400/simCUBIC400_h0.4_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.6s[id400,],'./400/simCUBIC400_h0.6_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.8s[id400,],'./400/simCUBIC400_h0.8_1000phe.txt',sep='\t',row.names=F,quote=F)

write.table(y1s[id600,],'./600/simCUBIC600_h1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.1s[id600,],'./600/simCUBIC600_h0.1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.2s[id600,],'./600/simCUBIC600_h0.2_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.4s[id600,],'./600/simCUBIC600_h0.4_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.6s[id600,],'./600/simCUBIC600_h0.6_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.8s[id600,],'./600/simCUBIC600_h0.8_1000phe.txt',sep='\t',row.names=F,quote=F)

write.table(y1s[id800,],'./800/simCUBIC800_h1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.1s[id800,],'./800/simCUBIC800_h0.1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.2s[id800,],'./800/simCUBIC800_h0.2_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.4s[id800,],'./800/simCUBIC800_h0.4_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.6s[id800,],'./800/simCUBIC800_h0.6_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.8s[id800,],'./800/simCUBIC800_h0.8_1000phe.txt',sep='\t',row.names=F,quote=F)

write.table(y1s[id1000,],'./1000/simCUBIC1000_h1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.1s[id1000,],'./1000/simCUBIC1000_h0.1_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.2s[id1000,],'./1000/simCUBIC1000_h0.2_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.4s[id1000,],'./1000/simCUBIC1000_h0.4_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.6s[id1000,],'./1000/simCUBIC1000_h0.6_1000phe.txt',sep='\t',row.names=F,quote=F)
write.table(y0.8s[id1000,],'./1000/simCUBIC1000_h0.8_1000phe.txt',sep='\t',row.names=F,quote=F)

cat('===== all phenos simlation done! =====\n')


