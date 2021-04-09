
# simu.R

# The script to assess the power of GS predicted phenotype on GWAS.
# 1-data process


rm(list=ls())


male='MG_1522' # Zheng58
# male='MG_1542' # Jing724

x=read.delim(paste0('./data/',male,'_pred1221.txt'))
x1=x[x$index=='Obs',]
x2=x[x$index=='Pred',]

x1=x1[,match(c('line','DTT','PH','EW','DTT_HV','PH_HV','EW_HV'),names(x1))]
x2=x2[,match(c('line','DTT','PH','EW','DTT_HV','PH_HV','EW_HV'),names(x2))]

y=read.delim('./data/pop6210.MZJ2442.Nor_gemma.txt')
y=y[,match(c('seqname','DTT','PH','EW','DTT_HV','PH_HV','EW_HV'),names(y))]
names(y)[1]='line'

# stopifnot(all(x1$line%in%y$line))
stopifnot(all(x$line%in%y$line))
y1=y[match(x1$line,y$line),]

stopifnot(all(x1$DTT==y1$DTT,na.rm=T))
stopifnot(all(x1$PH==y1$PH,na.rm=T))
stopifnot(all(x1$EW==y1$EW,na.rm=T))

z1=x1
z2=rbind(x1,x2)
z3=y[match(z2$line,y$line),]


id1428=read.table('./data/cubic1428_id',stringsAs=F)$V1
id1428s=paste(id1428,male,sep='_X_')

id207=read.delim('./data/cubic207_id_rs',stringsAs=F)[,1]
id207s=paste(id207,male,sep='_X_')

stopifnot(all(id207s%in%z1$line))
stopifnot(all(id1428s%in%z3$line))

# sort phenotype by genotype ids

z11=z1[match(id207s,z1$line),]
z21=z2[match(id1428s,z2$line),]
z31=z3[match(id1428s,z3$line),]

z11$line=id207
z21$line=z31$line=id1428

write.table(z11,'./data/cubic207_obs.txt',sep='\t',row.names=F,quote=F)
write.table(z21,'./data/cubic1428_pred.txt',sep='\t',row.names=F,quote=F)
write.table(z31,'./data/cubic1428_obs.txt',sep='\t',row.names=F,quote=F)


# 2-simulate phenotype by adding QTL effect per trait architecture

rm(list=ls())

x=read.delim('qtlgeno207.txt',stringsAs=F)
y=read.delim('qtlgeno1428.txt',stringsAs=F)

z1=read.delim('../data/cubic207_obs.txt',stringsAs=F)
z2=read.delim('../data/cubic1428_obs.txt',stringsAs=F)
z3=read.delim('../data/cubic1428_pred.txt',stringsAs=F)

trs=names(z1)[-1]

# ds=c(0.2,0.4,0.8,1,1.2,1.4,1.6,1.8,2,2.2)
ds=seq(0.1,2,0.1)
dstext=format(seq(0.1,2,0.1),nsmall=1)
# rn=50
# dn=5
# tn=3

rn=100
dn=length(ds)
tn=3


z10=z1[,tn+1]
z20=z2[,tn+1]
z30=z3[,tn+1]


for(i in 1:dn){

	d=ds[i]
	sz1=list()
	sz2=list()
	sz3=list()
	for(j in 1:rn){
		sz1[[j]]=x[,j]*d+z10
		sz2[[j]]=y[,j]*d+z20
		sz3[[j]]=y[,j]*d+z30
	}

	sz11=data.frame(z10,do.call(cbind,sz1))
	sz21=data.frame(z20,do.call(cbind,sz2))
	sz31=data.frame(z30,do.call(cbind,sz3))

	write.table(sz11,paste('../pheno/new1/',trs[tn],'/sim207_',trs[tn],'_',dstext[i],'.txt',sep=''),sep='\t',row.names=F,col.names=F,quote=F)
	write.table(sz21,paste('../pheno/new1/',trs[tn],'/sim1428obs_',trs[tn],'_',dstext[i],'.txt',sep=''),sep='\t',row.names=F,col.names=F,quote=F)
	write.table(sz31,paste('../pheno/new1/',trs[tn],'/sim1428pred_',trs[tn],'_',dstext[i],'.txt',sep=''),sep='\t',row.names=F,col.names=F,quote=F)
	
}







