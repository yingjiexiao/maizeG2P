
# pwr.R

rm(list=ls())

# for 207 pop

x01=read.delim(paste('./207/DTT/output/0.1/DTT_0.1_1.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
x02=read.delim(paste('./207/PH/output/0.1/PH_0.1_1.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
x03=read.delim(paste('./207/EW/output/0.1/EW_0.1_1.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]

dstext=format(seq(0.1,2,0.1),nsmall=1)

pd1=list()
pd2=list()
pd3=list()
for(j in 1:length(dstext)){
	d=dstext[j]
	cat(d,'\n')

	p1=vector()
	p2=vector()
	p3=vector()
	for(i in 1:100){

	x11=read.delim(paste('./207/DTT/output/',d,'/DTT_',d,'_',i+1,'.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
	x12=read.delim(paste('./207/PH/output/',d,'/PH_',d,'_',i+1,'.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
	x13=read.delim(paste('./207/EW/output/',d,'/EW_',d,'_',i+1,'.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]

	p1[i]=x11$p_wald[i]
	p2[i]=x12$p_wald[i]
	p3[i]=x13$p_wald[i]	

	}

	pd1[[j]]=p1
	pd2[[j]]=p2
	pd3[[j]]=p3

}

xd1=data.frame(x01,do.call(cbind,pd1))
xd2=data.frame(x02,do.call(cbind,pd2))
xd3=data.frame(x03,do.call(cbind,pd3))


pwr1=apply(xd1[,-1:-2],2,function(x)sum(x<1/3e6)/100)
pwr2=apply(xd2[,-1:-2],2,function(x)sum(x<1/3e6)/100)
pwr3=apply(xd3[,-1:-2],2,function(x)sum(x<1/3e6)/100)

dt1=data.frame(pop=207,trait=rep(c('DTT','PH','EW'),each=20),effect=rep(seq(0.1,2,0.1),3),pwr=c(pwr1,pwr2,pwr3))
write.table(dt1,'pwr_207.txt',sep='\t',row.names=F,quote=F)

# ggplot(dt1,aes(effect,pwr))+geom_line(aes(col=trait))+geom_point(aes(shape=trait,col=trait))


############################
# for 1428 prediction pop

rm(list=ls())

x01=read.delim(paste('./1428pred/DTT/output/0.1/DTT_0.1_1.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
x02=read.delim(paste('./1428pred/PH/output/0.1/PH_0.1_1.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
x03=read.delim(paste('./1428pred/EW/output/0.1/EW_0.1_1.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]

dstext=format(seq(0.1,2,0.1),nsmall=1)

pd1=list()
pd2=list()
pd3=list()
for(j in 1:length(dstext)){
	d=dstext[j]
	cat(d,'\n')

	p1=vector()
	p2=vector()
	p3=vector()
	for(i in 1:100){

	x11=read.delim(paste('./1428pred/DTT/output/',d,'/DTT_',d,'_',i+1,'.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
	x12=read.delim(paste('./1428pred/PH/output/',d,'/PH_',d,'_',i+1,'.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
	x13=read.delim(paste('./1428pred/EW/output/',d,'/EW_',d,'_',i+1,'.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]

	p1[i]=x11$p_wald[i]
	p2[i]=x12$p_wald[i]
	p3[i]=x13$p_wald[i]	

	}

	pd1[[j]]=p1
	pd2[[j]]=p2
	pd3[[j]]=p3

}

xd1=data.frame(x01,do.call(cbind,pd1))
xd2=data.frame(x02,do.call(cbind,pd2))
xd3=data.frame(x03,do.call(cbind,pd3))

pwr1=apply(xd1[,-1:-2],2,function(x)sum(x<1/3e6)/100)
pwr2=apply(xd2[,-1:-2],2,function(x)sum(x<1/3e6)/100)
pwr3=apply(xd3[,-1:-2],2,function(x)sum(x<1/3e6)/100)

dt1=data.frame(pop='1428pred',trait=rep(c('DTT','PH','EW'),each=20),effect=rep(seq(0.1,2,0.1),3),pwr=c(pwr1,pwr2,pwr3))
write.table(dt1,'pwr_1428pred.txt',sep='\t',row.names=F,quote=F)

############
# for 1428 obs pop

rm(list=ls())

x01=read.delim(paste('./1428obs/DTT/output/0.1/DTT_0.1_1.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
x02=read.delim(paste('./1428obs/PH/output/0.1/PH_0.1_1.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
x03=read.delim(paste('./1428obs/EW/output/0.1/EW_0.1_1.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]

dstext=format(seq(0.1,2,0.1),nsmall=1)

pd1=list()
pd2=list()
pd3=list()
for(j in 1:length(dstext)){
	d=dstext[j]
	cat(d,'\n')

	p1=vector()
	p2=vector()
	p3=vector()
	for(i in 1:100){

	x11=read.delim(paste('./1428obs/DTT/output/',d,'/DTT_',d,'_',i+1,'.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
	x12=read.delim(paste('./1428obs/PH/output/',d,'/PH_',d,'_',i+1,'.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]
	x13=read.delim(paste('./1428obs/EW/output/',d,'/EW_',d,'_',i+1,'.assoc.txt',sep=''),stringsAs=F)[,c(2,12)]

	p1[i]=x11$p_wald[i]
	p2[i]=x12$p_wald[i]
	p3[i]=x13$p_wald[i]	

	}

	pd1[[j]]=p1
	pd2[[j]]=p2
	pd3[[j]]=p3

}

xd1=data.frame(x01,do.call(cbind,pd1))
xd2=data.frame(x02,do.call(cbind,pd2))
xd3=data.frame(x03,do.call(cbind,pd3))

pwr1=apply(xd1[,-1:-2],2,function(x)sum(x<1/3e6)/100)
pwr2=apply(xd2[,-1:-2],2,function(x)sum(x<1/3e6)/100)
pwr3=apply(xd3[,-1:-2],2,function(x)sum(x<1/3e6)/100)

dt1=data.frame(pop='1428obs',trait=rep(c('DTT','PH','EW'),each=20),effect=rep(seq(0.1,2,0.1),3),pwr=c(pwr1,pwr2,pwr3))
write.table(dt1,'pwr_1428obs.txt',sep='\t',row.names=F,quote=F)

####################

# test power for three schemes

rm(list=ls())

dt1=read.delim('pwr_207.txt')
dt2=read.delim('pwr_1428pred.txt')
# dt3=read.delim('pwr_1428obs.txt')

require(ggplot2)

dt=rbind(dt1,dt2)

ggplot(dt,aes(effect,pwr))+
geom_line(aes(col=trait,linetype=pop))+
geom_point(aes(col=trait))+
# facet_wrap(~pop)+
theme(

	)+
scale_x_continuous(breaks=seq(0.1,2,0.1),label=seq(0.1,2,0.1))+
scale_linetype_discrete(label=c('207 lines+1221 predicted','207 lines'))+
labs(x='Genetic effect (times s.d.)',y='Power',col=NULL,linetype=NULL)

# ggsave('power_GS_GWAS.png')
ggsave('power_GS_GWAS.pdf')
































