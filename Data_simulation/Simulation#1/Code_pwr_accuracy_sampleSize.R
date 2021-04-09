#!/bin/sh

# this R script is used to draw plot for QTN detection power against phenotype prediction accuracy and sample size.


rm(list=ls())
require(ggplot2)
require(plyr)

inf=read.delim('./Sim20QTN_R2.txt',stringsAs=F)
inf1=inf[,3:(3+499)]
infr2=data.frame(id=1:20,r2=rowSums(inf1)/500)
infr2$grp=ifelse(infr2$r2>=0.07,1,ifelse(infr2$r2>=0.03,2,3))
x=readRDS('./power20QTN_h2_size.rds')
x1=x[,-1:-4]
x2=ifelse(x1<1e-5,1,0)

dt=data.frame(x[,1:4],pwr=rowSums(x2)/500)
dt$grp=infr2$grp[match(dt$id,infr2$id)]

# dt1=ddply(dt,.(h2,size,grp),summarise,pwr=mean(pwr))
dt1=ddply(dt,.(h2,size),summarise,pwr=mean(pwr))

ggplot(dt1,aes(size,pwr))+
geom_line(aes(col=factor(h2),linetype=factor(h2)),size=0.8)+
# facet_wrap(~grp)+
theme_bw()+
scale_x_continuous(breaks=c(207,400,600,800,1000,1428))+
scale_y_continuous(limits=c(0,1))+
scale_colour_manual(values=c("#8491B4FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#E64B35FF"),labels=c('r2=0.1','r2=0.2','r2=0.4','r2=0.6','r2=0.8','Simulated Phenotype'),guide=guide_legend(reverse=T))+
scale_linetype_manual(values=c(rep('longdash',5),'solid'),labels=c('r2=0.1','r2=0.2','r2=0.4','r2=0.6','r2=0.8','Simulated Phenotype'),guide=guide_legend(reverse=T))+
labs(x='Population size',y='Power',col=NULL,linetype=NULL)

ggsave('./power20QTN2.pdf')