
#!/bin/sh

# A script for drawing fdr plot based on a new definition: fdr=N1/(N1+N0), where N0 is predefined number of QTNs and the N1 is the number of signicant loci not included in QTNs as false positives.


# fdr
rm(list=ls())
require(ggplot2)
require(plyr)

x=readRDS('./fdn20QTN_h2_size.rds')
dt=data.frame(x[,1:2],fdn=rowSums(x[,-1:-2])/500)
dt$grp=ifelse(dt$h2==1,1,2)

ggplot(dt,aes(size,fdn/20))+
# geom_point(aes(col=factor(h2)))+
geom_line(aes(col=factor(h2),linetype=factor(h2)),size=0.8)+
scale_x_continuous(breaks=c(207,400,600,800,1000,1428))+
scale_y_continuous(limits=c(0,0.12),breaks=seq(0,0.12,0.02))+
scale_colour_manual(values=c("#8491B4FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#E64B35FF"),labels=c('r2=0.1','r2=0.2','r2=0.4','r2=0.6','r2=0.8','Simulated Phenotype'),guide=guide_legend(reverse=T))+
scale_linetype_manual(values=c(rep('longdash',5),'solid'),labels=c('r2=0.1','r2=0.2','r2=0.4','r2=0.6','r2=0.8','Simulated Phenotype'),guide=guide_legend(reverse=T))+
theme_bw()+
# theme(legend.direction='horizontal',
	# legend.box='')+
labs(x='Population size',y='False positive rate',col=NULL,linetype=NULL)

ggsave('./fdr20QTN2.pdf')










