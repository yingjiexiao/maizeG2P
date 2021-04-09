

# draw plots for power and fdr


# power

rm(list=ls())
require(ggplot2)
require(plyr)
require(grid)
require(scales)

# plot#1
inf=read.delim('./power/Sim20QTN_R2.txt',stringsAs=F)
inf1=inf[,3:(3+499)]
infr2=data.frame(id=1:20,r2=rowSums(inf1)/500)
x=readRDS('./power/power20QTN_h2_size.rds')
x1=x[,-1:-4]
x2=ifelse(x1<1e-5,1,0)

dt=data.frame(x[,1:4],pwr=rowSums(x2)/500)
dt1=ddply(dt,.(h2,size),summarise,pwr=mean(pwr))

p1=ggplot(dt1,aes(size,pwr))+
geom_line(aes(col=factor(h2),linetype=factor(h2)),size=0.8)+
# facet_wrap(~grp)+
theme_bw()+
theme(
	plot.margin=unit(c(1,4,1,2),'lines'),
	legend.position="none",
)+
scale_x_continuous(breaks=c(207,400,600,800,1000,1428))+
scale_y_continuous(limits=c(0,1))+
scale_colour_manual(values=c("#8491B4FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#E64B35FF"),labels=c('r2=0.1','r2=0.2','r2=0.4','r2=0.6','r2=0.8','Simulated Phenotype'),guide=guide_legend(reverse=T))+
scale_linetype_manual(values=c(rep('longdash',5),'solid'),labels=c('r2=0.1','r2=0.2','r2=0.4','r2=0.6','r2=0.8','Simulated Phenotype'),guide=guide_legend(reverse=T))+
labs(x='Population size',y='Power',col=NULL,linetype=NULL)

# plot #2
xx=readRDS('./fdr/fdn20QTN_h2_size.rds')
dtt=data.frame(xx[,1:2],fdn=rowSums(xx[,-1:-2])/500)

p2=ggplot(dtt,aes(size,fdn/(fdn+20)))+
# geom_point(aes(col=factor(h2)))+
geom_line(aes(col=factor(h2),linetype=factor(h2)),size=0.8)+
scale_x_continuous(breaks=c(207,400,600,800,1000,1428))+
scale_y_continuous(limits=c(0,0.12),breaks=seq(0,0.12,0.02))+
scale_colour_manual(values=c("#8491B4FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#E64B35FF"),labels=c('r2=0.1','r2=0.2','r2=0.4','r2=0.6','r2=0.8','Simulated Phenotype'),guide=guide_legend(reverse=T))+
scale_linetype_manual(values=c(rep('longdash',5),'solid'),labels=c('r2=0.1','r2=0.2','r2=0.4','r2=0.6','r2=0.8','Simulated Phenotype'),guide=guide_legend(reverse=T))+
theme_bw()+
theme(
	# legend.direction='horizontal',
	# legend.box='',
	plot.margin=unit(c(1,1,1,2),'lines'),
)+
labs(x='Population size',y='False positive rate',col=NULL,linetype=NULL)


# plot #3
inf=read.delim('./power/Sim20QTN_R2.txt',stringsAs=F)
inf1=inf[,3:(3+499)]
infr2=data.frame(id=1:20,r2=rowSums(inf1)/500)
infr2$grp=ifelse(infr2$r2>=0.07,3,ifelse(infr2$r2>=0.03,2,1))
x=readRDS('./power/power20QTN_h2_size.rds')
# x1=x[,-1:-4]

cuts=read.delim('./perm/cutGWAS.perm.txt',stringsAs=F)
cuts$para=paste(cuts$acc,cuts$size,sep='_')
x$para=paste(x$h2,x$size,sep='_')
para.uni=unique(x$para)

typeone=seq(0.05,0.95,0.05)

dt1=list()
for(i in 1:length(para.uni)){
    cat('=======',para.uni[i],'======\n')
    cut1=cuts$cuts1[cuts$para==para.uni[i]]
    x1=x[x$para==para.uni[i],-ncol(x)]

    dt=list()
    for(j in 1:length(cut1)){
        cat('==========',cut1[j],'=======\n')
        x2=ifelse(x1[,-1:-4]<10^(-cut1[j]),1,0)
        dt[[j]]=data.frame(typeone=typeone[j],x1[,1:4],pwr=rowSums(x2)/500)
    }

    dt1[[i]]=data.frame(para=para.uni[i],do.call(rbind,dt))
}

dts=do.call(rbind,dt1)
dts$grp=infr2$grp[match(dts$id,infr2$id)]

dts1=ddply(dts,.(typeone,h2,size,grp),summarise,pwr=mean(pwr))
dts1$grp1=ifelse(dts1$grp==1,'PVE<3%',ifelse(dts1$grp==2,'3%<PVE<7%','PVE>7%'))
dts1$grp1=factor(dts1$grp1,levels=c('PVE<3%','3%<PVE<7%','PVE>7%'))

p3=ggplot(dts1[dts1$size%in%c(207,600,1428),],aes(typeone,pwr))+
# ggplot(dts1[dts1$size%in%c(400,800,1000),],aes(typeone,pwr))+
geom_line(aes(col=factor(h2),linetype=factor(h2)),size=0.8)+
facet_grid(size~grp1)+
theme_bw()+
theme(
    legend.position="none",
    # legend.justification=c(-0.1,0.8),
    strip.text=element_text(size=12),
    legend.text=element_text(size=10),
    axis.text=element_text(size=10),
    axis.title=element_text(size=12),
    # legend.spacing=unit(0,'cm'),
    # legend.margin=margin(0,0,0,0,'cm'),
	plot.margin=unit(c(1,1,1,1),'lines'),
)+
scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2),label=format(seq(0,1,0.2),nsmall=1))+
scale_colour_manual(values=c("#8491B4FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#E64B35FF"),labels=c('r2=0.1','r2=0.2','r2=0.4','r2=0.6','r2=0.8','Simulated Phenotype'),guide=guide_legend(reverse=T))+
scale_linetype_manual(values=c(rep('longdash',5),'solid'),labels=c('r2=0.1','r2=0.2','r2=0.4','r2=0.6','r2=0.8','Simulated Phenotype'),guide=guide_legend(reverse=T))+
labs(y='Power',x='Type I error',col=NULL,linetype=NULL)



pdf(paste('./joint_pwr_fdr.pdf'),height=8,width=8)
grid.newpage()
pushViewport(viewport(layout=grid.layout(5,4)))
vplayout = function(x,y){
  viewport(layout.pos.row=x,layout.pos.col=y)
}

print(p1,vp=vplayout(1:2,1:2))
print(p2,vp=vplayout(1:2,3:4))
print(p3,vp=vplayout(3:5,1:4))

dev.off()

