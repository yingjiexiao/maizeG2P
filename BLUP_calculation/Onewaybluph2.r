### blup estimates for phenotypic traits across multiple envionments##
### Broad-sense heritability estimation across multiple envionments
### THE DATA SHOULD INCLUDE THREE COLUMNS AT LEAST,DENOTED AS 'LINE','LOC' AND 'TRAIT1...TRAITn'
### All comments, please contact Yingjie Xiao using Email: shanren0179@163.com

Onewaybluph2<-function(x){

##Calculate variance components
##require package lme4
library(lme4)
id<-as.character(x$LINE[x$LOC%in%unique(x$LOC)[1]])
line.blup<-list()
heritability<-list()

for(i in 3:ncol(x)){
	cat(names(x)[i],'start\n')
	varcomp<-lmer(x[,i]~(1|LINE)+(1|LOC),data=x) #solve total variance into components
	#Extract genotypic variance and residual variance and calculate broad-sense heritability
	var.trans = lme4::VarCorr(varcomp)
	# var<-as.data.frame(summary(varcomp)@REmat,row.names=T)
	var = data.frame(Groups=c('LINE','LOC','Residual'),Variance=c(as.numeric(var.trans$LINE),
			as.numeric(var.trans$LOC),attr(var.trans,'sc')^2),check.names=F)
	gvar<-as.numeric(as.character(var$Variance))[var$Groups%in%'LINE']
	evar<-as.numeric(as.character(var$Variance))[var$Groups%in%'Residual']
	heritability[[i-2]]<-matrix(gvar/(gvar+evar/length(unique(x$LOC))),1,1,dimnames=list(colnames(x)[i],'H2'))

	f<-fixef(varcomp)
	r<-ranef(varcomp)$LINE
	blup<-f+r
	line.blup[[i-2]]<-blup[match(id,rownames(blup)),]
	}
	heritability=do.call('rbind',heritability)
	line.blup=data.frame(id,do.call('cbind',line.blup))
	colnames(line.blup)<-c('line',names(x)[-c(1:2)])

# cat('----------------------------OUTPUT RESULTS OF EACH TRAITS INTO FILES:------------------------------------------------------------------','\n',
# '----------------------------------------------------------------------BroadSenseHeritability.txt-------------------------------------------','\n',
# '----------------------------------------------------------------------BlupEachlines.txt-----------------------------------------------------','\n')
# write.table(heritability,'BroadSenseHeritability.txt',sep='\t',quote=F)
# write.table(line.blup,'BlupEachlines.txt',sep='\t',quote=F,row.names=F)
return(line.blup)
}
