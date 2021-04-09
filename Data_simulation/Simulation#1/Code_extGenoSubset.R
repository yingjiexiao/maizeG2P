# Extract a simulation genotype dataset including 100K SNPs with minor genotype frequency around 0.3 in the CUBIC population of 1404 inbred lines.

rm(list=ls())

x=list()
for(i in 1:10){
x1=read.delim(paste('/public/home/yjxiao/cubic/allfreq/chr',i,'_freq_progeny.txt',sep=''),stringsAs=F)
snpid=read.delim(paste('/public/home/yjxiao/cubic/emaize/snpid/chr',i,'_snpid',sep=''),stringsAs=F)
stopifnot(all(snpid$rs%in%x1$snp))

x[[i]]=x1[x1$snp%in%snpid$rs,]
}

xs=do.call(rbind,x)
# write.table(xs,'cubic_freq4.5M.txt',sep='\t',row.names=F,quote=F)
# saveRDS(xs,file='cubic_freq4.5M.rds')
xs=xs[xs$het_count<5,]

set.seed(2020)
xs1=xs[sort(sample(1:nrow(xs),100000,replace=FALSE)),]
nc=sum(xs1$minor_count>400 & xs1$minor_count<450)
write.table(xs1,'cubic_freq100K.txt',sep='\t',row.names=F,quote=F)

cat('========= output done! ======\n')
cat('=== allele frequence ~30%:',nc,'\n')



