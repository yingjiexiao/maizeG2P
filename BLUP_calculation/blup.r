# blup.r

rm(list=ls())

source('Onewaybluph2.r')

x=read.delim('./data/Ncii_2015_5locs_parents.txt')

names(x)[1:2]=c('LOC','LINE')

x1=Onewaybluph2(x)

write.table(x1,'./res/Ncii_2015_blup_parents',sep='\t',row.names=F,quote=F)

