
# The script was used to normalize and combine types of data considering the population structure, multi-environment
# Notes:
# MJ1221: the blups values of Jing724 F1 hybrids in 2014
# MZ1221: the blups values of Zheng58 F1 hybrids in 2014
# pop6210: the blups values of 6210 F1 hybrids resulting from the crossing of 207 maternal lines with the 30 paternal testers in 2015
# Females of MJ1221 and MZ1221 are the same, but different with pop6210; the total number of maternal lines is 1428
# MPH: middle-parent heterosis, _F: female, _M: male
# The normalized data or its subsets were subsequently used in G2P predictions and gwas analysis
 

# input data, observed phenotype
MJ1221 <- read.table("./data/MJ1221.txt", header=TRUE)
MZ1221 <- read.table("./data/MZ1221.txt", header=TRUE)
pop6210 <- read.table("./data/pop6210.txt", header=TRUE)

# normalized BLUP values of MJ1221/MZ1221
MJ1221.Nor <- data.frame(line=MJ1221$seqname, apply(MJ1221[ ,-c(1,8:10)], 2, scale))
MZ1221.Nor <- data.frame(line=MZ1221$seqname, apply(MZ1221[ ,-c(1,8:10)], 2, scale))

# normalized BLUP values of parents from pop6210
pop6210$male <- gsub(".*_X_", "", pop6210$seqname)
pop6210_M <- pop6210[,c('male','DTT_M','PH_M','EW_M')][!duplicated(pop6210[,c('male','DTT_M','PH_M','EW_M')]),]
pop6210_M.Nor <- data.frame(male=pop6210_M$male, apply(pop6210_M[ ,-1], 2, scale))

# normalized BLUP values of pop6210
pop6210.Nor <- NULL
for(i in 1:length(unique(pop6210$male))){
  malei <- subset(pop6210, pop6210$male==unique(pop6210$male)[i])
  malei.Nor <- data.frame(line=malei$seqname, apply(malei[ ,-c(1,8:10,14)], 2, scale))
  pop6210.Nor <- rbind(pop6210.Nor, malei.Nor)
}

# combine 
pop8652.Nor <- rbind(pop6210.Nor, MZ1221.Nor, MJ1221.Nor)
pop8652_M.Nor <- pop6210_M.Nor[match(gsub(".*_X_", "", pop8652.Nor$line), pop6210_M.Nor$male), ]
res <- data.frame(pop8652.Nor[,-c(8:10)], pop8652_M.Nor[,-1], pop8652.Nor[,8:10])
write.table(res, "./res/pop8652.Nor.txt", row.names = F, quote=F, sep="\t")      #output data
