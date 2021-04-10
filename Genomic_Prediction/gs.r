
# The script was used to predict the phenotypes of testing samples and derive predictability by the genomic best linear unbiased prediction (gBLUP) model. 
# In this paper, gBLUP model was to infer the phenotypes of the 34,188 F1 combinations we did not phenotype in the field.
#   8652 F1 hybrids with measured phenotypes were used to train the G2P model, while remaining 1221 samples of each paternal line were the predicted population.
#   So we had done this procedure 28 times to obtain the predicted phenotypes of the 34,188 F1 combinations.
# Notes:
# geno_id: samples IDs in the genotypic data, including training and testing lines
# test_id: samples IDs in the testing dataset
# pop8652: normalized phenotypic data for 8652 samples which had observed phenotypes
# A.kin: the centered genomic relationship matrix (cGRM) deduced by the GEMMA (Genome-wide Efficient Mixed Model Association) software
# traits: those need to be predicted agronomic or derived traits in GS model

# load the dependency R packages
library("sommer")

# input data
geno_id <- read.table("./data/geno.id", header=F)
test_id <- read.table("./datatest.id", header=F)
pop8652 <- read.table("./data/pop8652.Nor.txt", header=T)
## kinship
A.kin <- read.table("./data/add.cXX.txt", header=F)
rownames(A.kin) <- geno_id$V1
colnames(A.kin) <- geno_id$V1
A.kin <- as.matrix(A.kin)
## predicted traits
traits <- c("DTT","PH","EW","MPH.DTT","MPH.PH","MPH.EW")

# prepare phenotypes required in GS model, including training and testing samples, in which, phenotypes of testing lines were NA
if(length(setdiff(geno_id$V1,pop8652$line))==0){
  pheno <- pop8652[match(geno_id$V1, pop8652$line),]
}else{
  test <- data.frame(line=test_id$V1)
  test$DTT <- NA
  test$PH <- NA
  test$EW <- NA
  test <- data.frame(test, pop8652[match(gsub("_X_.*", "", test$line), gsub("_X_.*", "", pop8652$line)), ][,c('DTT_F', 'PH_F', 'EW_F')])
  test <- data.frame(test, pop8652[match(gsub(".*_X_", "", test$line), gsub(".*_X_", "", pop8652$line)), ][,c('DTT_M', 'PH_M', 'EW_M')])
  test$MPH.DTT <- NA
  test$MPH.PH <- NA
  test$MPH.EW <- NA
  train_id <- setdiff(geno_id$V1, test_id$V1)
  train <- pop8652[pop8652$line %in% train_id,]
  pheno <- rbind(train, test)
}

# prepare data for gblup model
## an incidence matrix for fixed effects.
X11 <- model.matrix(~DTT_F+DTT_M, data=pheno)
X12 <- model.matrix(~PH_F+PH_M, data=pheno)
X13 <- model.matrix(~EW_F+EW_M, data=pheno)
X1 <- list(X11, X12, X13, X11, X12, X13)
## incidence matrices and var-cov matrices for random effects.
Za <- diag(dim(A.kin)[1])
ETA <- list(add=list(Z=Za,K=A.kin))
vv <- which(pheno$line %in% test_id$V1)

# fit the model
preds <- data.frame(line=pheno$line[vv])
cors <- NULL
for(j in 1:length(traits)){ 
  y <- pheno[,colnames(pheno)==traits[j]]
  y.trn <- y
  y.trn[vv] <- NA
  ## gblup model 
  ans <- mmer(Y=y.trn, X=X1[[j]], Z=ETA)
  ## predicted values
  predj <- data.frame(fitted=as.vector(ans$fitted.y)[vv])
  colnames(predj) <- traits[j]
  preds <- data.frame(preds, predj)
  if(length(setdiff(geno_id$V1,pop8652$line))==0){
    corj <- cor(as.vector(ans$fitted.y)[vv], y[vv], use="complete.obs")
  } else {
    corj <- NA
  }
  cors <- c(cors, corj)
}

# output data
write.table(preds, file="./res/predicts.txt", row.names=F, quote=F, sep="\t")
if(length(setdiff(geno_id$V1,pop8652$line))==0){
  write.table(data.frame(trait=traits, cor=cors), file="./res/cors.txt", row.names=F, quote=F, sep="\t")
}
