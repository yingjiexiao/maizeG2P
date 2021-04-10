# maizeG2P

######################################################################################
#####   Maize G2P project on a large-scale NCII design with 42K hybrids			######
######################################################################################

# Project Introduction:
1.  Heterosis has been a century-old mystery in biology, which has attracted successive generations of scientists to exploit. Many different explanations of heterosis have been raised most likely due to the fact that, heterosis differently expresses while the studied species, populations and traits are different. Therefore, a unified theoretical framework regarding the secret of heterosis is still being expected. In this work, we present multiple novel findings that may refresh the current view of heterosis.
2. We propose a new genetic design of thirty sets of half-sibling F1 populations, containing a total of 42,840 hybrids for us to study heterosis. These F1 populations were created by crossing the 1,428 CUBIC inbred lines (Liu et al., 2020, Genome Biology) as a maternal pool, with 30 paternal inbred lines world-widely collected. We collected flowering time, plant height and ear weight across multiple environments, up to ~2,570,400 phenotypic data records. This large dataset of F1 populations presents by far the most diverse genetic backgrounds to discover heterosis-determinant genes, exploit the genetic mechanism of heterosis and build the theoretical foundation to apply Big Data-driven genomic design breeding in near future.
3. The genotype (~ 4.5M SNPs) and phenotype of 1458 parents and 42820 F1 hybrids were deposited publicly in the ZEAMAP database (http://zeamap.hzau.edu.cn/ftp/99_MaizegoResources/01_CUBIC_related/). 
4. The scripts for statistical analyses such as phenotypic data normalization, BLUP calculation across multiple environments, genomic prediction, data simulation are provided. 

# 1. BLUP calculation
The best linear unbiased prediction (BLUP) is a widely acceptable approach to extract the genetic value from the multi-environment phenotypes. In the present study, we calculated the BLUP value of each parent and hybrid across 5 environments in each year. The BLUP values are treated to be phenotypes for genomic prediction (G2P) and GWAS analyses. Here, we provide one sample data "Ncii_2015_5locs_parents.txt" for testing the R scripts (blup.r and Onewayblueh2.r).

# 2. Phenotypic data normalization:
1. The "normalization.r" is a user-friendly R code which can be adapted to any R version without any packages dependencies.
2. We provide one R script and several sample data for testing the phenotypic normalization process. The "MJ1221.txt", "MZ1221.txt" files contain the BLUP values of phenotypic data collected in the years of 2014, while the "pop6210.txt" file contains the BLUP values of phenotypic data collected in the years of 2015. 
3. All things needed is to run the "normalization.r", in which, the program will import three sample data "MJ1221.txt", "MZ1221.txt", "pop6210.txt" and export the normalized phenotypes in "pop8652.Nor.txt" file.

# 3. Genomic Prediction:
1. The "gs.r" is a user-friendly R code which can be adapted to any R version
2. Our code need to install the ‘sommer’ package in R, in which, the 'mmer' function was used to fitted the gBLUP model. The version of ‘sommer’ package is 3.0.
3. We provide one R script and several sample data for testing. The "pop8652.Nor.txt" file contain the normalized phenotypic data for 8652 hybrids. The "geno.id" and "test.id" contain a vector of hybrid IDs for training and testing data for genomic prediction. The "add.cXX.zip" file contain a 8652x8652 centered kinship matrix for pairwise hybrid relationship derived from ~4.5M SNPs. Please uncompress the zip file before testing the R script, will obtain a kinship file "add.cXX.txt".
4. All things needed is to run the "gs.r", in which, the program will import four sample data "geno.id", "test.id", "pop8652.Nor.txt", "add.cXX.txt" and export the predicted phenotypes of testing hybrids. If the testing hybrids had measured phenotypes, the program will also output the Pearson correlation coefficient (r) between predicted phenotype and observed phenotype as prediction accuracy in "cors.txt" file.

# 4. Data simulation:
1) We design two simulation schemes for testing GWAS power based on G2P-predicted phenotype. The details on simulation procedures can be seen in the document "./Data_simulation/Simulation details.md".
2) simulation#1. We de novo simulated phenotype by additing polygenic effects of 20 QTNs and random residuals. We then simulated predicted phenotypes with six levels of prediction accuracy to evaluate the influence of G2P accuracy on GWAS power and FDR.
3) simulation#2. We simulated phenotype by adding one QTN effect to a realistic phenotype and predicted phenotype. The process was repeated by testing the impact of different effect sizes of QTN on GWAS power of detecting the specific QTNs in the context of real and predicted phenotypic architectures.
