## readme.txt

1. System requirements

The "gs.r" is a user-friendly R code which can be adapted to any R version

2. Intallation guide
 
Our code need to install the ¡®sommer¡¯ package in R, in which, the 'mmer' function was used to fitted the gBLUP model. The version of ¡®sommer¡¯ package is 3.0.

3 and 4. Demo and Intructions for use

We provide one R scripts: gs.r.
All things need to do is run the "gs.r", in which, the program will import four sample data "geno.id", "test.id", "pop8652.Nor.txt", "add.cXX.txt" and export the predicted phenotypes of testing samples in "pop8652.Nor.txt" file.
If the testing lines had measured phenotypes, the program will also import the pearson correlation coefficient (r) between predicted phenotype and observed phenotype as prediction accuracy in "cors.txt" file.
