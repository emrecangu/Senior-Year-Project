setwd("e:/bitirme_tezi/maincohort")
#1 Emre SKAT
library(SKAT)
Close_SSD()
File.Bed <- "emreobesity.bed"
File.Bim <- "emreobesity.bim"
File.Fam <- "emreobesity.fam"
File.SetID <- "emreobesitysetID.setID"
File.SSD <- "emreobesity.SSD"
File.Info <- "emreobesity.SSD.info"

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)

FAM <- Read_Plink_FAM(File.Fam, Is.binary=TRUE)

y <- FAM$Phenotype

SSD.INFO<-Open_SSD(File.SSD, File.Info)

obj<-SKAT_Null_Model(y ~ 1, out_type="D", n.Resampling = 1000, type.Resampling = "bootstrap")

out <- SKAT.SSD.All(SSD.INFO, obj)
warnings(10000000)
output.df = out$results
write.table(output.df, file = "./emreobesityresultsbinaryouttypeC")

# 2 - Burden with SKAT-O weights and dichotomous
library(SKAT)
Close_SSD()
File.Bed <- "emreobesity.bed"
File.Bim <- "emreobesity.bim"
File.Fam <- "emreobesity.fam"
File.SetID <- "emreobesitysetID.setID"
File.SSD <- "emreobesity.SSD"
File.Info <- "emreobesity.SSD.Info"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFO <- Open_SSD(File.SSD, File.Info)
FAM<-Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)
y <- FAM$Phenotype
obj <- SKAT_Null_Model(y ~ 1, data=NULL, out_type="D", type.Resampling="bootstrap", n.Resampling=0, Adjustment=TRUE)
out.skat <- SKAT.SSD.All(SSD.INFO, obj, kernel = "linear.weighted", 
                         method="Burden", weights.beta=c(1,25), weights=NULL, 
                         impute.method="random", r.corr=0, is_check_genotype=TRUE,
                         is_dosage = FALSE, missing_cutoff=0.9 , max_maf=1, estimate_MAF=1)
output.df = out$results
write.table(output.df, file = "./emreobesityresultsskatodichotomousemrehoca")

# 3 - SKAT binary
library(SKAT)
Close_SSD()
File.Bed <- "emreobesity.bed"
File.Bim <- "emreobesity.bim"
File.Fam <- "emreobesity.fam"
File.SetID <- "emreobesitysetID.setID"
File.SSD <- "emreobesity.SSD"
File.Info <- "emreobesity.SSD.Info"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFO <- Open_SSD(File.SSD, File.Info)
FAM<-Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)
y <- FAM$Phenotype
obj <- SKAT_Null_Model(y ~ 1, data=NULL, out_type="D", type.Resampling="bootstrap", Adjustment=FALSE)
out.skat <- SKATBinary.SSD.All(SSD.INFO, obj, kernel = "linear.weighted",
                               method="Burden", method.bin="ER", weights.beta=c(1,25), weights=NULL, 
                               impute.method="random", r.corr=1, is_check_genotype=TRUE,
                               is_dosage = FALSE, missing_cutoff=1 , max_maf=1, estimate_MAF=1)

output.df = out$results
write.table(output.df, file = "./emreobesityresultsbinaryouttypeDemrehoca")

# 4 - SKATO emmax
library(SKAT)
Close_SSD()
File.Bed <- "emreobesity.bed"
File.Bim <- "emreobesity.bim"
File.Fam <- "emreobesity.fam"
File.SetID <- "emreobesitysetID.setID"
File.SSD <- "emreobesity.SSD"
File.Info <- "emreobesity.SSD.Info"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFO <- Open_SSD(File.SSD, File.Info)
File.Kin <- "./aa.hBN.kinf"
FAM<-Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)
y <- FAM$Phenotype
obj <- SKAT_NULL_emmaX (y ~ 1, data=NULL, K=NULL, Kin.File=File.Kin, ngrids=100, 
                        llim=-10, ulim=10, esp=1e-10, Is.GetEigenResult=FALSE)
out.skato.emmax <- SKAT.SSD.All(SSD.INFO, obj, kernel = "linear.weighted",
                                method="Burden", weights.beta=c(1,25),
                                weights=NULL, impute.method="bestguess", r.corr=0,
                                is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=1)

output.df = out$results
write.table(output.df, file = "./emreobesityresultsskatoemmaxemrehoca")