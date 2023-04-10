
# 2 - Burden with SKAT-O weights and dichotomous Emre Hoca

library(SKAT)
Close_SSD()
File.Bed <- "obesitycohort.bed"
File.Bim <- "obesitycohort.bim"
File.Fam <- "obesitycohort.fam"
File.SetID <- "obesitycohortsetID.setID.temp"
File.SSD <- "obesitycohort.SSD"
File.Info <- "obesitycohort.SSD.info"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFO <- Open_SSD(File.SSD, File.Info)
FAM<-Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)
y <- FAM$Phenotype
obj <- SKAT_Null_Model(y ~ 1, data=NULL, out_type="D", type.Resampling="bootstrap", n.Resampling=0, Adjustment=TRUE)
out.skat <- SKAT.SSD.All(SSD.INFO, obj, kernel = "linear.weighted", 
                         method="Burden", weights.beta=c(1,25), weights=NULL, 
                         impute.method="random", r.corr=0, is_check_genotype=TRUE,
                         is_dosage = FALSE, missing_cutoff=0.9 , max_maf=1, estimate_MAF=1)
output.df = out.skat$results
write.table(output.df, file = "./resultsinorder/2obesitycohortnonmissingburdenskatodichotomous")

# QQ plot
Get_EffectiveNumberTest(out.skat$results$MAP, alpha=0.05)
QQPlot_Adj(out.skat$results$P.value, out.skat$results$MAP)

library(SKAT)
Close_SSD()
File.Bed <- "obesitycohortMissing20.bed"
File.Bim <- "obesitycohortMissing20.bim"
File.Fam <- "obesitycohort.fam"
File.SetID <- "obesitycohortsetID.setID"
File.SSD <- "obesitycohortMissing20.SSD"
File.Info <- "obesitycohortMissing20.SSD.info"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFO <- Open_SSD(File.SSD, File.Info)
FAM<-Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)
y <- FAM$Phenotype
obj <- SKAT_Null_Model(y ~ 1, data=NULL, out_type="D", type.Resampling="bootstrap", n.Resampling=0, Adjustment=TRUE)
out.skat <- SKAT.SSD.All(SSD.INFO, obj, kernel = "linear.weighted", 
                         method="Burden", weights.beta=c(1,25), weights=NULL, 
                         impute.method="random", r.corr=0, is_check_genotype=TRUE,
                         is_dosage = FALSE, missing_cutoff=0.9 , max_maf=1, estimate_MAF=1)
output.df = out.skat$results
write.table(output.df, file = "./resultsinorder/2obesitycohortmissing20burdenskatodichotomous")

# QQ plot

QQPlot_Adj(out.skat$results$P.value, out.skat$results$MAP)