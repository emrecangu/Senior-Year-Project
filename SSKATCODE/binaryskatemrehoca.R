
## Binary SKAT Emre Hoca

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
obj <- SKAT_Null_Model(y ~ 1, data=NULL, out_type="D", type.Resampling="bootstrap", Adjustment=FALSE)
out.skat <- SKATBinary.SSD.All(SSD.INFO, obj, kernel = "linear.weighted",
                               method="Burden", method.bin="ER", weights.beta=c(1,25), weights=NULL, 
                               impute.method="random", r.corr=1, is_check_genotype=TRUE,
                               is_dosage = FALSE, missing_cutoff=1 , max_maf=1, estimate_MAF=1)

output.df = out.skat$results
write.table(output.df, file = "./resultsinorder/3obesitycohortbinaryskat")


# QQ plot

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
obj <- SKAT_Null_Model(y ~ 1, data=NULL, out_type="D", type.Resampling="bootstrap", Adjustment=FALSE)
out.skat <- SKATBinary.SSD.All(SSD.INFO, obj, kernel = "linear.weighted",
                               method="Burden", method.bin="ER", weights.beta=c(1,25), weights=NULL, 
                               impute.method="random", r.corr=1, is_check_genotype=TRUE,
                               is_dosage = FALSE, missing_cutoff=1 , max_maf=1, estimate_MAF=1)

output.df = out.skat$results
write.table(output.df, file = "./resultsinorder/3obesitycohortmissing20binaryskat")

# QQ plot

QQPlot_Adj(out.skat$results$P.value, out.skat$results$MAP)
