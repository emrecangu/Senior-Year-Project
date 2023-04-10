setwd("e:/bitirme_tezi/maincohort")

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
output.df = out$OUT.snp.mac
df <- data.frame(x1=x1, x2=x2, y=y)
write.table(output.df, file = "./emreobesityresultssnpmacbinaryouttypeD")
