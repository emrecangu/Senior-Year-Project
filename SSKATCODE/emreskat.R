setwd("e:/bitirme_tezi/maincohort")

library(SKAT)

File.Bed <- "emreobesity.bed"
File.Bim <- "emreobesity.bim"
File.Fam <- "emreobesity.fam"
File.SetID <- "emreobesitysetID.setID"
File.SSD <- "emreobesity.SSD"
File.Info <- "emreobesity.SSD.info"

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)

SSD.INFO<-Open_SSD(File.SSD, File.Info)
data(emreobesitysetID.setID)
obj<-SKAT_Null_Model(y.b ~ X, out_type="D", data=NULL)
SKAT(Z, obj, kernel = "linear.weighted")$p.value

SKAT_Null_Model(formula, data=File.SSD, out_type="C", n.Resampling=0
                , type.Resampling="bootstrap", Adjustment=TRUE)
SKAT.SSD.All(SSD.INFO, obj, obj.SNPWeight=NULL)
