
## SKATBinary Robust LW with logistic weights

## objx stands for emmax null model, when obj = skat null model

## File.Cov = File_Cov

library(SKAT)
Close_SSD()
File.Bed <- "obesitycohort2.bed"
File.Bim <- "obesitycohort2.bim"
File.Fam <- "obesitycohort.fam"
File.Cov <- "cov5age.txt"
File.Kin <- "emmax.hBN.kinf"
File.SetID <- "setID2.temp"
File.SSD <- "obesitycohort2.SSD"
File.Info <- "obesitycohort2.SSD.info"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFO <- Open_SSD(File.SSD, File.Info)

FAM<-Read_Plink_FAM_Cov(File.Cov, File.Cov, Is.binary=TRUE, flag1=0, cov_header = TRUE)
FAM<-Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)

y <- FAM$Phenotype

obj <- SKAT_Null_Model(y ~ 1, data=NULL, out_type="D", Adjustment=TRUE, n.Resampling =  0, type.Resampling = "bootstrap")

SSD.INFO_num <- as.numeric(unlist(SSD.INFO))
unlist(SSD.INFO_num)
list(SSD.INFO_num)

a <- Get_Genotypes_SSD(SSD.INFO_num, Set_Index, is_ID = TRUE)
is.atomic(SSD.INFO_num)
data_ZX <- as.data.frame(t("SKATBinary_Robust_SSD_ALL-AGECOV-Linear-Burden-0505-50-bestguess"))
data_ZX
weights <- Get_Logistic_Weights(data_ZX, par1 = 0.5, par2 = 0.5)

write.table(weights,file = "./resultsinorder/weights00705")
out.skat <- SKATBinary_Robust.SSD.All(SSD.INFO, objx, kernel = "linear",method="Burden",
                                      weights = NULL, weights.beta = c(1,50),
                                      impute.method="bestguess", method.bin = "ER", r.corr  = NULL, is_check_genotype = TRUE,
                                      is_dosage = FALSE, missing_cutoff = 1, max_maf = 1, estimate_MAF,
                                      N.Resampling = 2* 10^6, seednum = 100, epsilon = 10^-6, SetID = NULL,
                                      obj.SNPWeight = NULL)

out.skat.binaryx <- SKATBinary_Robust.SSD.All(SSD.INFO, objx, kernel = "linear.weighted", method="SKAT"
                                            , r.corr=NULL, weights.beta=c(1,25), weights = NULL
                                            , impute.method = "bestguess", is_check_genotype=TRUE
                                            , is_dosage = FALSE, missing_cutoff=0.15, max_maf=1
                                            , estimate_MAF=1)

#out.skat <- ((SSD.INFO, objx, kernel = "linear.weighted", method = "SKATO",
#            is_dosage = FALSE, weights = NULL, weights.beta = c(1,50), 
#           r.corr = 0, impute.method = "bestguess", is_check_genotype = TRUE ))

out.b <- Power_Logistic_R(Haplotypes = NULL, SNP.Location = NULL, SubRegion.Length = -1, Prevalence = 0.01, Case.Prop = 0.5,
                          Causal.Percent = 5, Causal.MAF.Cutoff = 0.03, alpha = c(0.01, 10^(-3), 10^(-6)),N.Sample.ALL = 500 *(1:10))
Resampling_FWER(out.skat,FWER=0.5)
warnings()
output.df = out.skat.binaryx$results
write.table(output.df, file="./resultsinorder/SKATBinary_Robust_SSD_ALL-LW-Burden-25-bestguess-015-new.txt", col.names=TRUE, row.names=FALSE)
QQPlot_Adj(out.skat.binaryx$results$P.value, out.skat.binaryx$results$MAP)

out.skat <- data_ZX

library(fdrtool)
qvals <- fdrtool(out.skat.binary$results$P.value, statistic = "pvalue", cutoff.method="fndr", plot = FALSE)
out.skat.binary$results$Q.value <- qvals$qval
out.skat.binary$results <-out.skat.binary$results[order(out.skat.binary$results$P.value),]

t<-out.skat.binary$results;
t<-t[with(t,order(P.value)),];
d<-dim(t)[1];
t$rank<-seq(1,d);
t$bonf<-0.05/d;
t$bonf.sig<-t$P.value<t$bonf;
fn<-"./resultsinorder/SKAT_SSD_ALL-LW-Burden-50-bestguess.txt";
write.table(t,file=fn,sep="\t",quote=F,row.names=F,col.names=T);
fastqq2 <- function(pvals, ...) { np <- length(pvals); thin.idx <- 1:np; thin.logp.exp <- -log10(thin.idx/(np+1)); thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]]); plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])), ...); abline(0, 1, col='gray', lty=2); thin.idx <- c((0.9)^(5:1), thin.idx); logp.cint.95 <- -log10(qbeta(0.95, thin.idx, np - thin.idx + 1)); logp.cint.05 <- -log10(qbeta(0.05, thin.idx, np - thin.idx + 1)); thin.logp.exp <- -log10(thin.idx/(np+1)); lines(thin.logp.exp, logp.cint.95, lty=2, col='red'); lines(thin.logp.exp, logp.cint.05, lty=2, col='red'); }
l <- 1; fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); mt<-0.05/length(t$P.value); abline(h=-log10(mt), lty=3); png(paste(fn,'QQ','png',sep='.')); fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); dev.off();

a