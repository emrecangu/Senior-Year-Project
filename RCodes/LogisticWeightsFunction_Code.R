
## SKATBinary SKAT LW with logistic weights

## objx stands for emmax null model, when obj = skat null model

library(SKAT)
Close_SSD()
File.Bed <- "obesitycohort2.bed"
File.Bim <- "obesitycohort2.bim"
File.Fam <- "obesitycohort.fam"
File.SetID <- "setID2.temp"
File.SSD <- "obesitycohort2.SSD"
File.Info <- "obesitycohort2.SSD.info"
File.Kin <- "emmax.hBN.kinf"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFOX <- Open_SSD(File.SSD, File.Info)
FAM<-Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)
y <- FAM$Phenotype
objx <- SKAT_NULL_emmaX(y ~ 1, data=NULL, K=NULL, Kin.File=File.Kin, ngrids=100, 
                        llim=-10, ulim=10, esp=1e-10, Is.GetEigenResult=TRUE)

SSD.INFOX_num <- as.numeric(unlist(SSD.INFOX))
unlist(SSD.INFOX_num)
list(SSD.INFOX_num)

a <- Get_Genotypes_SSD(SSD.INFOX_num, Set_Index, is_ID = TRUE)
is.atomic(SSD.INFO_num)
data_ZX <- as.data.frame(t(SSD.INFOX_num))
data_ZX
weights <- Get_Logistic_Weights(data_ZX, par1 = 0.07, par2 = 150)

cool_function(SSD.INFO)
Open_SSD(File.SSD, File.Info)
out.skato.emmax  <- SKAT.SSD.All(SSD.INFOX, objx, kernel = "linear.weighted", method = "Burden",
                         is_dosage = FALSE, weights = weights, weights.beta = c(1,50), 
                         r.corr = 0, impute.method = "bestguess", is_check_genotype = TRUE )
Resampling_FWER(out.skato.emmax,FWER=0.5)
warnings()
output.df = out.skato.emmax$results
write.table(output.df, file="./resultsinorder/SKATSSD_ALL-EigenTrue-LW-SKATO-Burden-bestguess.txt", col.names=TRUE, row.names=FALSE)
QQPlot_Adj(out.skato.emmax$results$P.value, out.skato.emmax$results$MAP)

library(fdrtool)
qvals <- fdrtool(out.skato.emmax$results$P.value, statistic = "pvalue", cutoff.method="fndr", plot = FALSE)
out.skato.emmax$results$Q.value <- qvals$qval
out.skato.emmax$results <-out.skato.emmax$results[order(out.skato.emmax$results$P.value),]

t<-out.skato.emmax$results;
t<-t[with(t,order(P.value)),];
d<-dim(t)[1];
t$rank<-seq(1,d);
t$bonf<-0.05/d;
t$bonf.sig<-t$P.value<t$bonf;
fn<-"./resultsinorder/SKATSSD_ALL-EigenTrue-LW-SKATO-50-bestguess";
write.table(t,file=fn,sep="\t",quote=F,row.names=F,col.names=T);
fastqq2 <- function(pvals, ...) { np <- length(pvals); thin.idx <- 1:np; thin.logp.exp <- -log10(thin.idx/(np+1)); thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]]); plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])), ...); abline(0, 1, col='gray', lty=2); thin.idx <- c((0.9)^(5:1), thin.idx); logp.cint.95 <- -log10(qbeta(0.95, thin.idx, np - thin.idx + 1)); logp.cint.05 <- -log10(qbeta(0.05, thin.idx, np - thin.idx + 1)); thin.logp.exp <- -log10(thin.idx/(np+1)); lines(thin.logp.exp, logp.cint.95, lty=2, col='red'); lines(thin.logp.exp, logp.cint.05, lty=2, col='red'); }
l <- 1; fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); mt<-0.05/length(t$P.value); abline(h=-log10(mt), lty=3); png(paste(fn,'QQ','png',sep='.')); fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); dev.off();

a