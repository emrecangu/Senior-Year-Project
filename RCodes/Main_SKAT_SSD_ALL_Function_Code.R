
## Binary SKAT Emre Hoca

## REPEAT 50Kmis20 Variant SKAT - Burden - SKATO qqplot


library(SKAT)
Close_SSD()
File.Bed <- "./obesity2/obesitycohort2.bed"
File.Bim <- "./obesity2/obesitycohort2.bim"
File.Fam <- "./obesity2/obesitycohort.fam"
File.SetID <- "./setID2"
File.SSD <- "./obesity/obesitycohort.SSD"
File.Info <- "./obesity/obesitycohort.SSD.info"
File.Kin <- "./obesity/emmax.hBN.kinf"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFOX50k <- Open_SSD(File.SSD, File.Info)
FAM<-Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)
y <- FAM$Phenotype
objx <- SKAT_NULL_emmaX (y ~ 1, data=NULL, K=NULL, Kin.File=File.Kin, ngrids=100, 
                         llim=-10, ulim=10, esp=1e-10, Is.GetEigenResult=FALSE)
out.skato.emmax <- SKAT.SSD.All(SSD.INFOX50k, objx, kernel = "linear.weighted",
                                method="SKATO", weights.beta=c(1,40),
                                weights=NULL, impute.method="bestguess", r.corr=0,
                                is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=1)
out.power <- Power_Logistic(Haplotypes = NULL, SNP.Location = NULL, SubRegion.Length=-1
                            , Prevalence=0.01, Case.Prop=0.5, Causal.Percent=100, Causal.MAF.Cutoff=0.03
                            , alpha =c(0.01,10^(-3),10^(-6)), N.Sample.ALL = 3824 * (1:10)
                            , Weight.Param=c(1,40), N.Sim=100, OR.Type = "Log"
                            , MaxOR=5, Negative.Percent=0)
Get_RequiredSampleSize(out.power, Power=0.8)
warnings()
output.df = out.skato.emmax$results
write.table(output.df, file = "./obesity2/resultsinorder2/SKATSSD_ALL-BN-LW-SKATO-40-bestguess-0-100")

Get_EffectiveNumberTest(out.skato.emmax$results$MAP, alpha=0.05)

# QQ plot
source('binaryskatemrehoca.R')
T <- read.table('output.df',header=TRUE)
pdf('')
qq.conf.beta(T$P)
dev.off()
QQPlot_Adj(out.skato.emmax$results$P.value, out.skato.emmax$results$MAP)
QQPlot_Adj(out.skato.emmax$results$P.value, out.skato.emmax$P.value.Resampling, main="QQ plot", ntry=500, confidence=0.95, Is.unadjsted=FALSE
           , Is.legend=TRUE, xlab="Expected Quantiles (-log10 P-values)"
           , ylab="Observed Quantiles (-log10 P-values)")



install.packages("fdrtool")

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
fn<-"./obesity2/resultsinorder2/SKATSSD_ALL-BN-LW-SKATO-40-bestguess-0-100";
write.table(t,file=fn,sep="\t",quote=F,row.names=F,col.names=T);
fastqq2 <- function(pvals, ...) { np <- length(pvals); thin.idx <- 1:np; thin.logp.exp <- -log10(thin.idx/(np+1)); thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]]); plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])), ...); abline(0, 1, col='gray', lty=2); thin.idx <- c((0.9)^(5:1), thin.idx); logp.cint.95 <- -log10(qbeta(0.95, thin.idx, np - thin.idx + 1)); logp.cint.05 <- -log10(qbeta(0.05, thin.idx, np - thin.idx + 1)); thin.logp.exp <- -log10(thin.idx/(np+1)); lines(thin.logp.exp, logp.cint.95, lty=2, col='red'); lines(thin.logp.exp, logp.cint.05, lty=2, col='red'); }
l <- 1; fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); mt<-0.05/length(t$P.value); abline(h=-log10(mt), lty=3); png(paste(fn,'QQ','png',sep='.')); fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); dev.off();


library(SKAT)
Close_SSD()
File.Bed <- "./obesity/obesitycohortMissing20.bed"
File.Bim <- "./obesity/obesitycohortMissing20.bim"
File.Fam <- "./obesity/obesitycohort.fam"
File.SetID <- "./obesity/obesitycohortsetID.setID"
File.SSD <- "./obesity/obesitycohortMissing20.SSD"
File.Info <- "./obesity/obesitycohortMissing20.SSD.info"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFOX50kmis20 <- Open_SSD(File.SSD, File.Info)
FAM<-Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)
y <- FAM$Phenotype
objx <- SKAT_NULL_emmaX (y ~ 1, data=NULL, K=NULL, Kin.File=File.Kin, ngrids=100, 
                         llim=-10, ulim=10, esp=1e-10, Is.GetEigenResult=FALSE)
out.skato.emmax <- SKAT.SSD.All(SSD.INFOX50kmis20, objx, kernel = "linear.weighted",
                                method="Burden", weights.beta=c(1,40),
                                weights=NULL, impute.method="bestguess", r.corr=0,
                                is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=1)

output.df = out.skat$results
write.table(output.df, file = "./obesity/resultsinorder1/50kmis20-SKATSSD_ALL-BN-LW-Burden-40-bestguess-0-100")

# QQ plot

QQPlot_Adj(out.skat$results$P.value, out.skat$results$MAP)

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
fn<-"./obesity/resultsinorder1/50kmis20-SKATSSD_ALL-BN-LW-Burden-40-bestguess-0-100";
write.table(t,file=fn,sep="\t",quote=F,row.names=F,col.names=T);
fastqq2 <- function(pvals, ...) { np <- length(pvals); thin.idx <- 1:np; thin.logp.exp <- -log10(thin.idx/(np+1)); thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]]); plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])), ...); abline(0, 1, col='gray', lty=2); thin.idx <- c((0.9)^(5:1), thin.idx); logp.cint.95 <- -log10(qbeta(0.95, thin.idx, np - thin.idx + 1)); logp.cint.05 <- -log10(qbeta(0.05, thin.idx, np - thin.idx + 1)); thin.logp.exp <- -log10(thin.idx/(np+1)); lines(thin.logp.exp, logp.cint.95, lty=2, col='red'); lines(thin.logp.exp, logp.cint.05, lty=2, col='red'); }
l <- 1; fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); mt<-0.05/length(t$P.value); abline(h=-log10(mt), lty=3); png(paste(fn,'QQ','png',sep='.')); fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); dev.off();

a