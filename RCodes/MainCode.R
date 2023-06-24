library(SKAT)
Close_SSD()
File.Bed <- "Path/to/.bed"
File.Bim <- "Path/to/.bim"
File.Fam <- "Path/to/.fam"
File.SetID <- "Path/to/setID"
File.SSD <- "Path/to/.SSD"
File.Info <- "Path/to/.SSD.info"
File.Kin <- "Path/to/emmax.hBN.kinf"
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info, Is.FlipGenotype=TRUE)
SSD.INFO <- Open_SSD(File.SSD, File.Info)
FAM<-Read_Plink_FAM(File.Fam, Is.binary=TRUE, flag1=0)
y <- FAM$Phenotype

## Continous Traits - use C, Dichotomous Traits - use D, to obtain null model ##
obj <- SKAT_Null_Model(formula, data=NULL, out_type="D", n.Resampling=0
                            , type.Resampling="bootstrap", Adjustment=TRUE)

## Run the code below for non-kinship 
out.skat <- SKAT.SSD.All(SSD.INFO, obj, kernel = "linear or linear.weighted", 
                         method="Burden, SKAT and SKATO", weights.beta=c(1,25 or 40 or 50), weights=NULL, 
                         impute.method="bestguess or random or fixed", r.corr=0, is_check_genotype=TRUE,
                         is_dosage = FALSE, missing_cutoff=1 , max_maf=1, estimate_MAF=1)

## If the kinship model is used, use the below code for null model ##
objx <- SKAT_NULL_emmaX (y ~ 1, data=NULL, K=NULL, Kin.File=File.Kin, ngrids=100, 
                         llim=-10, ulim=10, esp=1e-10, Is.GetEigenResult=FALSE)

## Run the code below if kinship file is used
out.skato.emmax <- SKAT.SSD.All(SSD.INFO, objx, kernel = "linear.weighted",
                                method="SKATO", weights.beta=c(1,40),
                                weights=NULL, impute.method="bestguess", r.corr=0,
                                is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=1)

## To write the results, as a data frame, use the code below, it will generate tab-delimited text file ##
output.df = out.skat$results ## if you performed out.skato.emmax, change it ##
write.table(output.df, file = "path/to/write")

Get_EffectiveNumberTest(out.skato.emmax$results$MAP, alpha=0.05)

## QQ plot generation, adding bonf.sig to the data frame, ranking the associated genes
install.packages("fdrtool")

library(fdrtool)
qvals <- fdrtool(out.skat$results$P.value, statistic = "pvalue", cutoff.method="fndr", plot = FALSE)
out.skat$results$Q.value <- qvals$qval
out.skat$results <-out.skat$results[order(out.skat$results$P.value),]

t<-out.skat$results;
t<-t[with(t,order(P.value)),];
d<-dim(t)[1];
t$rank<-seq(1,d);
t$bonf<-0.05/d;
t$bonf.sig<-t$P.value<t$bonf;
fn<-"./obesity2/resultsinorder2/SKATSSD_ALL-BN-LW-SKATO-40-bestguess-0-100";
write.table(t,file=fn,sep="\t",quote=F,row.names=F,col.names=T);
fastqq2 <- function(pvals, ...) { np <- length(pvals); thin.idx <- 1:np; thin.logp.exp <- -log10(thin.idx/(np+1)); thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]]); plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])), ...); abline(0, 1, col='gray', lty=2); thin.idx <- c((0.9)^(5:1), thin.idx); logp.cint.95 <- -log10(qbeta(0.95, thin.idx, np - thin.idx + 1)); logp.cint.05 <- -log10(qbeta(0.05, thin.idx, np - thin.idx + 1)); thin.logp.exp <- -log10(thin.idx/(np+1)); lines(thin.logp.exp, logp.cint.95, lty=2, col='red'); lines(thin.logp.exp, logp.cint.05, lty=2, col='red'); }
l <- 1; fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); mt<-0.05/length(t$P.value); abline(h=-log10(mt), lty=3); png(paste(fn,'QQ','png',sep='.')); fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); dev.off();
