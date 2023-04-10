############## Farvat #################
#### http://statgen.snu.ac.kr/wisard/?act=aa_gnfam
## 아래의 3개 SNP은 discovery set에서 minor allele이 0개
#3	chr3195475923
#3	chr3195507271
#3	chr3195513076
##

# cd /home2/wjkim/project/gastric_cancer/3.FARVAT

# make binary file & exclude SNPs with MAF of 0
# plink --file ../0.data/discovery --make-bed --exclude 0.data/maf0_SNPs.txt --out 0.data/discovery
# plink --file ../0.data/validation --make-bed --out 0.data/validation	# 0 maf 일단은 제거하지않음

# Make kinship correlation matrix
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT')
library(kinship2)
fam <- read.table('0.data/discovery.fam',head=F,stringsAsFactor=F)
ped <- pedigree(famid=fam$V1,id=fam$V2,dadid=fam$V3,momid=fam$V4,sex=fam$V5,affected=fam$V6,missid=0)
KM <- 2*as.matrix(kinship(ped))
write.table(KM,'0.data/discovery_KM.txt',row.names=F,quote=F)

fam <- read.table('0.data/validation.fam',head=F,stringsAsFactor=F)
ped <- pedigree(famid=fam$V1,id=fam$V2,dadid=fam$V3,momid=fam$V4,sex=fam$V5,affected=fam$V6,missid=0)
KM <- 2*as.matrix(kinship(ped))
write.table(KM,'0.data/validation_KM.txt',row.names=F,quote=F)

# HDGC group
setwd("/home2/wjkim/project/gastric_cancer/3.FARVAT")
fam <- read.table('0.data/discovery.fam',head=F)
pheno <- read.table('../0.data/new_phenofile.txt',head=T)

fam$id <- paste(fam$V1,fam$V2,sep='_')
pheno$id <- paste(pheno$FID,pheno$IID,sep='_')

HDGC0 <- pheno$id[which(pheno$HDGC==0)]
HDGC1 <- pheno$id[which(pheno$HDGC==1)]

fam0 <- fam1 <- fam
fam0$V6[fam0$id%in%HDGC1] <- -9
fam1$V6[fam1$id%in%HDGC0] <- -9

write.table(fam0,'0.data/discovery_HDGC0.fam',row.names=F,col.names=F,quote=F)
write.table(fam1,'0.data/discovery_HDGC1.fam',row.names=F,col.names=F,quote=F)

system("cp 0.data/discovery.bed 0.data/discovery_HDGC0.bed")
system("cp 0.data/discovery.bed 0.data/discovery_HDGC1.bed")
system("cp 0.data/discovery.bim 0.data/discovery_HDGC0.bim")
system("cp 0.data/discovery.bim 0.data/discovery_HDGC1.bim")

### Farvat
# Model 1
farvat --bed 0.data/discovery.bed --set 0.data/MUC4_SNPset_ver1.txt --genetest --cor 0.data/discovery_KM.txt --sampvar ../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC --mispheno NA --out 1.output/1.model1/res_cov --ignoreparent --skato
# Model 2
farvat --bed 0.data/discovery.bed --set 0.data/MUC4_SNPset_ver1.txt --genetest --cor 0.data/discovery_KM.txt --sampvar ../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/2.model2/res_cov --ignoreparent --skato
# Model 3 
farvat --bed 0.data/discovery_HDGC0.bed --set 0.data/MUC4_SNPset_ver1.txt --genetest --cor 0.data/discovery_KM.txt --sampvar ../0.data/new_phenofile.txt --cname Sex,Smoking01 --mispheno NA --out 1.output/3.model3/res_cov --ignoreparent --skato
# Model 4 
farvat --bed 0.data/discovery_HDGC1.bed --set 0.data/MUC4_SNPset_ver1.txt --genetest --cor 0.data/discovery_KM.txt --sampvar ../0.data/new_phenofile.txt --cname Sex,Smoking01 --mispheno NA --out 1.output/4.model4/res_cov --ignoreparent --skato
# Model 5 
farvat --bed 0.data/discovery_HDGC0.bed --set 0.data/MUC4_SNPset_ver1.txt --genetest --cor 0.data/discovery_KM.txt --sampvar ../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/5.model5/res_cov --ignoreparent --skato
# Model 6 
farvat --bed 0.data/discovery_HDGC1.bed --set 0.data/MUC4_SNPset_ver1.txt --genetest --cor 0.data/discovery_KM.txt --sampvar ../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/6.model6/res_cov --ignoreparent --skato



### Results
setwd("/home2/wjkim/project/gastric_cancer/3.FARVAT/1.output")
loadfile <- function(i){
  aa <- read.table(paste0(i,'.model',i,'/res_cov.gene.res'),head=T,stringsAsFactor=F)
  return(aa)
}

res <- lapply(1:6,loadfile)
res.I <- do.call(rbind,res)
res.II <- cbind(rep(paste0('Model',1:6),each=2),res.I)
colnames(res.II)[1] <- 'Model'
write.csv(res.II,'results_20180717.csv',row.names=F,quote=F)


################## Binomial data

### Basic Covariates : Sex, Smoking
## Exposure : MUC4 and MUC4_blue
# using HDGC as a covariate
Model 1 : basic covariates + HDGC
Model 2 : basic covariates + HDGC + MUC1

# Divide dataset into HDGC groups
Model 3 : basic covariates + HDGC_0
Model 4 : basic covariates + HDGC_1
Model 5 : basic covariates + HDGC_0 + MUC1
Model 6 : basic covariates + HDGC_1 + MUC1



################# Using whole genotype data

#### merge separated vcf files for each family into one file
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data')
vcflist <- system('ls 1.vcf',intern=T)
for(i in vcflist){
  # only SNP
  system(paste0('vcftools --vcf 1.vcf/',i,' --plink --remove-indels --out 2.plink/',gsub('\\.vcf','',i)))
  system(paste0('plink --file 2.plink/',gsub('\\.vcf','',i),' --make-bed --out 2.plink/',gsub('\\.vcf','',i)))
}

dis_fams <- c('family10','sub_family11','sub_family1','sub_family2','family3','family5','family6','family7','family8','sub_family9')
val_fams <- c('family12_new','family13_new','family14','family15')

### discovery
## overlapped & biallelic SNP
bim <- read.table(paste0('2.plink/',dis_fams[1],'.bim'),head=F,stringsAsFactor=F)
bim$V7 <- apply(bim[,c('V5','V6')],1,function(x) paste(x,collapse='/'))
bim <- bim[,-c(5,6)]
print(c(dis_fams[1],nrow(bim)))
for(i in dis_fams[-1]){
  bim <- cbind(bim,do.call(rbind,strsplit(bim[,'V7'],'/')))
  colnames(bim)[(ncol(bim)-1):ncol(bim)] <- c('A1','A2')
  bim.1 <- read.table(paste0('2.plink/',i,'.bim'),head=F,stringsAsFactor=F)
  mer <- merge(bim,bim.1,by='V2',sort=F,all=F)
  aa <- apply(mer[,c('A1','A2','V5','V6')],1,function(x) paste(paste0(unique(x)[unique(x)!=0],collapse='/'),sum(unique(x)!=0),sep=','))
  bb <- do.call(rbind,strsplit(aa,','))
  bb[which(bb[,2]==1),1] <- paste0('0/',bb[which(bb[,2]==1),1])
  valid.snp <- which(bb[,2]<3)
  bim <- cbind(V2=as.character(mer[valid.snp,1]),V7=bb[valid.snp,1])
  print(c(i,nrow(bim)))
}
# Finally, 107049 SNPs remained.
fin.bim <- bim.1[bim.1$V2%in%bim[,1],]
write.table(fin.bim,'2.plink/common_SNPs_discovery.txt',col.names=F,row.names=F,quote=F)

for(i in dis_fams){
  # only overlapped & biallelic SNP
  system(paste0('plink --bfile 2.plink/',i,' --extract 2.plink/common_SNPs_discovery.txt --make-bed --out 2.plink/',i,'_common'))
}

## merge list
merlist.dis <- paste0('2.plink/',dis_fams[-1],'_common')
write.table(merlist.dis,'2.plink/merlist_discovery.txt',col.names=F,row.names=F,quote=F)

## merge
system(paste0('plink --bfile 2.plink/',dis_fams[1],'_common --merge-list 2.plink/merlist_discovery.txt --make-bed --out 2.plink/discovery'))

### validation
## overlapped & biallelic SNP
bim <- read.table(paste0('2.plink/',val_fams[1],'.bim'),head=F,stringsAsFactor=F)
bim$V7 <- apply(bim[,c('V5','V6')],1,function(x) paste(x,collapse='/'))
bim <- bim[,-c(5,6)]
print(c(val_fams[1],nrow(bim)))
for(i in val_fams[-1]){
  bim <- cbind(bim,do.call(rbind,strsplit(bim[,'V7'],'/')))
  colnames(bim)[(ncol(bim)-1):ncol(bim)] <- c('A1','A2')
  bim.1 <- read.table(paste0('2.plink/',i,'.bim'),head=F,stringsAsFactor=F)
  mer <- merge(bim,bim.1,by='V2',sort=F,all=F)
  aa <- apply(mer[,c('A1','A2','V5','V6')],1,function(x) paste(paste0(unique(x)[unique(x)!=0],collapse='/'),sum(unique(x)!=0),sep=','))
  bb <- do.call(rbind,strsplit(aa,','))
  bb[which(bb[,2]==1),1] <- paste0('0/',bb[which(bb[,2]==1),1])
  valid.snp <- which(bb[,2]<3)
  bim <- cbind(V2=as.character(mer[valid.snp,1]),V7=bb[valid.snp,1])
  print(c(i,nrow(bim)))
}
# Finally, 107049 SNPs remained.
fin.bim <- bim.1[bim.1$V2%in%bim[,1],]
write.table(fin.bim,'2.plink/common_SNPs_validation.txt',col.names=F,row.names=F,quote=F)

for(i in val_fams){
  # only overlapped & biallelic SNP
  system(paste0('plink --bfile 2.plink/',i,' --extract 2.plink/common_SNPs_validation.txt --make-bed --out 2.plink/',i,'_common'))
}

## merge list
merlist.dis <- paste0('2.plink/',val_fams[-1],'_common')
write.table(merlist.dis,'2.plink/merlist_validation.txt',col.names=F,row.names=F,quote=F)

## merge
system(paste0('plink --bfile 2.plink/',val_fams[1],'_common --merge-list 2.plink/merlist_validation.txt --make-bed --out 2.plink/validation'))

## maf0 SNPs
bim.dis <- read.table("2.plink/discovery.bim",head=F,stringsAsFactor=F)
rm.dis <- bim.dis[bim.dis[,5]==0,]
write.table(rm.dis,'2.plink/maf0_discovery.txt',row.names=F,col.names=F,quote=F)

bim.val <- read.table("2.plink/validation.bim",head=F,stringsAsFactor=F)
rm.val <- bim.val[bim.val[,5]==0,]
write.table(rm.val,'2.plink/maf0_validation.txt',row.names=F,col.names=F,quote=F)

### QC
## discovery
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink
plink --bfile discovery --geno 0.05 --hwe 1e-5 --exclude maf0_discovery.txt --make-bed --out qc_snp_dis
# 107049(Total)-15889(geno<0.05)-3700(hwe<1e-5)-24435(0 maf)=63025
plink --bfile qc_snp_dis --mind 0.05 --make-bed --out clean_discovery
# No subjects were removed

## validation
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink
plink --bfile validation --geno 0.05 --hwe 1e-5 --exclude maf0_validation.txt --make-bed --out qc_snp_val
# 168973(Total)-43092(geno<0.05)-1088(hwe<1e-5)-46914(0 maf)=77879
plink --bfile qc_snp_val --mind 0.05 --make-bed --out clean_validation
# No subjects were removed
rm qc_snp*
  rm family*
  rm sub*
  
  ### Kinship coefficient
  setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/')
## discovery
KM <- read.table('0.data/discovery_KM.txt',head=F,stringsAsFactor=F)
id <- KM[1,]
fam <- read.table('2.whole_gene/0.data/2.plink/clean_discovery.fam',head=F
                  ,stringsAsFactor=F)
new.KM <- KM[match(fam[,2],id)+1,match(fam[,2],id)]
colnames(new.KM) <- fam[,2]
write.table(new.KM,'2.whole_gene/0.data/2.plink/discovery_KM.txt',row.names=F,quote=F)

## validation
KM <- read.table('0.data/validation_KM.txt',head=F,stringsAsFactor=F)
fam <- read.table('2.whole_gene/0.data/2.plink/clean_validation.fam',head=F
                  ,stringsAsFactor=F)
id <- gsub('-0','-',KM[1,])
id.fam <- gsub('-0','-',fam[,2])
new.KM <- KM[match(id.fam,id)+1,match(id.fam,id)]
colnames(new.KM) <- fam[,2]
write.table(new.KM,'2.whole_gene/0.data/2.plink/validation_KM.txt',row.names=F,quote=F)

### Set file
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data')
nm <- readLines("NM_position.txt")
get.geneinfo <- function(i,bim){
  if(i%%100==0)print(i)
  genes <- strsplit(nm[2*i-1],">|\t")[[1]]
  tc <- strsplit(nm[2*i],"\t")[[1]]
  chr <- ifelse(tc[1]=='chrX',23,ifelse(tc[1]=='chrY',24,gsub('chr','',tc[1])))
  strand <- tc[2]
  
  TC <- do.call(rbind,strsplit(tc[3:length(tc)],';'))
  command_str <- paste0('SNPs <- bim[bim$V1==',chr,'& (',paste(apply(TC[,1:2,drop=F],1,function(xx) paste0('(bim$V4 >=', paste(xx,collapse=' & bim$V4 <='),')')),collapse='|'),'),2]')
  eval(parse(text=command_str))
  snps <- ifelse(length(SNPs)==0,'',paste0(SNPs,collapse=';'))
  
  return(data.frame(NM_ID=genes[2],Gene_symbol=genes[3],CHR=chr,SNPs=snps,stringsAsFactors=F))
}

library(parallel)
# discovery
bim <- read.table('2.plink/clean_discovery.bim',head=F,stringsAsFactor=F)
res <- mclapply(1:(length(nm)/2),get.geneinfo,bim=bim,mc.cores=15)
res.dis <- do.call(rbind,res)
res.dis.indata <- res.dis[res.dis$SNPs!='',]	# 6497 genes & 8548 transcriptome

# validation
bim <- read.table('2.plink/clean_validation.bim',head=F,stringsAsFactor=F)
res <- mclapply(1:(length(nm)/2),get.geneinfo,bim=bim,mc.cores=15)
res.val <- do.call(rbind,res)
res.val.indata <- res.val[res.val$SNPs!='',]	# 6452 genes & 10632 transcriptome

## transcriptome
get.setfile.tc <- function(i,indata){
  tc <- indata[i,1]
  snps <- strsplit(indata[i,4],';')[[1]]
  return(data.frame(gene=tc,SNP=snps,stringsAsFactors=F))
}

# discovery
sf.dis <- mclapply(1:nrow(res.dis.indata),get.setfile.tc,indata=res.dis.indata,mc.cores=15)
setfile.dis <- do.call(rbind,sf.dis)
write.table(setfile.dis,'2.plink/discovery_tc.setfile',row.names=F,col.names=F,quote=F)

# validataion
sf.val <- mclapply(1:nrow(res.val.indata),get.setfile.tc,indata=res.val.indata,mc.cores=15)
setfile.val <- do.call(rbind,sf.val)
write.table(setfile.val,'2.plink/validation_tc.setfile',row.names=F,col.names=F,quote=F)


## gene
get.setfile.gene <- function(gene,indata){
  tc <- indata[indata[,2]==gene,]
  snps <- unique(do.call(c,strsplit(tc[,4],';')))
  return(data.frame(gene=gene,SNP=snps,stringsAsFactors=F))
}

# discovery
sf.dis <- mclapply(unique(res.dis.indata[,2]),get.setfile.gene,indata=res.dis.indata,mc.cores=15)
setfile.dis <- do.call(rbind,sf.dis)
write.table(setfile.dis,'2.plink/discovery_gene.setfile',row.names=F,col.names=F,quote=F)

# validataion
sf.val <- mclapply(unique(res.val.indata[,2]),get.setfile.gene,indata=res.val.indata,mc.cores=15)
setfile.val <- do.call(rbind,sf.val)
write.table(setfile.val,'2.plink/validation_gene.setfile',row.names=F,col.names=F,quote=F)

#### HDGC group, family id
setwd("/home2/wjkim/project/gastric_cancer/3.FARVAT")
pheno <- read.table('../0.data/new_phenofile.txt',head=T)

# discovery
fam <- read.table('2.whole_gene/0.data/2.plink/clean_discovery.fam',head=F,stringsAsFactor=F)
fam[,1] <- pheno[match(fam$V2,pheno$IID),'FID']
fam[,6] <- pheno[match(fam$V2,pheno$IID),'GC']

HDGC0 <- pheno$IID[which(pheno$HDGC==0)]
HDGC1 <- pheno$IID[which(pheno$HDGC==1)]

fam0 <- fam1 <- fam
fam0$V6[fam0$V2%in%HDGC1] <- NA
fam1$V6[fam1$V2%in%HDGC0] <- NA

setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/')
write.table(fam,'clean_discovery.fam',row.names=F,col.names=F,quote=F)
write.table(fam0,'clean_discovery_HDGC0.fam',row.names=F,col.names=F,quote=F)
write.table(fam1,'clean_discovery_HDGC1.fam',row.names=F,col.names=F,quote=F)

system("cp clean_discovery.bed clean_discovery_HDGC0.bed")
system("cp clean_discovery.bed clean_discovery_HDGC1.bed")
system("cp clean_discovery.bim clean_discovery_HDGC0.bim")
system("cp clean_discovery.bim clean_discovery_HDGC1.bim")

# validation
setwd("/home2/wjkim/project/gastric_cancer/3.FARVAT")
fam <- read.table('2.whole_gene/0.data/2.plink/clean_validation.fam',head=F,stringsAsFactor=F)
fam$V2 <- gsub('-0','-',fam$V2)
fam[,1] <- pheno[match(fam$V2,pheno$IID),'FID']
fam[,6] <- pheno[match(fam$V2,pheno$IID),'GC']

HDGC0 <- pheno$IID[which(pheno$HDGC==0)]
HDGC1 <- pheno$IID[which(pheno$HDGC==1)]

fam0 <- fam1 <- fam
fam0$V6[fam$V2%in%HDGC1] <- NA
fam1$V6[fam$V2%in%HDGC0] <- NA

setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/')
write.table(fam,'clean_validation.fam',row.names=F,col.names=F,quote=F)
write.table(fam0,'clean_validation_HDGC0.fam',row.names=F,col.names=F,quote=F)
write.table(fam1,'clean_validation_HDGC1.fam',row.names=F,col.names=F,quote=F)

system("cp clean_validation.bed clean_validation_HDGC0.bed")
system("cp clean_validation.bed clean_validation_HDGC1.bed")
system("cp clean_validation.bim clean_validation_HDGC0.bim")
system("cp clean_validation.bim clean_validation_HDGC1.bim")


############ FARVAT
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/
  ### discovery 
  ## transcriptome
  # Model 2 (HDGC : all, MUC1 included)
  farvat --bed 0.data/2.plink/clean_discovery.bed --set 0.data/2.plink/discovery_tc.setfile --genetest --cor 0.data/2.plink/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/1.discovery/1.tc/2.model2/res_cov --ignoreparent --skato --thread 15
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/clean_discovery_HDGC0.bed --set 0.data/2.plink/discovery_tc.setfile --genetest --cor 0.data/2.plink/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/1.discovery/1.tc/5.model5/res_cov --ignoreparent --skato --thread 15
# Model 6 (HDGC : 1, MUC1 included)
farvat --bed 0.data/2.plink/clean_discovery_HDGC1.bed --set 0.data/2.plink/discovery_tc.setfile --genetest --cor 0.data/2.plink/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/1.discovery/1.tc/6.model6/res_cov --ignoreparent --skato --thread 15

## gene
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/clean_discovery.bed --set 0.data/2.plink/discovery_gene.setfile --genetest --cor 0.data/2.plink/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/1.discovery/2.gene/2.model2/res_cov --ignoreparent --skato --thread 15
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/clean_discovery_HDGC0.bed --set 0.data/2.plink/discovery_gene.setfile --genetest --cor 0.data/2.plink/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/1.discovery/2.gene/5.model5/res_cov --ignoreparent --skato --thread 15
# Model 6 (HDGC : 1, MUC1 included)
farvat --bed 0.data/2.plink/clean_discovery_HDGC1.bed --set 0.data/2.plink/discovery_gene.setfile --genetest --cor 0.data/2.plink/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/1.discovery/2.gene/6.model6/res_cov --ignoreparent --skato --thread 15


### validation 
## transcriptome
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/clean_validation.bed --set 0.data/2.plink/validation_tc.setfile --genetest --cor 0.data/2.plink/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/2.validation/1.tc/2.model2/res_cov --ignoreparent --skato --thread 15
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/clean_validation_HDGC0.bed --set 0.data/2.plink/validation_tc.setfile --genetest --cor 0.data/2.plink/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/2.validation/1.tc/5.model5/res_cov --ignoreparent --skato --thread 15
# Model 6 (HDGC : 1, MUC1 included) 너무 적어서 불가 4명
farvat --bed 0.data/2.plink/clean_validation_HDGC1.bed --set 0.data/2.plink/validation_tc.setfile --genetest --cor 0.data/2.plink/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/2.validation/1.tc/6.model6/res_cov --ignoreparent --skato --thread 15

## gene
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/clean_validation.bed --set 0.data/2.plink/validation_gene.setfile --genetest --cor 0.data/2.plink/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/2.validation/2.gene/2.model2/res_cov --ignoreparent --skato --thread 15
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/clean_validation_HDGC0.bed --set 0.data/2.plink/validation_gene.setfile --genetest --cor 0.data/2.plink/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/2.validation/2.gene/5.model5/res_cov --ignoreparent --skato --thread 15
# Model 6 (HDGC : 1, MUC1 included) 너무 적어서 불가 4명
farvat --bed 0.data/2.plink/clean_validation_HDGC1.bed --set 0.data/2.plink/validation_gene.setfile --genetest --cor 0.data/2.plink/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/2.validation/2.gene/6.model6/res_cov --ignoreparent --skato --thread 15

### Summary
## get plotdata
Mode <- c('tc','gene')
get.plotdata <- function(ii,jj){
  aa <- read.table(paste0(ii,'.',Mode[ii],'/',jj,'.model',jj,'/res_cov.gene.res'),head=T,stringsAsFactor=F)
  plotdata <- na.omit(aa[,c(1,2,ncol(aa),3,4,5,7,9,10)])
  colnames(plotdata)[1:3] <- c('CHR','GENE','p.val')
  plotdata <- plotdata[with(plotdata,order(CHR,START)),]
  write.table(plotdata,paste0(ii,".",Mode[ii],"/SKATO_",Mode[ii],"_model",jj,"_plotdata.txt"),row.names=F,quote=F)
}

# discovery
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/1.discovery')
for(ii in 1:2) for(jj in c(2,5,6)) get.plotdata(ii,jj)

# validation
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/2.validation')
for(ii in 1:2) for(jj in c(2,5)) get.plotdata(ii,jj)

## results
Mode <- c('tc','gene')
get.results <- function(ii,jj){
  plotdata <- read.table(paste0(ii,".",Mode[ii],"//SKATO_",Mode[ii],"_model",jj,"_plotdata.txt"),head=T,stringsAsFactor=F)
  sig.level <- 0.05/nrow(plotdata)
  
  ## significant results
  sig.res <- plotdata[plotdata$p.val<=sig.level,]
  
  ## QQ plot
  png(paste(ii,".",Mode[ii],"//SKATO_",Mode[ii],"_model",jj,"_QQ.png",sep=""), width=900, height=1000, units = "px", bg = "white", res = 200)
  
  logit <- plotdata
  p.val <- logit[,'p.val']
  y <- -log(p.val,10)
  v <- -log10(0.05/sum(is.finite(y)))
  o.y <- sort(y[is.finite(y)],decreasing=T)
  xx<- (1:length(o.y))/length(o.y)
  x <- -log(xx,10)
  maxy <- max(o.y)
  miny <- min(o.y)
  plot(x, x, xlim=c(miny,maxy), ylim=c(miny,maxy), type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")))
  N=length(o.y)
  c95 <- rep(0,N)
  c05 <- rep(0,N)
  for(i in 1:N){
    c95[i] <- qbeta(0.95,i,N-i+1)
    c05[i] <- qbeta(0.05,i,N-i+1)
  }
  abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
  polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
  abline(a=0,b=1,col="black", lty=1)
  points(x, o.y, cex = .5, col = "dark red")
  
  dev.off()
  
  
  ## MHT plot
  png(paste(ii,".",Mode[ii],"//SKATO_",Mode[ii],"_model",jj,"_MHT.png",sep=""), width=2000, height=1000, units = "px", bg = "white", res = 200)
  
  in_data <- plotdata[,c('CHR','GENE','p.val')]
  
  kare <- table(in_data$CHR)
  CM <- cumsum(kare)
  n.markers <- sum(kare)
  n.chr <- length(kare)
  
  pos.1 <- 1
  pos <- c()
  for(i in 1:length(kare)){
    diff <- max(kare) - kare[i]
    d <- ifelse(diff%%2==0,diff/2,(diff-1)/2)
    pos <- c(pos,(pos.1+d):(pos.1+d+kare[i]-1))
    pos.1 <- ifelse(diff%%2==0,pos[length(pos)]+d+1,pos[length(pos)]+d+2)
  }
  
  y <- -log(in_data$p.val,10)
  max.y <- ifelse(max(y)<6,6,max(y))
  
  colors <- rainbow(n.chr)
  par(xaxt = "n", yaxt = "n")
  plot(pos, y, type = "n", xlab = "Chromosome", ylab = expression(paste(-log[10],"(p-value)")), axes = FALSE, ylim=c(0,max.y+1),cex.main=1.5, cex.lab=1.2)
  axis(1, tick = FALSE)
  axis(2, tick = FALSE)
  for (i in 1:n.chr) {
    u <- CM[i]
    l <- CM[i] - kare[i] + 1
    cat("Plotting points ", l, "-", u, "\n")
    chr <- l:u
    tran.y <- y[chr]
    points(pos[chr], tran.y, col = colors[i],pch=20,cex=0.5)
  }
  par(xaxt = "s", yaxt = "s")
  AT <- cumsum(c(0,rep(max(kare),length(kare)-1)))+max(kare)/2
  LB <- 1:n.chr
  LB[LB==23] <- 'X'
  LB[LB==24] <- 'Y'
  axis(1, at = AT, labels = LB)
  if(max.y+1 < 10) {
    axis(2, at = 0:(max.y+1))
  } else {
    axis(2, at = seq(0,max.y+1,by=2))
  }
  abline(h=-log10(sig.level),lty=2,col='gray')
  dev.off()
}

# discovery
setwd('C:\\Users\\김원지\\Desktop\\GC\\farvat\\whole_gene\\1.discovery')	
for(ii in 1:2) for(jj in c(2,5,6)) get.results(ii,jj)


cd ~/project/gastric_cancer/0.data
head -n 1 new_phenofile.txt > head.txt

<R>
  setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/1.vcf')

head <- readLines('~/project/gastric_cancer/0.data/head.txt')
head.1 <- strsplit(head,' ')[[1]]
snps <- head.1[grep('chr',head.1)]
pos <- gsub['chr3','',snps)
vcfs <- system('ls',intern=T)

bb <- function(p){
  aa <- sapply(vcfs,function(vcf){
    res <- system(paste0('grep ',p,' ',vcf),intern=T)
    ifelse(length(res)==0, return(NULL), return(res))
  })
  return(aa)
}



################# Using whole genotype data ------ monomorphic variants도 고려

#### Remove indel & multi-alleic
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/1.vcf')
dis_fams <- c('family10','sub_family11','sub_family1','sub_family2','family3','family5','family6','family7','family8','sub_family9')
val_fams <- c('family12_new','family13_new','family14','family15')
vcflist <- paste0(c(dis_fams,val_fams),".vcf")

for(i in vcflist){
  # only bi-allelic SNP
  system(paste0('vcftools --vcf ',i,' --recode --remove-indels --max-alleles 2 --out biSNP_',gsub('\\.vcf','',i)))
}

## allele 정보 추출
for(i in vcflist){
  system(paste0("awk -F '\t' '{if($1!~/#/) print $1,$2,$4,$5}' biSNP_",gsub('\\.vcf','',i),".recode.vcf > SNPlist_",gsub('\\.vcf','',i),".txt"))
}

varinfo <- read.table(paste0("SNPlist_",gsub('\\.vcf','',vcflist[1]),".txt"),head=F,stringsAsFactor=F)
varinfo$chrpos <- gsub('[ ]+','',apply(varinfo[,1:2],1,function(x) paste(x,collapse=':')))
for(i in vcflist[-1]){
  varinfo.1 <- read.table(paste0("SNPlist_",gsub('\\.vcf','',i),".txt"),head=F,stringsAsFactor=F)
  varinfo.1$chrpos <- gsub('[ ]+','',apply(varinfo.1[,1:2],1,function(x) paste(x,collapse=':')))
  mer <- merge(varinfo,varinfo.1,by='chrpos',all=F,sort=F)
  if(sum(mer[,4]!=mer[,8])==0){
    if(sum(mer[,5]!=mer[,9])!=0){
      rm.list <- mer[mer[,5]!=mer[,9],1]
    } else {
      rm.list <- NULL
    }
  } else {
    rm.list <- mer[mer[,4]!=mer[,8],1]
  }
  browser()
  ## exclude multi-allelic SNP
  varinfo.2 <- rbind(varinfo,varinfo.1)
  varinfo.2 <- varinfo.2[!duplicated(varinfo.2[,'chrpos']),]
  if(!is.null(rm.list)) varinfo.2 <- varinfo.2[!varinfo.2[,'chrpos']%in%rm.list,]
  varinfo <- varinfo.2
}
# 이렇게 되면 아직 multi SNP이 남아있게 됨. (rm.list로 제거 된 다음 또 나오면 그때는 제거되지 않기 때문. 당장 코드수정이 어려우니, 뒤에서 missnp이용하여 제거하기로 하자.)
colnames(varinfo) <- c('CHR','POSITION','REF','ALT','SNP')
varinfo <- varinfo[,c('CHR','SNP','POSITION','REF','ALT')]
varinfo[grep('GL|hs',varinfo[,1]),1] <- 0
#xchr
xchr <- which(varinfo[,1]=='X')
varinfo[xchr,1] <- 23
varinfo[xchr,2] <- gsub('X',23,varinfo[xchr,2])
#ychr
ychr <- which(varinfo[,1]=='Y')
varinfo[ychr,1] <- 24
varinfo[ychr,2] <- gsub('Y',24,varinfo[ychr,2])
#MTchr
mtchr <- which(varinfo[,1]=='MT')
varinfo[mtchr,1] <- 26
varinfo[mtchr,2] <- gsub('MT',26,varinfo[mtchr,2])
write.table(varinfo,'All_biSNP.list',row.names=F,col.names=F,quote=F)

#### plink format 변환
### only biSNP
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/1.vcf')
vcflist <- system('ls *.recode.vcf',intern=T)
for(i in vcflist){
  # only bi-allelic SNP
  system(paste0('vcftools --vcf ',i,' --plink --out ../2.plink/',gsub('\\.recode\\.vcf','',i)))
  system(paste0('plink --file ../2.plink/',gsub('\\.recode\\.vcf','',i),' --extract All_biSNP.list --make-bed --out ../2.plink/final_',gsub('\\.recode\\.vcf','',i)))
}

#### monomorphic SNP 
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink')
files <- gsub('\\.log','',system('ls final_biSNP*.log',intern=T))
biSNPs <- read.table('../1.vcf/All_biSNP.list',head=F,stringsAsFactor=F)

get.newGT <- function(kk){
  fam <- read.table(paste0(files[kk],'.fam'),head=F,stringsAsFactor=F)
  bim <- read.table(paste0(files[kk],'.bim'),head=F,stringsAsFactor=F)
  
  new <- biSNPs[!biSNPs[,2]%in%bim[,2],]
  get.GT <- function(i,nn,GT){
    if(i==0){
      return(fam)
    } else {
      if(i%%10000==0) print(i)
      res <- matrix(GT[i,4],nrow=nn,ncol=2)
      return(res)
    }
  }
  for(chr in unique(new$V1)){
    new.chr <- new[new$V1==chr,]
    res.newGT <- lapply(0:nrow(new.chr),get.GT,nn=nrow(fam),GT=new.chr)
    ped <- do.call(cbind,res.newGT)
    map <- data.frame(new.chr[,c(1,2)],0,new.chr[,3])
    
    write.table(ped,paste0('chr',chr,'_',files[kk],'.ped'),row.names=F,col.names=F,quote=F)
    write.table(map,paste0('chr',chr,'_',files[kk],'.map'),row.names=F,col.names=F,quote=F)
  }
  merfiles <- matrix(paste0(rep(paste0('chr',unique(new$V1),'_',files[kk]),each=2),c('.ped','.map')),ncol=2,byrow=T)
  write.table(merfiles,paste0('newGT_',files[kk],'.txt'),col.names=F,row.names=F,quote=F)
  system(paste0('plink --bfile ',files[kk],' --merge-list newGT_',files[kk],'.txt --make-bed --out whole/whole_',files[kk]))	
}

library(parallel)
RES <- mclapply(1:length(files),get.newGT,mc.cores=14)


#### merge separated vcf files for each family into one file
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink')
dis_fams <- c('family10','sub_family11','sub_family1','sub_family2','family3','family5','family6','family7','family8','sub_family9')
val_fams <- c('family12_new','family13_new','family14','family15')

### discovery
## merge list
merlist.dis <- paste0('whole/whole_final_biSNP_',dis_fams[-1])
write.table(merlist.dis,'whole/whole_merlist_discovery.txt',col.names=F,row.names=F,quote=F)

## merge
system(paste0('plink --bfile whole/whole_final_biSNP_',dis_fams[1],' --merge-list whole/whole_merlist_discovery.txt --make-bed --out whole/whole_discovery'))
for(i in dis_fams){
  system(paste0('plink --bfile whole/whole_final_biSNP_',i,' --exclude whole/whole_discovery-merge.missnp --make-bed --out whole/whole_final_realbiSNP_',i))
}
merlist.dis <- paste0('whole/whole_final_realbiSNP_',dis_fams[-1])
write.table(merlist.dis,'whole/whole_merlist_discovery_realbiSNP.txt',col.names=F,row.names=F,quote=F)
system(paste0('plink --bfile whole/whole_final_realbiSNP_',dis_fams[1],' --merge-list whole/whole_merlist_discovery_realbiSNP.txt --make-bed --out whole/whole_discovery'))
# Finally, 4139078 SNPs remained.

### validation
## merge list
merlist.val <- paste0('whole/whole_final_biSNP_',val_fams[-1])
write.table(merlist.val,'whole/whole_merlist_validation.txt',col.names=F,row.names=F,quote=F)

## merge
system(paste0('plink --bfile whole/whole_final_biSNP_',val_fams[1],' --merge-list whole/whole_merlist_validation.txt --make-bed --out whole/whole_validation'))
for(i in val_fams){
  system(paste0('plink --bfile whole/whole_final_biSNP_',i,' --exclude whole/whole_validation-merge.missnp --make-bed --out whole/whole_final_realbiSNP_',i))
}
merlist.val <- paste0('whole/whole_final_realbiSNP_',val_fams[-1])
write.table(merlist.val,'whole/whole_merlist_validation_realbiSNP.txt',col.names=F,row.names=F,quote=F)
system(paste0('plink --bfile whole/whole_final_realbiSNP_',val_fams[1],' --merge-list whole/whole_merlist_validation_realbiSNP.txt --make-bed --out whole/whole_validation'))
# Finally, 4140227 SNPs remained.


## maf0 SNPs
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole')

bim.dis <- read.table("whole_discovery.bim",head=F,stringsAsFactor=F)
rm.dis <- bim.dis[bim.dis[,5]==0,]	# 822634 SNPs
write.table(rm.dis,'maf0_discovery.txt',row.names=F,col.names=F,quote=F)

bim.val <- read.table("whole_validation.bim",head=F,stringsAsFactor=F)
rm.val <- bim.val[bim.val[,5]==0,]	# 1994820 SNPs
write.table(rm.val,'maf0_validation.txt',row.names=F,col.names=F,quote=F)

### QC
## discovery
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole
plink --bfile whole_discovery --geno 0.05 --hwe 1e-5 --exclude maf0_discovery.txt --make-bed --out qc_snp_dis
# 4139078(Total)-2329672(geno<0.05)-7741(hwe<1e-5)-822634(0 maf)=979031
plink --bfile qc_snp_dis --mind 0.05 --make-bed --out clean_discovery
# No subjects were removed

## validation
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole
plink --bfile whole_validation --geno 0.05 --hwe 1e-5 --exclude maf0_validation.txt --make-bed --out qc_snp_val
# 4140227(Total)-1748457(geno<0.05)-1982(hwe<1e-5)-1994820(0 maf)=394968
plink --bfile qc_snp_val --mind 0.05 --make-bed --out clean_validation
# No subjects were removed
rm qc_snp*
  
  ### Kinship coefficient
  setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/')
## discovery
KM <- read.table('0.data/discovery_KM.txt',head=F,stringsAsFactor=F)
id <- KM[1,]
fam <- read.table('2.whole_gene/0.data/2.plink/whole/clean_discovery.fam',head=F,stringsAsFactor=F)
new.KM <- KM[match(fam[,2],id)+1,match(fam[,2],id)]
colnames(new.KM) <- fam[,2]
write.table(new.KM,'2.whole_gene/0.data/2.plink/whole/discovery_KM.txt',col.names=F,row.names=F,quote=F)


## validation
KM <- read.table('0.data/validation_KM.txt',head=F,stringsAsFactor=F)
fam <- read.table('2.whole_gene/0.data/2.plink/whole/clean_validation.fam',head=F
                  ,stringsAsFactor=F)
id <- gsub('-0','-',KM[1,])
id.fam <- gsub('-0','-',fam[,2])
new.KM <- KM[match(id.fam,id)+1,match(id.fam,id)]
colnames(new.KM) <- fam[,2]
write.table(new.KM,'2.whole_gene/0.data/2.plink/whole/validation_KM.txt',col.names=F,row.names=F,quote=F)

### Set file
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data')
nm <- readLines("NM_position.txt")
get.geneinfo <- function(i,bim){
  if(i%%100==0)print(i)
  genes <- strsplit(nm[2*i-1],">|\t")[[1]]
  tc <- strsplit(nm[2*i],"\t")[[1]]
  chr <- ifelse(tc[1]=='chrX',23,ifelse(tc[1]=='chrY',24,gsub('chr','',tc[1])))
  strand <- tc[2]
  
  TC <- do.call(rbind,strsplit(tc[3:length(tc)],';'))
  command_str <- paste0('SNPs <- bim[bim$V1==',chr,'& (',paste(apply(TC[,1:2,drop=F],1,function(xx) paste0('(bim$V4 >=', paste(xx,collapse=' & bim$V4 <='),')')),collapse='|'),'),2]')
  eval(parse(text=command_str))
  snps <- ifelse(length(SNPs)==0,'',paste0(SNPs,collapse=';'))
  
  return(data.frame(NM_ID=genes[2],Gene_symbol=genes[3],CHR=chr,SNPs=snps,stringsAsFactors=F))
}

library(parallel)
# discovery
bim <- read.table('2.plink/whole/clean_discovery.bim',head=F,stringsAsFactor=F)
res <- mclapply(1:(length(nm)/2),get.geneinfo,bim=bim,mc.cores=14)
res.dis <- do.call(rbind,res)
res.dis.indata <- res.dis[grep(';',res.dis$SNPs),]	# 10270 genes & 16804 transcriptome

# validation
bim <- read.table('2.plink/whole/clean_validation.bim',head=F,stringsAsFactor=F)
res <- mclapply(1:(length(nm)/2),get.geneinfo,bim=bim,mc.cores=15)
res.val <- do.call(rbind,res)
res.val.indata <- res.val[grep(';',res.val$SNPs),]	# 8582 genes & 13991 transcriptome

## transcriptome
get.setfile.tc <- function(i,indata){
  tc <- indata[i,1]
  snps <- strsplit(indata[i,4],';')[[1]]
  return(data.frame(gene=tc,SNP=snps,stringsAsFactors=F))
}

# discovery
sf.dis <- mclapply(1:nrow(res.dis.indata),get.setfile.tc,indata=res.dis.indata,mc.cores=15)
setfile.dis <- do.call(rbind,sf.dis)
setfile.dis <- setfile.dis[order(setfile.dis[,1]),]
write.table(setfile.dis,'discovery_tc.setfile',row.names=F,col.names=F,quote=F)

# validataion
sf.val <- mclapply(1:nrow(res.val.indata),get.setfile.tc,indata=res.val.indata,mc.cores=15)
setfile.val <- do.call(rbind,sf.val)
setfile.val <- setfile.val[order(setfile.val[,1]),]
write.table(setfile.val,'validation_tc.setfile',row.names=F,col.names=F,quote=F)


## gene
get.setfile.gene <- function(gene,indata){
  tc <- indata[indata[,2]==gene,]
  snps <- unique(do.call(c,strsplit(tc[,4],';')))
  return(data.frame(gene=gene,SNP=snps,stringsAsFactors=F))
}

# discovery
sf.dis <- mclapply(unique(res.dis.indata[,2]),get.setfile.gene,indata=res.dis.indata,mc.cores=15)
setfile.dis <- do.call(rbind,sf.dis)
setfile.dis <- setfile.dis[order(setfile.dis[,1]),]
write.table(setfile.dis,'discovery_gene.setfile',row.names=F,col.names=F,quote=F)

# validataion
sf.val <- mclapply(unique(res.val.indata[,2]),get.setfile.gene,indata=res.val.indata,mc.cores=15)
setfile.val <- do.call(rbind,sf.val)
setfile.val <- setfile.val[order(setfile.val[,1]),]
write.table(setfile.val,'validation_gene.setfile',row.names=F,col.names=F,quote=F)

#### HDGC group, family id
setwd("/home2/wjkim/project/gastric_cancer/3.FARVAT")
pheno <- read.table('../0.data/new_phenofile.txt',head=T)

# discovery
fam <- read.table('2.whole_gene/0.data/2.plink/whole/clean_discovery.fam',head=F,stringsAsFactor=F)
fam[,1] <- pheno[match(fam$V2,pheno$IID),'FID']
fam[,6] <- pheno[match(fam$V2,pheno$IID),'GC']

HDGC0 <- pheno$IID[which(pheno$HDGC==0)]
HDGC1 <- pheno$IID[which(pheno$HDGC==1)]

fam0 <- fam1 <- fam
fam0$V6[fam0$V2%in%HDGC1] <- NA
fam1$V6[fam1$V2%in%HDGC0] <- NA

setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole/')
write.table(fam,'clean_discovery.fam',row.names=F,col.names=F,quote=F)
write.table(fam0,'clean_discovery_HDGC0.fam',row.names=F,col.names=F,quote=F)
write.table(fam1,'clean_discovery_HDGC1.fam',row.names=F,col.names=F,quote=F)

system("cp clean_discovery.bed clean_discovery_HDGC0.bed")
system("cp clean_discovery.bed clean_discovery_HDGC1.bed")
system("cp clean_discovery.bim clean_discovery_HDGC0.bim")
system("cp clean_discovery.bim clean_discovery_HDGC1.bim")

# validation
setwd("/home2/wjkim/project/gastric_cancer/3.FARVAT")
fam <- read.table('2.whole_gene/0.data/2.plink/whole/clean_validation.fam',head=F,stringsAsFactor=F)
fam$V2 <- gsub('-0','-',fam$V2)
fam[,1] <- pheno[match(fam$V2,pheno$IID),'FID']
fam[,6] <- pheno[match(fam$V2,pheno$IID),'GC']

HDGC0 <- pheno$IID[which(pheno$HDGC==0)]
HDGC1 <- pheno$IID[which(pheno$HDGC==1)]

fam0 <- fam1 <- fam
fam0$V6[fam$V2%in%HDGC1] <- NA
fam1$V6[fam$V2%in%HDGC0] <- NA

setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole/')
write.table(fam,'clean_validation.fam',row.names=F,col.names=F,quote=F)
write.table(fam0,'clean_validation_HDGC0.fam',row.names=F,col.names=F,quote=F)
write.table(fam1,'clean_validation_HDGC1.fam',row.names=F,col.names=F,quote=F)

system("cp clean_validation.bed clean_validation_HDGC0.bed")
system("cp clean_validation.bed clean_validation_HDGC1.bed")
system("cp clean_validation.bim clean_validation_HDGC0.bim")
system("cp clean_validation.bim clean_validation_HDGC1.bim")


############ FARVAT - covariate
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/
  ### discovery 
  ## transcriptome
  # Model 2 (HDGC : all, MUC1 included)
  farvat --bed 0.data/2.plink/whole/clean_discovery.bed --set 0.data/discovery_tc.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/1.discovery/1.tc/2.model2/res_cov --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC0.bed --set 0.data/discovery_tc.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/1.discovery/1.tc/5.model5/res_cov --ignoreparent --skato --thread 1 --1case5
# Model 6 (HDGC : 1, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC1.bed --set 0.data/discovery_tc.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/1.discovery/1.tc/6.model6/res_cov --ignoreparent --skato --thread 15 --1case

## gene
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery.bed --set 0.data/discovery_gene.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/1.discovery/2.gene/2.model2/res_cov --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC0.bed --set 0.data/discovery_gene.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/1.discovery/2.gene/5.model5/res_cov --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC1.bed --set 0.data/discovery_gene.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/1.discovery/2.gene/6.model6/res_cov --ignoreparent --skato --thread 15 --1case


### validation 
## transcriptome
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation.bed --set 0.data/validation_tc.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/2.validation/1.tc/2.model2/res_cov --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC0.bed --set 0.data/validation_tc.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/2.validation/1.tc/5.model5/res_cov --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included) 
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC1.bed --set 0.data/validation_tc.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/2.validation/1.tc/6.model6/res_cov --ignoreparent --skato --thread 15 --1case

## gene
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation.bed --set 0.data/validation_gene.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/2.validation/2.gene/2.model2/res_cov --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC0.bed --set 0.data/validation_gene.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/2.validation/2.gene/5.model5/res_cov --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included) 
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC1.bed --set 0.data/validation_gene.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/2.validation/2.gene/6.model6/res_cov --ignoreparent --skato --thread 15 --1case


############ FARVAT - prevalence
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/
  ### discovery 
  ## transcriptome
  # Model 2 (HDGC : all, MUC1 included)
  farvat --bed 0.data/2.plink/whole/clean_discovery.bed --set 0.data/discovery_tc.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/1.discovery/1.tc/2.model2/res_prev --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC0.bed --set 0.data/discovery_tc.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/1.discovery/1.tc/5.model5/res_prev --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC1.bed --set 0.data/discovery_tc.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/1.discovery/1.tc/6.model6/res_prev --ignoreparent --skato --thread 15 --1case

## gene
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery.bed --set 0.data/discovery_gene.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/1.discovery/2.gene/2.model2/res_prev --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC0.bed --set 0.data/discovery_gene.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/1.discovery/2.gene/5.model5/res_prev --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC1.bed --set 0.data/discovery_gene.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/1.discovery/2.gene/6.model6/res_prev --ignoreparent --skato --thread 15 --1case


### validation 
## transcriptome
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation.bed --set 0.data/validation_tc.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/2.validation/1.tc/2.model2/res_prev --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC0.bed --set 0.data/validation_tc.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/2.validation/1.tc/5.model5/res_prev --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included) 
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC1.bed --set 0.data/validation_tc.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/2.validation/1.tc/6.model6/res_prev --ignoreparent --skato --thread 15 --1case

## gene
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation.bed --set 0.data/validation_gene.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/2.validation/2.gene/2.model2/res_prev --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC0.bed --set 0.data/validation_gene.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/2.validation/2.gene/5.model5/res_prev --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included) 
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC1.bed --set 0.data/validation_gene.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/2.validation/2.gene/6.model6/res_prev --ignoreparent --skato --thread 15 --1case


### Summary
## get plotdata
Mode <- c('tc','gene')
get.plotdata <- function(ii,jj){
  # covariate
  aa <- read.table(paste0(ii,'.',Mode[ii],'/',jj,'.model',jj,'/res_cov.gene.res'),head=T,stringsAsFactor=F)
  plotdata <- na.omit(aa[,c(1,2,ncol(aa),3,4,5,7,9,10)])
  colnames(plotdata)[1:3] <- c('CHR','GENE','p.val')
  plotdata <- plotdata[with(plotdata,order(CHR,START)),]
  write.table(plotdata,paste0(ii,".",Mode[ii],"/SKATO_cov_",Mode[ii],"_model",jj,"_plotdata.txt"),row.names=F,quote=F)
  
  # prevalence
  aa <- read.table(paste0(ii,'.',Mode[ii],'/',jj,'.model',jj,'/res_prev.gene.res'),head=T,stringsAsFactor=F)
  plotdata <- na.omit(aa[,c(1,2,ncol(aa),3,4,5,7,9,10)])
  colnames(plotdata)[1:3] <- c('CHR','GENE','p.val')
  plotdata <- plotdata[with(plotdata,order(CHR,START)),]
  write.table(plotdata,paste0(ii,".",Mode[ii],"/SKATO_prev_",Mode[ii],"_model",jj,"_plotdata.txt"),row.names=F,quote=F)
}

# discovery
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/1.discovery')
for(ii in 1:2) for(jj in c(2,5,6)) get.plotdata(ii,jj)

# validation
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/2.validation')
for(ii in 1:2) for(jj in c(2,5,6)) get.plotdata(ii,jj)

## results
Mode <- c('tc','gene')
get.results <- function(ii,jj){
  for(type in c('cov','prev')){
    plotdata <- read.table(paste0(ii,".",Mode[ii],"/SKATO_",type,"_",Mode[ii],"_model",jj,"_plotdata.txt"),head=T,stringsAsFactor=F)
    
    sig.level <- 0.05/nrow(plotdata)
    
    ## significant results
    sig.res <- plotdata[plotdata$p.val<=sig.level,]
    
    ## QQ plot
    png(paste(ii,".",Mode[ii],"/SKATO_",type,"_",Mode[ii],"_model",jj,"_QQ.png",sep=""), width=900, height=1000, units = "px", bg = "white", res = 200)
    
    logit <- plotdata
    p.val <- logit[,'p.val']
    y <- -log(p.val,10)
    v <- -log10(0.05/sum(is.finite(y)))
    o.y <- sort(y[is.finite(y)],decreasing=T)
    xx<- (1:length(o.y))/length(o.y)
    x <- -log(xx,10)
    maxy <- max(o.y)
    miny <- min(o.y)
    plot(x, x, xlim=c(miny,maxy), ylim=c(miny,maxy), type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")))
    N=length(o.y)
    c95 <- rep(0,N)
    c05 <- rep(0,N)
    for(i in 1:N){
      c95[i] <- qbeta(0.95,i,N-i+1)
      c05[i] <- qbeta(0.05,i,N-i+1)
    }
    abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
    polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
    abline(a=0,b=1,col="black", lty=1)
    points(x, o.y, cex = .5, col = "dark red")
    
    dev.off()
    
    
    ## MHT plot
    png(paste(ii,".",Mode[ii],"/SKATO_",type,"_",Mode[ii],"_model",jj,"_MHT.png",sep=""), width=2000, height=1000, units = "px", bg = "white", res = 200)
    
    in_data <- plotdata[,c('CHR','GENE','p.val')]
    
    kare <- table(in_data$CHR)
    CM <- cumsum(kare)
    n.markers <- sum(kare)
    n.chr <- length(kare)
    
    pos.1 <- 1
    pos <- c()
    for(i in 1:length(kare)){
      diff <- max(kare) - kare[i]
      d <- ifelse(diff%%2==0,diff/2,(diff-1)/2)
      pos <- c(pos,(pos.1+d):(pos.1+d+kare[i]-1))
      pos.1 <- ifelse(diff%%2==0,pos[length(pos)]+d+1,pos[length(pos)]+d+2)
    }
    
    y <- -log(in_data$p.val,10)
    max.y <- ifelse(max(y)<6,6,max(y))
    
    colors <- rainbow(n.chr)
    par(xaxt = "n", yaxt = "n")
    plot(pos, y, type = "n", xlab = "Chromosome", ylab = expression(paste(-log[10],"(p-value)")), axes = FALSE, ylim=c(0,max.y+1),cex.main=1.5, cex.lab=1.2)
    axis(1, tick = FALSE)
    axis(2, tick = FALSE)
    for (i in 1:n.chr) {
      u <- CM[i]
      l <- CM[i] - kare[i] + 1
      cat("Plotting points ", l, "-", u, "\n")
      chr <- l:u
      tran.y <- y[chr]
      points(pos[chr], tran.y, col = colors[i],pch=20,cex=0.5)
    }
    par(xaxt = "s", yaxt = "s")
    AT <- cumsum(c(0,rep(max(kare),length(kare)-1)))+max(kare)/2
    LB <- 1:n.chr
    LB[LB==23] <- 'X'
    LB[LB==24] <- 'Y'
    axis(1, at = AT, labels = LB)
    if(max.y+1 < 10) {
      axis(2, at = 0:(max.y+1))
    } else {
      axis(2, at = seq(0,max.y+1,by=2))
    }
    abline(h=-log10(sig.level),lty=2,col='gray')
    dev.off()
  }
}

### discovery
#setwd('C:\\Users\\김원지\\Desktop\\GC\\farvat\\whole_gene\\1.discovery')
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/1.discovery')
for(ii in 1:2) for(jj in c(2,5,6)) get.results(ii,jj)

### validataion
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/2.validation')
for(ii in 1:2) for(jj in c(2,5,6)) get.results(ii,jj)


###### MUC4 ######
#### setfile
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/0.data')
muc4 <- read.table('MUC4_SNPset_ver1.txt',head=F,stringsAsFactor=F)
muc4 <- muc4[muc4[,1]=='MCU4',]
muc4[,1] <- 'MUC4'
muc4[,2] <- gsub('chr3','3:',muc4[,2])
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data')
write.table(muc4,'MUC4.setfile',col.names=F,row.names=F,quote=F)

#### farvat
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/
  
  ### prevalence
  ## discovery 
  # Model 2 (HDGC : all, MUC1 included)
  farvat --bed 0.data/2.plink/whole/clean_discovery.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/1.discovery/3.MUC4/model2_MUC4_prev --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC0.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/1.discovery/3.MUC4/model5_MUC4_prev --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC1.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/1.discovery/3.MUC4/model6_MUC4_prev --ignoreparent --skato --thread 15 --1case

## validation 
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/2.validation/3.MUC4/model2_MUC4_prev --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC0.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/2.validation/3.MUC4/model5_MUC4_prev --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included) 
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC1.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/2.validation/3.MUC4/model6_MUC4_prev --ignoreparent --skato --thread 15 --1case


### covariate
## discovery 
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/1.discovery/3.MUC4/model2_MUC4_cov --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC0.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/1.discovery/3.MUC4/model5_MUC4_cov --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_discovery_HDGC1.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/discovery_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/1.discovery/3.MUC4/model6_MUC4_cov --ignoreparent --skato --thread 15 --1case

## validation 
# Model 2 (HDGC : all, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,HDGC,MUC1 --mispheno NA --out 1.output/2.validation/3.MUC4/model2_MUC4_cov --ignoreparent --skato --thread 15 --1case
# Model 5 (HDGC : 0, MUC1 included)
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC0.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/2.validation/3.MUC4/model5_MUC4_cov --ignoreparent --skato --thread 15 --1case
# Model 6 (HDGC : 1, MUC1 included) 
farvat --bed 0.data/2.plink/whole/clean_validation_HDGC1.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/validation_KM.txt --sampvar ../../0.data/new_phenofile.txt --cname Sex,Smoking01,MUC1 --mispheno NA --out 1.output/2.validation/3.MUC4/model6_MUC4_cov --ignoreparent --skato --thread 15 --1case


## results
Mode <- c('tc','gene')
get.results <- function(ii,jj){
  for(type in c('cov','prev')){
    plotdata <- read.table(paste0(ii,".",Mode[ii],"/SKATO_",type,"_",Mode[ii],"_model",jj,"_plotdata.txt"),head=T,stringsAsFactor=F)
    muc4 <- read.table(paste0('3.MUC4/model',jj,'_MUC4_',type,'.gene.res'),head=T,stringsAsFactor=F)
    muc4 <- muc4[,c(1,2,ncol(muc4)-4,3,4,5,7,9,10)]
    colnames(muc4) <- colnames(plotdata)
    plotdata <- rbind(plotdata,muc4)
    plotdata <- plotdata[with(plotdata,order(CHR,START)),]
    
    sig.level <- 0.05/nrow(plotdata)
    
    ## significant results
    sig.res <- plotdata[plotdata$p.val<=7e-5,]
    write.csv(sig.res,paste("3.MUC4/SKATO_",type,"_",Mode[ii],"_model",jj,".csv",sep=""),row.names=F,quote=F)
    
    ## QQ plot
    png(paste("3.MUC4/SKATO_",type,"_",Mode[ii],"_model",jj,"_QQ.png",sep=""), width=900, height=1000, units = "px", bg = "white", res = 200)
    
    logit <- plotdata
    p.val <- logit[,'p.val']
    y <- -log(p.val,10)
    v <- -log10(0.05/sum(is.finite(y)))
    o.y <- sort(y[is.finite(y)],decreasing=T)
    xx<- (1:length(o.y))/length(o.y)
    x <- -log(xx,10)
    maxy <- max(o.y)
    miny <- min(o.y)
    plot(x, x, xlim=c(miny,maxy), ylim=c(miny,maxy), type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")))
    N=length(o.y)
    c95 <- rep(0,N)
    c05 <- rep(0,N)
    for(i in 1:N){
      c95[i] <- qbeta(0.95,i,N-i+1)
      c05[i] <- qbeta(0.05,i,N-i+1)
    }
    abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
    polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
    abline(a=0,b=1,col="black", lty=1)
    points(x, o.y, cex = .5, col = "dark red")
    
    dev.off()
    
    
    ## MHT plot
    png(paste("3.MUC4/SKATO_",type,"_",Mode[ii],"_model",jj,"_MHT.png",sep=""), width=2000, height=1000, units = "px", bg = "white", res = 200)
    
    in_data <- plotdata[,c('CHR','GENE','p.val')]
    
    kare <- table(in_data$CHR)
    CM <- cumsum(kare)
    n.markers <- sum(kare)
    n.chr <- length(kare)
    
    pos.1 <- 1
    pos <- c()
    for(i in 1:length(kare)){
      diff <- max(kare) - kare[i]
      d <- ifelse(diff%%2==0,diff/2,(diff-1)/2)
      pos <- c(pos,(pos.1+d):(pos.1+d+kare[i]-1))
      pos.1 <- ifelse(diff%%2==0,pos[length(pos)]+d+1,pos[length(pos)]+d+2)
    }
    
    y <- -log(in_data$p.val,10)
    max.y <- ifelse(max(y)<6,6,max(y))
    
    colors <- rainbow(n.chr)
    par(xaxt = "n", yaxt = "n")
    plot(pos, y, type = "n", xlab = "Chromosome", ylab = expression(paste(-log[10],"(p-value)")), axes = FALSE, ylim=c(0,max.y+1),cex.main=1.5, cex.lab=1.2)
    axis(1, tick = FALSE)
    axis(2, tick = FALSE)
    for (i in 1:n.chr) {
      u <- CM[i]
      l <- CM[i] - kare[i] + 1
      cat("Plotting points ", l, "-", u, "\n")
      chr <- l:u
      tran.y <- y[chr]
      points(pos[chr], tran.y, col = colors[i],pch=20,cex=0.5)
    }
    par(xaxt = "s", yaxt = "s")
    AT <- cumsum(c(0,rep(max(kare),length(kare)-1)))+max(kare)/2
    LB <- 1:n.chr
    LB[LB==23] <- 'X'
    LB[LB==24] <- 'Y'
    axis(1, at = AT, labels = LB)
    if(max.y+1 < 10) {
      axis(2, at = 0:(max.y+1))
    } else {
      axis(2, at = seq(0,max.y+1,by=2))
    }
    abline(h=-log10(sig.level),lty=2,col='gray')
    dev.off()
  }
}


### discovery
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/1.discovery')
for(ii in 1:2) for(jj in c(2,5,6)) get.results(ii,jj)

### validataion
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/2.validation')
for(ii in 1:2) for(jj in c(2,5,6)) get.results(ii,jj)

### MUC4
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output')
dis <- read.table('1.discovery/3.MUC4/model2_MUC4_prev.gene.res',head=T,stringsAsFactor=F)
val <- read.table('2.validation/3.MUC4/model2_MUC4_prev.gene.res',head=T,stringsAsFactor=F)
muc4 <- rbind(dis,val)
muc4 <- muc4[,c(1,2,ncol(muc4)-4,ncol(muc4),3,4,5,7,9,10)]
rownames(muc4) <- c('discovery','validation')
write.csv(muc4,'MUC4_results.csv',quote=F)


######## only Prevalence #############

## results
Mode <- c('tc','gene')
get.results <- function(ii,jj){
  for(type in c('prev')){
    plotdata <- read.table(paste0(ii,".",Mode[ii],"\\SKATO_",type,"_",Mode[ii],"_model",jj,"_plotdata.txt"),head=T,stringsAsFactor=F)
    muc4 <- read.table(paste0('3.MUC4\\model',jj,'_MUC4_',type,'.gene.res'),head=T,stringsAsFactor=F)
    muc4 <- muc4[,c(1,2,ncol(muc4),3,4,5,7,9,10)]
    colnames(muc4) <- colnames(plotdata)
    plotdata <- rbind(plotdata,muc4)
    plotdata <- plotdata[with(plotdata,order(CHR,START)),]
    
    sig.level <- 0.05/nrow(plotdata)
    
    ## significant results
    sig.res <- plotdata[plotdata$p.val<=7e-5,]
    write.csv(sig.res,paste("3.MUC4\\SKATO_",type,"_",Mode[ii],"_model",jj,".csv",sep=""),row.names=F,quote=F)
    
    ## QQ plot
    png(paste("3.MUC4\\SKATO_",type,"_",Mode[ii],"_model",jj,"_QQ.png",sep=""), width=900, height=1000, units = "px", bg = "white", res = 200)
    
    logit <- plotdata
    p.val <- logit[,'p.val']
    y <- -log(p.val,10)
    v <- -log10(0.05/sum(is.finite(y)))
    o.y <- sort(y[is.finite(y)],decreasing=T)
    xx<- (1:length(o.y))/length(o.y)
    x <- -log(xx,10)
    maxy <- max(o.y)
    miny <- min(o.y)
    plot(x, x, xlim=c(miny,maxy), ylim=c(miny,maxy), type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")))
    N=length(o.y)
    c95 <- rep(0,N)
    c05 <- rep(0,N)
    for(i in 1:N){
      c95[i] <- qbeta(0.95,i,N-i+1)
      c05[i] <- qbeta(0.05,i,N-i+1)
    }
    abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
    polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
    abline(a=0,b=1,col="black", lty=1)
    points(x, o.y, cex = .5, col = "dark red")
    
    dev.off()
    
    
    ## MHT plot
    png(paste("3.MUC4\\SKATO_",type,"_",Mode[ii],"_model",jj,"_MHT.png",sep=""), width=2000, height=1000, units = "px", bg = "white", res = 200)
    
    in_data <- plotdata[,c('CHR','GENE','p.val')]
    
    kare <- table(in_data$CHR)
    CM <- cumsum(kare)
    n.markers <- sum(kare)
    n.chr <- length(kare)
    
    pos.1 <- 1
    pos <- c()
    for(i in 1:length(kare)){
      diff <- max(kare) - kare[i]
      d <- ifelse(diff%%2==0,diff/2,(diff-1)/2)
      pos <- c(pos,(pos.1+d):(pos.1+d+kare[i]-1))
      pos.1 <- ifelse(diff%%2==0,pos[length(pos)]+d+1,pos[length(pos)]+d+2)
    }
    
    y <- -log(in_data$p.val,10)
    max.y <- ifelse(max(y)<6,6,max(y))
    
    colors <- rainbow(n.chr)
    par(xaxt = "n", yaxt = "n")
    plot(pos, y, type = "n", xlab = "Chromosome", ylab = expression(paste(-log[10],"(p-value)")), axes = FALSE, ylim=c(0,max.y+1),cex.main=1.5, cex.lab=1.2)
    axis(1, tick = FALSE)
    axis(2, tick = FALSE)
    for (i in 1:n.chr) {
      u <- CM[i]
      l <- CM[i] - kare[i] + 1
      cat("Plotting points ", l, "-", u, "\n")
      chr <- l:u
      tran.y <- y[chr]
      points(pos[chr], tran.y, col = colors[i],pch=20,cex=0.5)
    }
    par(xaxt = "s", yaxt = "s")
    AT <- cumsum(c(0,rep(max(kare),length(kare)-1)))+max(kare)/2
    LB <- 1:n.chr
    LB[LB==23] <- 'X'
    LB[LB==24] <- 'Y'
    axis(1, at = AT, labels = LB)
    if(max.y+1 < 10) {
      axis(2, at = 0:(max.y+1))
    } else {
      axis(2, at = seq(0,max.y+1,by=2))
    }
    abline(h=-log10(sig.level),lty=2,col='gray')
    dev.off()
  }
}


### discovery
setwd('C:\\Users\\김원지\\Downloads\\whole_gene\\1.discovery')
for(ii in 1:2) get.results(ii,2)

### validataion
setwd('C:\\Users\\김원지\\Downloads\\whole_gene\\2.validation')
for(ii in 1:2) get.results(ii,2)

### MUC4
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output')
dis <- read.table('1.discovery/3.MUC4/model2_MUC4_prev.gene.res',head=T,stringsAsFactor=F)
val <- read.table('2.validation/3.MUC4/model2_MUC4_prev.gene.res',head=T,stringsAsFactor=F)
muc4 <- rbind(dis,val)
muc4 <- muc4[,c(1,2,ncol(muc4)-4,ncol(muc4),3,4,5,7,9,10)]
rownames(muc4) <- c('discovery','validation')
write.csv(muc4,'MUC4_results.csv',quote=F)


################# discovery + validation #############################3
#### merge separated vcf files for each family into one file
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink')
dis_fams <- c('family10','sub_family11','sub_family1','sub_family2','family3','family5','family6','family7','family8','sub_family9')
val_fams <- c('family12_new','family13_new','family14','family15')
all_fams <- c(dis_fams,val_fams)

## merge list
merlist.dis <- paste0('whole/whole_final_biSNP_',all_fams[-1])
write.table(merlist.dis,'whole/whole_allmerlist_discovery.txt',col.names=F,row.names=F,quote=F)

## merge
system(paste0('plink --bfile whole/whole_final_biSNP_',all_fams[1],' --merge-list whole/whole_allmerlist_discovery.txt --make-bed --out whole/allmerge/allmerge'))
for(i in all_fams){
  system(paste0('plink --bfile whole/whole_final_biSNP_',i,' --exclude whole/allmerge/allmerge-merge.missnp --make-bed --out whole/allmerge/biSNP_',i))
}
merlist.dis <- paste0('biSNP_',all_fams[-1])
write.table(merlist.dis,'whole/allmerge/allmerge_realbiSNP.txt',col.names=F,row.names=F,quote=F)
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole/allmerge')
system(paste0('plink --bfile biSNP_',all_fams[1],' --merge-list allmerge_realbiSNP.txt --make-bed --out allmerge'))
# Finally, 4138346 SNPs remained.

## maf0 SNPs
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole/allmerge')

bim.dis <- read.table("allmerge.bim",head=F,stringsAsFactor=F)
rm.dis <- bim.dis[bim.dis[,5]==0,]	# 248367 SNPs
write.table(rm.dis,'maf0.txt',row.names=F,col.names=F,quote=F)

### QC
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole/allmerge
plink --bfile allmerge --geno 0.05 --hwe 1e-5 --exclude maf0.txt --make-bed --out qc_snp
# 4138346(Total)-2593882(geno<0.05)-13606(hwe<1e-5)-248367(0 maf)=1282491
plink --bfile qc_snp --mind 0.05 --make-bed --out clean_allmerge
# No subjects were removed

### Kinship coefficient
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole/')
fam.dis <- read.table('clean_discovery.fam',head=F,stringsAsFactor=F)
fam.val <- read.table('clean_validation.fam',head=F,stringsAsFactor=F)

KM.dis <- as.matrix(read.table('discovery_KM.txt',head=F,stringsAsFactor=F))
KM.val <- as.matrix(read.table('validation_KM.txt',head=F,stringsAsFactor=F))

fam <- read.table('allmerge/clean_allmerge.fam',head=F,stringsAsFactor=F)
fam[,2] <- gsub('-0','-',fam[,2])
new.KM <- matrix(0,nrow(fam),nrow(fam))
new.KM[match(fam.dis[,2],fam[,2]),match(fam.dis[,2],fam[,2])] <- KM.dis
new.KM[match(fam.val[,2],fam[,2]),match(fam.val[,2],fam[,2])] <- KM.val

write.table(new.KM,'allmerge/allmerge_KM.txt',col.names=F,row.names=F,quote=F)


### Set file
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole/allmerge')
nm <- readLines("/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/NM_position.txt")
get.geneinfo <- function(i,bim){
  if(i%%100==0)print(i)
  genes <- strsplit(nm[2*i-1],">|\t")[[1]]
  tc <- strsplit(nm[2*i],"\t")[[1]]
  chr <- ifelse(tc[1]=='chrX',23,ifelse(tc[1]=='chrY',24,gsub('chr','',tc[1])))
  strand <- tc[2]
  
  TC <- do.call(rbind,strsplit(tc[3:length(tc)],';'))
  command_str <- paste0('SNPs <- bim[bim$V1==',chr,'& (',paste(apply(TC[,1:2,drop=F],1,function(xx) paste0('(bim$V4 >=', paste(xx,collapse=' & bim$V4 <='),')')),collapse='|'),'),2]')
  eval(parse(text=command_str))
  snps <- ifelse(length(SNPs)==0,'',paste0(SNPs,collapse=';'))
  
  return(data.frame(NM_ID=genes[2],Gene_symbol=genes[3],CHR=chr,SNPs=snps,stringsAsFactors=F))
}

library(parallel)
bim <- read.table('clean_allmerge.bim',head=F,stringsAsFactor=F)
res <- mclapply(1:(length(nm)/2),get.geneinfo,bim=bim,mc.cores=20)
res.dis <- do.call(rbind,res)
res.dis.indata <- res.dis[grep(';',res.dis$SNPs),]	# 11105 genes & 18293 transcriptome
write.table(res.dis,'gene_SNP_matching.txt',row.names=F,quote=F)

## transcriptome
get.setfile.tc <- function(i,indata){
  tc <- indata[i,1]
  snps <- strsplit(indata[i,4],';')[[1]]
  return(data.frame(gene=tc,SNP=snps,stringsAsFactors=F))
}

# discovery
sf.dis <- mclapply(1:nrow(res.dis.indata),get.setfile.tc,indata=res.dis.indata,mc.cores=15)
setfile.dis <- do.call(rbind,sf.dis)
setfile.dis <- setfile.dis[order(setfile.dis[,1]),]
write.table(setfile.dis,'allmerge_tc.setfile',row.names=F,col.names=F,quote=F)


## gene
get.setfile.gene <- function(gene,indata){
  tc <- indata[indata[,2]==gene,]
  snps <- unique(do.call(c,strsplit(tc[,4],';')))
  return(data.frame(gene=gene,SNP=snps,stringsAsFactors=F))
}

# discovery
sf.dis <- mclapply(unique(res.dis.indata[,2]),get.setfile.gene,indata=res.dis.indata,mc.cores=15)
setfile.dis <- do.call(rbind,sf.dis)
setfile.dis <- setfile.dis[order(setfile.dis[,1]),]
write.table(setfile.dis,'allmerge_gene.setfile',row.names=F,col.names=F,quote=F)

#### HDGC group, family id
setwd("/home2/wjkim/project/gastric_cancer/3.FARVAT")
pheno <- read.table('../0.data/new_phenofile.txt',head=T)

# discovery
fam <- read.table('2.whole_gene/0.data/2.plink/whole/allmerge/clean_allmerge.fam',head=F,stringsAsFactor=F)
fam[,1] <- pheno[match(fam$V2,pheno$IID),'FID']
fam[,6] <- pheno[match(fam$V2,pheno$IID),'GC']

HDGC0 <- pheno$IID[which(pheno$HDGC==0)]
HDGC1 <- pheno$IID[which(pheno$HDGC==1)]

fam0 <- fam1 <- fam
fam0$V6[fam0$V2%in%HDGC1] <- NA
fam1$V6[fam1$V2%in%HDGC0] <- NA

setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole/allmerge')
write.table(fam,'clean_allmerge.fam',row.names=F,col.names=F,quote=F)
write.table(fam0,'clean_allmerge_HDGC0.fam',row.names=F,col.names=F,quote=F)
write.table(fam1,'clean_allmerge_HDGC1.fam',row.names=F,col.names=F,quote=F)

system("cp clean_allmerge.bed clean_allmerge_HDGC0.bed")
system("cp clean_allmerge.bed clean_allmerge_HDGC1.bed")
system("cp clean_allmerge.bim clean_allmerge_HDGC0.bim")
system("cp clean_allmerge.bim clean_allmerge_HDGC1.bim")


############ FARVAT - prevalence
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/0.data/2.plink/whole/allmerge/
  ## transcriptome
  # Model 2 (HDGC : all, MUC1 included)
  farvat --bed clean_allmerge.bed --set allmerge_tc.setfile --genetest --cor allmerge_KM.txt --prevalence 0.0046 --mispheno NA --out /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/3.allmerge/1.tc/2.model2/res_prev --ignoreparent --skato --thread 20 --1case

## gene
# Model 2 (HDGC : all, MUC1 included)
farvat --bed clean_allmerge.bed --set allmerge_gene.setfile --genetest --cor allmerge_KM.txt --prevalence 0.0046 --mispheno NA --out /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/3.allmerge/2.gene/2.model2/res_prev --ignoreparent --skato --thread 20 --1case

## Summary
Mode <- c('tc','gene')
get.plotdata <- function(ii,jj){
  # prevalence
  aa <- read.table(paste0(ii,'.',Mode[ii],'/',jj,'.model',jj,'/res_prev.gene.res'),head=T,stringsAsFactor=F)
  plotdata <- na.omit(aa[,c(1,2,ncol(aa),3,4,5,7,9,10)])
  colnames(plotdata)[1:3] <- c('CHR','GENE','p.val')
  plotdata <- plotdata[with(plotdata,order(CHR,START)),]
  write.table(plotdata,paste0(ii,".",Mode[ii],"/SKATO_prev_",Mode[ii],"_model",jj,"_plotdata.txt"),row.names=F,quote=F)
}

setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/3.allmerge')
for(ii in 1:2) get.plotdata(ii,2)

###### MUC4 ######
#### farvat
cd /home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/
  farvat --bed 0.data/2.plink/whole/allmerge/clean_allmerge.bed --set 0.data/MUC4.setfile --genetest --cor 0.data/2.plink/whole/allmerge/allmerge_KM.txt --prevalence 0.0046 --mispheno NA --out 1.output/3.allmerge/3.MUC4/model2_MUC4_prev --ignoreparent --skato --thread 20 --1case

###### Plots ######
Mode <- c('tc','gene')
get.results <- function(ii,jj){
  for(type in c('prev')){
    plotdata <- read.table(paste0(ii,".",Mode[ii],"\\SKATO_",type,"_",Mode[ii],"_model",jj,"_plotdata.txt"),head=T,stringsAsFactor=F)
    chrpos <- apply(plotdata[,c('CHR','START','END')],1,function(x) paste0(x,collapse=':'))
    plotdata <- plotdata[!duplicated(chrpos),]
    
    muc4 <- read.table(paste0('3.MUC4\\model',jj,'_MUC4_',type,'.gene.res'),head=T,stringsAsFactor=F)
    muc4 <- muc4[,c(1,2,ncol(muc4),3,4,5,7,9,10)]
    colnames(muc4) <- colnames(plotdata)
    plotdata <- rbind(plotdata,muc4)
    plotdata <- plotdata[with(plotdata,order(CHR,START)),]
    
    sig.level <- 0.05/nrow(plotdata)
    
    ## significant results
    sig.res <- plotdata[plotdata$p.val<=7e-5,]
    write.csv(sig.res,paste("3.MUC4\\SKATO_",type,"_",Mode[ii],"_model",jj,".csv",sep=""),row.names=F,quote=F)
    
    ## QQ plot
    png(paste("3.MUC4\\SKATO_",type,"_",Mode[ii],"_model",jj,"_QQ.png",sep=""), width=900, height=1000, units = "px", bg = "white", res = 200)
    
    logit <- plotdata
    p.val <- logit[,'p.val']
    y <- -log(p.val,10)
    v <- -log10(0.05/sum(is.finite(y)))
    o.y <- sort(y[is.finite(y)],decreasing=T)
    xx<- (1:length(o.y))/length(o.y)
    x <- -log(xx,10)
    maxy <- max(o.y)
    miny <- min(o.y)
    plot(x, x, xlim=c(miny,maxy), ylim=c(miny,maxy), type = "n", xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")))
    N=length(o.y)
    c95 <- rep(0,N)
    c05 <- rep(0,N)
    for(i in 1:N){
      c95[i] <- qbeta(0.95,i,N-i+1)
      c05[i] <- qbeta(0.05,i,N-i+1)
    }
    abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
    polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
    abline(a=0,b=1,col="black", lty=1)
    points(x, o.y, cex = .5, col = "dark red")
    
    dev.off()
    
    
    ## MHT plot
    png(paste("3.MUC4\\SKATO_",type,"_",Mode[ii],"_model",jj,"_MHT.png",sep=""), width=2000, height=1000, units = "px", bg = "white", res = 200)
    
    in_data <- plotdata[,c('CHR','GENE','p.val')]
    
    kare <- table(in_data$CHR)
    CM <- cumsum(kare)
    n.markers <- sum(kare)
    n.chr <- length(kare)
    
    pos.1 <- 1
    pos <- c()
    for(i in 1:length(kare)){
      diff <- max(kare) - kare[i]
      d <- ifelse(diff%%2==0,diff/2,(diff-1)/2)
      pos <- c(pos,(pos.1+d):(pos.1+d+kare[i]-1))
      pos.1 <- ifelse(diff%%2==0,pos[length(pos)]+d+1,pos[length(pos)]+d+2)
    }
    
    y <- -log(in_data$p.val,10)
    max.y <- ifelse(max(y)<6,6,max(y))
    
    colors <- rainbow(n.chr)
    par(xaxt = "n", yaxt = "n")
    plot(pos, y, type = "n", xlab = "Chromosome", ylab = expression(paste(-log[10],"(p-value)")), axes = FALSE, ylim=c(0,max.y+1),cex.main=1.5, cex.lab=1.2)
    axis(1, tick = FALSE)
    axis(2, tick = FALSE)
    for (i in 1:n.chr) {
      u <- CM[i]
      l <- CM[i] - kare[i] + 1
      cat("Plotting points ", l, "-", u, "\n")
      chr <- l:u
      tran.y <- y[chr]
      points(pos[chr], tran.y, col = colors[i],pch=20,cex=0.5)
    }
    par(xaxt = "s", yaxt = "s")
    AT <- cumsum(c(0,rep(max(kare),length(kare)-1)))+max(kare)/2
    LB <- 1:n.chr
    LB[LB==23] <- 'X'
    LB[LB==24] <- 'Y'
    axis(1, at = AT, labels = LB)
    if(max.y+1 < 10) {
      axis(2, at = 0:(max.y+1))
    } else {
      axis(2, at = seq(0,max.y+1,by=2))
    }
    abline(h=-log10(sig.level),lty=2,col='gray')
    dev.off()
    return(paste0(Mode[ii],':',nrow(plotdata)))
  }
}


### discovery
setwd('C:\\Users\\김원지\\Downloads\\whole_gene\\3.allmerge')
for(ii in 1:2) get.results(ii,2)


#################### 2/11/2019 : QQ plots for different number of variants in gene or transcriptome ##########################
setwd('/home2/wjkim/project/gastric_cancer/3.FARVAT/2.whole_gene/1.output/3.allmerge')
Mode <- c('tc','gene')
get.results <- function(ii,jj){
  for(type in c('prev')){
    plotdata <- read.table(paste0(ii,".",Mode[ii],"/SKATO_",type,"_",Mode[ii],"_model",jj,"_plotdata.txt"),head=T,stringsAsFactor=F)
    chrpos <- apply(plotdata[,c('CHR','START','END')],1,function(x) paste0(x,collapse=':'))
    plotdata <- plotdata[!duplicated(chrpos),]
    plotdata$class <- ifelse(plotdata$NVARIANT<10,1,ifelse(plotdata$NVARIANT<30,2,ifelse(plotdata$NVARIANT<50,3,4)))
    
    muc4 <- read.table(paste0('3.MUC4/model',jj,'_MUC4_',type,'.gene.res'),head=T,stringsAsFactor=F)
    muc4 <- muc4[,c(1,2,ncol(muc4),3,4,5,7,9,10)]
    muc4$class <- ifelse(muc4$NVARIANT<10,1,ifelse(muc4$NVARIANT<30,2,ifelse(muc4$NVARIANT<50,3,4)))
    colnames(muc4) <- colnames(plotdata)
    plotdata <- rbind(plotdata,muc4)
    plotdata <- plotdata[with(plotdata,order(CHR,START)),]
    
    sig.level <- 0.05/nrow(plotdata)
    mains <- c('# variants < 10','10 <= # variants <30 ','30 <= # variants <50','# variants >= 50')
    ## QQ plot
    png(paste("3.MUC4/SKATO_",type,"_",Mode[ii],"_model",jj,"_QQ_all.png",sep=""), width=1800, height=1900, units = "px", bg = "white", res = 200)
    par(mfrow=c(2,2))
    for(k in 1:4){
      logit <- plotdata[plotdata$class==k,]
      p.val <- logit[,'p.val']
      y <- -log(p.val,10)
      v <- -log10(0.05/sum(is.finite(y)))
      o.y <- sort(y[is.finite(y)],decreasing=T)
      xx<- (1:length(o.y))/length(o.y)
      x <- -log(xx,10)
      maxy <- max(o.y)
      miny <- min(o.y)
      plot(x, x, xlim=c(miny,maxy), ylim=c(miny,maxy), type = "n", main=mains[k], xlab = expression(paste("Expected ",-log[10],"(p-value)")), ylab = expression(paste("Observed ",-log[10],"(p-value)")))
      N=length(o.y)
      c95 <- rep(0,N)
      c05 <- rep(0,N)
      for(i in 1:N){
        c95[i] <- qbeta(0.95,i,N-i+1)
        c05[i] <- qbeta(0.05,i,N-i+1)
      }
      abline(h = c(0:max(max(x),max(o.y[is.finite(o.y)]))), v =c(0:max(max(x),max(o.y[is.finite(o.y)]))), col = "darkgray", lty=3)
      polygon(c(x,sort(x)),c(-log(c95,10),sort(-log(c05,10))),col=c("gray"),border=NA)
      abline(a=0,b=1,col="black", lty=1)
      points(x, o.y, cex = .5, col = 'dark red')
    }
    
    dev.off()
    return(paste0(Mode[ii],':',nrow(plotdata)))
  }
}

for(ii in 1:2) get.results(ii,2)



