library("poolfstat")
library("geosphere")
library("vegan")
library("pcadapt")
library("car")
library("DescTools")
library("lfmm")
library("reshape2")

##Read population names and coordinates
pool.names<-scan("~/Dropbox/Spruce/Spruce_reseq/Helenas_project/PoolSeq/new_analyses_aug22/samples.txt",what="character")
coords<-read.table("~/Dropbox/Spruce/Spruce_reseq/Helenas_project/PoolSeq/new_analyses_aug22/coords.csv",head=T,sep=",")
##Population type
idx<-(coords$Type=="Old")
cols<-c("tan2","darkgreen")
col.vector<-sapply(idx,function(x){return(cols[x+1])})


##Select climate variables
clim.var<-c("annualPET","aridityIndexThornthwaite","climaticMoistureIndex","continentality","embergerQ","growingDegDays0","growingDegDays5","maxTempColdest","minTempWarmest","monthCountByTemp10","PETColdestQuarter","PETDriestQuarter","PETseasonality","PETWarmestQuarter","PETWettestQuarter","thermicityIndex")


##Read SNP data
pools.vscn<-vcf2pooldata(vcf.file="/Users/pron0005/Dropbox/Spruce/Spruce_reseq/Helenas_project/PoolSeq/new_analyses_aug22/All_Pops_sorted_full_varscan_no_chloroplast.vcf",
                         poolsizes=rep(60,45),poolnames=pool.names,min.cov.per.pool = -1,remove.indels = T,
                         min.dist.from.indels=5)


##Filter SNPs
##All populations
pools.varscan<-pooldata.subset(pools.vscn,min.cov.per.pool = 60 ,cov.qthres.per.pool = c(0.001,0.999), min.maf = 0.001)
##Old populations
pools.varscan.old<-pooldata.subset(pools.varscan,pool.index=c(1:45)[idx])#,min.cov.per.pool = 60 ,cov.qthres.per.pool = c(0.001,0.999), min.maf = 0.001)
##Planted populations
pools.varscan.young<-pooldata.subset(pools.varscan,pool.index=c(1:45)[!idx])#,min.cov.per.pool = 60 ,cov.qthres.per.pool = c(0.001,0.999), min.maf = 0.001)


##Calculate Fst for SNPs
fst.varscan<-computeFST(pools.varscan,method="Identity")
fst.varscan.old<-computeFST(pools.varscan.old,method="Identity")
fst.varscan.young<-computeFST(pools.varscan.young,method="Identity")


##Calculate pairwise Fst among populations
fst.p.varscan<-compute.pairwiseFST(pools.varscan,method="Identity")
fst.p.varscan.old<-compute.pairwiseFST(pools.varscan.old,method="Identity")
fst.p.varscan.young<-compute.pairwiseFST(pools.varscan.young,method="Identity")


##Convert read counts to allele frequencies
pools.af<-(pools.varscan@refallele.readcount/pools.varscan@readcoverage)


##Distance matrix for populations from latitude and longitude
Dmat<-matrix(rep(NA,45*45),nrow=45,ncol=45)
for(i in 1:45){
  for(j in 1:45){
    Dmat[i,j]<-distm(c(coords[i,"Long"],coords[i,"Lat"]),c(coords[j,"Long"],coords[j,"Lat"]),fun=distHaversine)
  }
}
Dmat<-Dmat/1000


##Isolation by distance
##Old
par(mfrow=c(1,2))
par(mgp=c(3,0.7,0))
plot(fst.p.varscan.old@PairwiseFSTmatrix~Dmat[idx,idx],las=1,col="forestgreen",pch=19,xlab="Geographic distance (km)",ylab="Genetic distance (Fst)",ylim=c(0.02,0.08))
mantel(fst.p.varscan.old@PairwiseFSTmatrix,Dmat[idx,idx])
lm1<-lm(as.vector(fst.p.varscan.old@PairwiseFSTmatrix)~as.vector(Dmat[idx,idx]))
summary(lm1)
abline(lm1,lty=3)
#Planted
plot(fst.p.varscan.young@PairwiseFSTmatrix~Dmat[!idx,!idx],las=1,col="tan2",pch=19,xlab="Geographic distance (km)",ylab="Genetic distance (Fst)",ylim=c(0.02,0.08))
mantel(fst.p.varscan.young@PairwiseFSTmatrix,Dmat[!idx,!idx])
lm1<-lm(as.vector(fst.p.varscan.young@PairwiseFSTmatrix)~as.vector(Dmat[!idx,!idx]))
summary(lm1)
abline(lm1,lty=3)


##Heterozygosity
f.stats.varscan=compute.fstats(pools.varscan,computeDstat = FALSE,return.F4.blockjackknife.samples = FALSE)

par(mfrow=c(1,3))
par(mar=c(5,5,3,1))
par(mgp=c(3,0.7,0))
##Number of non-variable sites per population and type
plot(47552-apply(pools.varscan@refallele.readcount>0,2,sum)~as.factor(idx+1),outline=F,col=c("tan2","forestgreen"),xlab="Population",xaxt="n",las=1,ylab="Non-variable sites")
summary(lm(47552-apply(pools.varscan@refallele.readcount>0,2,sum)~as.factor(idx+1)))
leveneTest(lm(47552-apply(pools.varscan@refallele.readcount>0,2,sum)~as.factor(idx+1)))
axis(1,at=c(1,2),c("Planted","Old"))

##Proportion of rare allels (p<0.05)
plot(apply(pools.af<0.05,2,mean)~as.factor(idx+1),outline=F,col=c("tan2","forestgreen"),xlab="Population",xaxt="n",ylab="Proportion of rare alleles",las=1)
summary(lm(apply(pools.af<0.05,2,mean)~as.factor(idx+1)))
leveneTest(lm(apply(pools.af<0.05,2,mean)~as.factor(idx+1)))
axis(1,at=c(1,2),c("Planted","Old"))

##Within-population diversities Old vs Planted
plot(f.stats.varscan@heterozygosities$Estimate~as.factor((idx+1)),xaxt="n",yaxt="n",xlab="Population",bty="n",ylab="Heterozygosity",outline=F,col=c("tan2","forestgreen"))
summary(lm(f.stats.varscan@heterozygosities$Estimate~as.factor((idx+1))))
leveneTest(f.stats.varscan@heterozygosities$Estimate~as.factor((idx+1)))
axis(1,at=c(1,2),c("Planted","Old"))
axis(2,las=1,at=seq(0.25,0.27,by=0.005))
t.test(f.stats.varscan@heterozygosities$Estimate~as.factor(idx))
var.test(f.stats.varscan@heterozygosities$Estimate~as.factor(idx))
par(mgp=c(3,1,0))


##Compare Fst between population classes
##1 - planted vs planted, 2 - planted vs old, 4 - old vs old
#pop.pairs<-outer(as.numeric(idx+1),as.numeric(idx+1))
#plot(as.vector(fst.p.varscan@PairwiseFSTmatrix)~as.factor(pop.pairs),las=1,xlab="",ylab=expression("Pairwise F"[st]),xaxt="n",outline=FALSE,ylim=c(0,0.08),col=c("tan2","lightseagreen","forestgreen"))
#axis(1,at=c(1,2,3),labels=c("P vs P","O v P","O v O"))
#lm1<-aov(as.vector(fst.p.varscan@PairwiseFSTmatrix)~as.factor(pop.pairs))
#summary(lm1)
#TukeyHSD(lm1)
##Homogeneity of variances
#leveneTest(as.vector(fst.p.varscan@PairwiseFSTmatrix)~as.factor(pop.pairs),center="median")


par(mfrow=c(1,2))
##F2
idx.f2<-sapply(strsplit(f.stats.varscan@comparisons$F2[,1],"\\."),"[[",1)==sapply(strsplit(f.stats.varscan@comparisons$F2[,2],"\\."),"[[",1)
boxplot(f.stats.varscan@f2.values[idx.f2,]~rep(c(1,1,2),15),names=c("O vs P","P vs P"),xlab="",ylab=expression("Genetic differentiation, F"["ST"]),las=1,col=c("forestgreen","tan2"),outline=F)
kruskal.test(f.stats.varscan@f2.values[idx.f2,]~rep(c(1,1,2),15))

##Genetic similarity
idx.f3star<-sapply(strsplit(f.stats.varscan@comparisons$F3star[,1],"\\."),"[[",1)==sapply(strsplit(f.stats.varscan@comparisons$F3star[,2],"\\."),"[[",1)&sapply(strsplit(f.stats.varscan@comparisons$F3star[,1],"\\."),"[[",1)==sapply(strsplit(f.stats.varscan@comparisons$F3star[,3],"\\."),"[[",1)
idx.f3.Old<-grepl("\\.1;",dimnames(f.stats.varscan@comparisons$F3star)[[1]][idx.f3star])

boxplot(f.stats.varscan@f3star.values$Estimate[idx.f3star]~idx.f3.Old,names=c("P vs O","P vs P"),xlab="",ylab=expression("Genetic similarity, f"[3]),las=1,col=c("forestgreen","tan2"),outline=F)
kruskal.test(f.stats.varscan@f3star.values$Estimate[idx.f3star]~idx.f3.Old)



par(mfrow=c(1,2))
par(mgp=c(3,0.6,0))
##Heterozygosity versus Latitude
##Old
lm1.o<-lm(f.stats.varscan@heterozygosities$Estimate[idx]~coords$Lat[idx])
summary(lm1.o)
plot(f.stats.varscan@heterozygosities$Estimate[idx]~coords$Lat[idx],pch=19,col="forestgreen",xlab="Latitude",ylab="Heterozygosity",main="",las=1,ylim=c(0.250,0.270))
abline(lm1.o,lty=3)
##Young
lm1.y<-lm(f.stats.varscan@heterozygosities$Estimate[!idx]~coords$Lat[!idx])
summary(lm1.y)
plot(f.stats.varscan@heterozygosities$Estimate[!idx]~coords$Lat[!idx],pch=19,col="tan2",xlab="Latitude",ylab="Heterozygosity",main="",las=1,ylim=c(0.250,0.270))
abline(lm1.y,lty=3)



##Climate distance versus genetic distance (Fst)
par(mfrow=c(1,2))
par(mgp=c(3,1,0))
##Outlier SNPs
##Old
clim.idx<-c(7:22)
clim.D.old<-vegdist(samples[idx,c(clim.idx)],method="euclidean")
quantile(fst.varscan.old$snp.FST,0.995,na.rm=T)
fst.snp.idx<-(fst.varscan.old$snp.FST>0.345)
fst.snp.idx[is.na(fst.snp.idx)]<-FALSE
pools.varscan.old.outlier<-pooldata.subset(pools.varscan.old,snp.index=c(1:47552)[fst.snp.idx])
fst.p.varscan.old.outlier<-compute.pairwiseFST(pools.varscan.old.outlier,method="Identity")
fst.vec<-fst.p.varscan.old.outlier@PairwiseFSTmatrix[lower.tri(fst.p.varscan.old.outlier@PairwiseFSTmatrix)]
plot(fst.vec~sqrt(as.vector(clim.D.old)),main="",xlab="Climate distance",ylab=expression("Genetic distance, F"[ST]),pch=19,col="forestgreen",las=1,ylim=c(0.1,0.7))
lm.outlier<-lm(fst.vec~sqrt(as.vector(clim.D.old)))
summary(lm.outlier)
mantel(as.dist(fst.p.varscan.old.outlier@PairwiseFSTmatrix),clim.D.old)
abline(lm.outlier,lty=3)

##Planted
clim.D.young<-vegdist(samples[!idx,c(clim.idx)],method="euclidean",upper=T)
quantile(fst.varscan.young$snp.FST,0.995)
#fst.snp.idx<-(fst.varscan.young$snp.FST>0.357)
pools.varscan.young.outlier<-pooldata.subset(pools.varscan.young,snp.index=c(1:47552)[fst.snp.idx])
fst.p.varscan.young.outlier<-compute.pairwiseFST(pools.varscan.young.outlier,method="Identity")

fst.vec<-fst.p.varscan.young.outlier@PairwiseFSTmatrix[lower.tri(fst.p.varscan.young.outlier@PairwiseFSTmatrix)]
plot(fst.vec~sqrt(as.vector(clim.D.young)),main="",xlab="Climate distance",ylab=expression("Genetic distance, F"[ST]),pch=19,col="tan2",las=1,ylim=c(0.0,0.7))
lm.outlier<-lm(fst.vec~sqrt(as.vector(clim.D.young)))
summary(lm.outlier)
mantel(as.dist(fst.p.varscan.young.outlier@PairwiseFSTmatrix),clim.D.young)
abline(lm.outlier,lty=3)


##Null correlations for environmental variables
null.stats<-NULL
for(i in 1:100000){
  cr<-cor.test(runif(15,0,1),scale(samples[idx,sample(c(7:22),size=1)]))
  #cr<-cor.test(runif(30,0,1),samples[!idx,sample(c(7:22),size=1)])
  null.stats<-rbind(null.stats,c(i,cr$estimate,cr$p.value))
}
null.stats<-as.data.frame(null.stats)
names(null.stats)<-c("repl","corr","pval")
quantile(null.stats$corr,0.999,na.rm=T)

##SNP associations using correlations for all SNPs
##Old and Young populations separately
snp.stats.old.r<-NULL
snp.stats.young.r<-NULL
for(j in 7:22){
  temp.o<-NULL
  temp.y<-NULL
  print(j)
  for(i in 1:47552){
    cr.o<-cor.test(pools.af[i,idx],samples[idx,j])
    cr.y<-cor.test(pools.af[i,!idx],samples[!idx,j])
    temp.o<-rbind(temp.o,c(cr.o$estimate))
    temp.y<-rbind(temp.y,c(cr.y$estimate))
  }
  snp.stats.old.r<-cbind(snp.stats.old.r,temp.o)
  snp.stats.young.r<-cbind(snp.stats.young.r,temp.y)
}
snp.stats.old.r<-as.data.frame(snp.stats.old.r)
names(snp.stats.old.r)<-names(samples)[7:22]

snp.stats.young.r<-as.data.frame(snp.stats.young.r)
names(snp.stats.young.r)<-names(samples)[7:22]

##PCA for population structure
par(mfrow=c(1,2))
pca.varscan<-prcomp(t(pools.af))
plot(pca.varscan$x[,1],pca.varscan$x[,2],pch=19,col=col.vector,xlab="PC1 (10.1%)",ylab="PC2 (6.9%)",las=1,cex=1.75)
plot(pca.varscan$sdev^2,ylab="Eigenvalues",xlab="Principal components",pch=19)

##SNP associations using LFMM for all SNPs 
##Old and Young populations separately
##LFMM PCA
Y<-t(pools.af)
#Y<-t(pools.af)
pc<-prcomp(Y)
plot(pc$sdev[1:20]^2)
##Climate variable

##Climate variables to use
names(samples[,c(7:22)])[c(2,3,4,6,7,8,11)]

sig.snps<-NULL
snp.stats.old.lfmm<-NULL
snp.stats.young.lfmm<-NULL
##for(i in 6+c(2,3,4,6,7,8,11)){
for(i in 6+c(1:16)){
  Y.old<-t(pools.af)[idx,]
  Y.young<-t(pools.af)[!idx,]
  X.old<-scale(samples[idx,i],scale=T,center=T)
  X.young<-scale(samples[!idx,i],scale=T,center=T)
  
  ##Set up model
  mod.lfmm.old <- lfmm_ridge(Y = Y.old, 
                             X = X.old, 
                             K = 4)
  mod.lfmm.young <- lfmm_ridge(Y = Y.young, 
                               X = X.young, 
                               K = 4)
  
  ##ridge regression
  pv.o <- lfmm_test(Y = Y.old, 
                    X = X.old, 
                    lfmm = mod.lfmm.old, 
                    calibrate = "gif")
  pv.y <- lfmm_test(Y = Y.young, 
                    X = X.young, 
                    lfmm = mod.lfmm.young, 
                    calibrate = "gif")
  
  ##QQ-plot
  #pvalues <- pv$calibrated.pvalue 
  #qqplot(rexp(length(pvalues), rate = log(10)),
  #     -log10(pvalues), xlab = "Expected quantile",
  #     pch = 19, cex = .4)
  #abline(0,1)
  
  snp.stats.old.lfmm<-cbind(snp.stats.old.lfmm,pv.o$calibrated.pvalue)
  snp.stats.young.lfmm<-cbind(snp.stats.young.lfmm,pv.y$calibrated.pvalue)
}
##Old
snp.stats.old.lfmm<-as.data.frame(snp.stats.old.lfmm)
names(snp.stats.old.lfmm)<-names(samples)[7:22]
##Young
snp.stats.young.lfmm<-as.data.frame(snp.stats.young.lfmm)
names(snp.stats.young.lfmm)<-names(samples)[7:22]


##Partial Mantel correlations
traits<-c(8,9,10,12,13,14,17)
mstat<-NULL
for(i in traits){
  clim.idx<-c(i)

  ##Current climate
  #clim.D.old<-vegdist(samples[idx,c(clim.idx)],method="euclidean")
  #clim.D.young<-vegdist(samples[!idx,c(clim.idx)],method="euclidean",upper=T)
  
  ##RCP4.5
  clim.D.old<-vegdist(samples4.5[idx,c(clim.idx)],method="euclidean")
  clim.D.young<-vegdist(samples4.5[!idx,c(clim.idx)],method="euclidean",upper=T)
  
  ##Select SNPs using the top 0.1% (quantile 0.999) for Old populations
  #cor.snp.idx<-snp.stats.old.r[,i-6]>0.730
  #cor.snp.idx[is.na(cor.snp.idx)]<-FALSE
 
  
  ##Select SNPs using BF-corrected pvalues (1.05e-6) for Old populations
  cor.snp.idx<-snp.stats.old.lfmm[,i-6]<1.05e-6
  cor.snp.idx[is.na(cor.snp.idx)]<-FALSE
  #print(c(snp.stats.old.lfmm[,i-6]>0.730,snp.stats.young.lfmm[,i-6]>0.730))
  
  pools.varscan.old.outlier<-pooldata.subset(pools.varscan.old,snp.index=c(1:47552)[cor.snp.idx])
  fst.p.varscan.old.outlier<-compute.pairwiseFST(pools.varscan.old.outlier,method="Identity")
  pools.varscan.young.outlier<-pooldata.subset(pools.varscan.young,snp.index=c(1:47552)[cor.snp.idx])
  fst.p.varscan.young.outlier<-compute.pairwiseFST(pools.varscan.young.outlier,method="Identity")
  
  ##Old
  lm1<-mantel.partial(fst.p.varscan.old.outlier@PairwiseFSTmatrix,as.matrix(clim.D.old),fst.p.varscan@PairwiseFSTmatrix[idx,idx],method="pear",perm=9999)
  ##Young
  lm2<-mantel.partial(fst.p.varscan.young.outlier@PairwiseFSTmatrix,as.matrix(clim.D.young),fst.p.varscan@PairwiseFSTmatrix[!idx,!idx],method="pear",perm=9999)
  
  #mstat<-rbind(mstat,c(i,sum(snp.stats.old.r[,i-6]>0.730,na.rm=T),lm1$statistic,lm1$signif,sum(snp.stats.young.r[,i-6]>0.730,na.rm=T),lm2$statistic,lm2$signif))
  mstat<-rbind(mstat,c(i,sum(snp.stats.old.lfmm[,i-6]<1.05e-6,na.rm=T),lm1$statistic,lm1$signif,sum(snp.stats.young.lfmm[,i-6]<1.05e-6,na.rm=T),lm2$statistic,lm2$signif))
}
mstat<-as.data.frame(mstat)
names(mstat)<-c("clim_var","n.old","old.r","old.p","n.pl","pl.r","pl.p")


##
##Plot climate vs genetic distance (Fst) based on SNPs selected using correlations
##
traits<-c(8,9,10,12,13,14,17)

##Correlations
cor.snp.idx<-apply(snp.stats.old.r>0.730,1,any)
cor.snp.idx[is.na(cor.snp.idx)]<-FALSE

clim.D.old<-vegdist(samples[idx,traits],method="euclidean")
clim.D.young<-vegdist(samples[!idx,traits],method="euclidean",upper=T)
pools.varscan.old.outlier<-pooldata.subset(pools.varscan.old,snp.index=c(1:47552)[cor.snp.idx])
fst.p.varscan.old.outlier<-compute.pairwiseFST(pools.varscan.old.outlier,method="Identity")
pools.varscan.young.outlier<-pooldata.subset(pools.varscan.young,snp.index=c(1:47552)[cor.snp.idx])
fst.p.varscan.young.outlier<-compute.pairwiseFST(pools.varscan.young.outlier,method="Identity")

par(mfrow=c(2,2))
fst.vec<-fst.p.varscan.old.outlier@PairwiseFSTmatrix[lower.tri(fst.p.varscan.old.outlier@PairwiseFSTmatrix)]
plot(fst.vec~sqrt(as.vector(clim.D.old)),main="",xlab="Climate distance",ylab=expression("Genetic distance, F"[ST]),pch=19,col="forestgreen",las=1,ylim=c(0.0,0.15))
lm.outlier<-lm(fst.vec~sqrt(as.vector(clim.D.old)))
summary(lm.outlier)
mantel(as.dist(fst.p.varscan.old.outlier@PairwiseFSTmatrix),clim.D.old)
abline(lm.outlier,lty=3)

fst.vec<-fst.p.varscan.young.outlier@PairwiseFSTmatrix[lower.tri(fst.p.varscan.young.outlier@PairwiseFSTmatrix)]
plot(fst.vec~sqrt(as.vector(clim.D.young)),main="",xlab="Climate distance",ylab=expression("Genetic distance, F"[ST]),pch=19,col="tan2",las=1,ylim=c(0.0,0.15))
lm.outlier<-lm(fst.vec~sqrt(as.vector(clim.D.young)))
summary(lm.outlier)
mantel(as.dist(fst.p.varscan.young.outlier@PairwiseFSTmatrix),clim.D.young)
abline(lm.outlier,lty=3)

##LFMM
cor.snp.idx<-apply(snp.stats.old.lfmm<1.05e-6,1,any)
cor.snp.idx[is.na(cor.snp.idx)]<-FALSE

clim.D.old<-vegdist(samples[idx,traits],method="euclidean")
clim.D.young<-vegdist(samples[!idx,traits],method="euclidean",upper=T)
pools.varscan.old.outlier<-pooldata.subset(pools.varscan.old,snp.index=c(1:47552)[cor.snp.idx])
fst.p.varscan.old.outlier<-compute.pairwiseFST(pools.varscan.old.outlier,method="Identity")
pools.varscan.young.outlier<-pooldata.subset(pools.varscan.young,snp.index=c(1:47552)[cor.snp.idx])
fst.p.varscan.young.outlier<-compute.pairwiseFST(pools.varscan.young.outlier,method="Identity")

fst.vec<-fst.p.varscan.old.outlier@PairwiseFSTmatrix[lower.tri(fst.p.varscan.old.outlier@PairwiseFSTmatrix)]
plot(fst.vec~sqrt(as.vector(clim.D.old)),main="",xlab="Climate distance",ylab=expression("Genetic distance, F"[ST]),pch=19,col="forestgreen",las=1,ylim=c(0.0,0.15))
lm.outlier<-lm(fst.vec~sqrt(as.vector(clim.D.old)))
summary(lm.outlier)
mantel(as.dist(fst.p.varscan.old.outlier@PairwiseFSTmatrix),clim.D.old)
abline(lm.outlier,lty=3)

fst.vec<-fst.p.varscan.young.outlier@PairwiseFSTmatrix[lower.tri(fst.p.varscan.young.outlier@PairwiseFSTmatrix)]
plot(fst.vec~sqrt(as.vector(clim.D.young)),main="",xlab="Climate distance",ylab=expression("Genetic distance, F"[ST]),pch=19,col="tan2",las=1,ylim=c(0.0,0.15))
lm.outlier<-lm(fst.vec~sqrt(as.vector(clim.D.young)))
summary(lm.outlier)
mantel(as.dist(fst.p.varscan.young.outlier@PairwiseFSTmatrix),clim.D.young)
abline(lm.outlier,lty=3)


 ##Climate PCA and biplot
##Identify variables with high correlation (>0.8)
#clim.var.idx<-6+c(1:16)[!(c(1:16) %in% FindCorr(cor(samples[,c(7:22)]),cutoff=0.8))]
clim.var.idx<-traits
par(mfrow=c(1,1))
clim.pca<-prcomp(samples[idx,clim.var.idx],scale=T)
biplot(clim.pca,xlabs=rep("",15),xlim=c(-0.2,0.2),ylim=c(-0.2,0.2),expand=0.26,ylabs=clim.var[clim.var.idx-6],cex=0.75,xlab="PC1 (61.7%)",ylab="PC2 (27.4%)",las=1)
points(clim.pca$x[,c(1,2)],pch=21,bg=col.vector[idx],cex=1.5)



##Estimate RONA across old and planted populations
source("~/Downloads/doi_10.5061_dryad.866t1g1pm__v3/scripts/rona.R")
##Trait to analyse
traits<-c(8,9,10,12,13,14,17)

RONA<-NULL
for(i in traits){
  
  ##
  ##Select SNPs using correlations
  ##
  cor.snp.idx<-snp.stats.old.r[,i-6]>0.730
  cor.snp.idx[is.na(cor.snp.idx)]<-FALSE
  print(sum(cor.snp.idx))
    
  ##
  ##Select SNPs using LFMM
  ##
  ##Bonferonni correction
  #cor.snp.idx<-snp.stats.old.lfmm[,i-6]<1.05e-6
  ##Top 0.1% SNPs
  #cor.snp.idx<-snp.stats.old.lfmm[,i-6]<quantile(snp.stats.old.lfmm[,i-6],0.001,na.rm=T)
  #cor.snp.idx[is.na(cor.snp.idx)]<-FALSE
  #print(sum(cor.snp.idx))
 
  ##Calculate RONA using selected SNPs
  RONA.temp<-rona_pred(pools.af.t[,cor.snp.idx],scale(samples[,i],scale=T,center=T),
                  scale(samples4.5[,i],scale=T,center=T),significance=F,weighted=T)
  RONA<-cbind(RONA,RONA.temp[[3]]$Avg_weighted_RONA)
}
RONA<-as.data.frame(RONA)
names(RONA)<-names(samples)[traits]
RONA$Pop<-(c("Planted","Old")[idx+1])


plot_signif<-function(x,y,ydiff,nsig=1){
  lines(x=c(x,x),y=c(y-ydiff,y))
  lines(x=c(x,x+1),y=c(y,y))
  lines(x=c(x+1,x+1),y=c(y-ydiff,y))
  text(x=x+0.5,y=y+ydiff,labels=paste(rep("*",nsig),collapse=""),adj=0.5)
}

##Plot RONA values 
par(mfrow=c(2,1))
par(mar=c(4,4,1,1))
boxplot(value~Pop+variable,data=melt(RONA.r),col=rep(c("forestgreen","tan2"),7),
        xlab="Climate variables",ylab="RONA",las=1,outline=F,names=c("AIT","","CMI","","CONT","","DegD0","","DegD5","","MTC","","PCQ",""),
        ylim=c(0,0.15))
anova(lm(value~Pop+variable,data=melt(RONA.r)))

anova(lm(RONA.r[,1]~RONA.r$Pop))
plot_signif(1,0.1,0.005,2)
anova(lm(RONA.r[,2]~RONA.r$Pop))
plot_signif(3,0.09,0.005,2)
anova(lm(RONA.r[,3]~RONA.r$Pop))
plot_signif(5,0.09,0.005,2)
anova(lm(RONA.r[,4]~RONA.r$Pop))
plot_signif(7,0.1,0.005,1)
anova(lm(RONA.r[,5]~RONA.r$Pop))
plot_signif(9,0.13,0.005,2)
anova(lm(RONA.r[,6]~RONA.r$Pop))
plot_signif(11,0.09,0.005,3)
anova(lm(RONA.r[,7]~RONA.r$Pop))
plot_signif(13,0.08,0.005,3)

boxplot(value~Pop+variable,data=melt(RONA.lfmm),col=rep(c("forestgreen","tan2"),7),
        xlab="Climate variables",ylab="RONA",las=1,outline=F,names=c("AIT","","CMI","","CONT","","DegD0","","DegD5","","MTC","","PCQ",""),
        ulim=c(0,0.20))
anova(lm(value~Pop+variable,data=melt(RONA.lfmm)))

anova(lm(RONA.lfmm[,1]~RONA.lfmm$Pop))
plot_signif(1,0.17,0.005,1)
anova(lm(RONA.lfmm[,2]~RONA.lfmm$Pop))
plot_signif(3,0.13,0.005,2)
anova(lm(RONA.lfmm[,3]~RONA.lfmm$Pop))
plot_signif(5,0.11,0.005,1)
anova(lm(RONA.lfmm[,4]~RONA.lfmm$Pop))
plot_signif(7,0.1,0.005,1)
anova(lm(RONA.lfmm[,5]~RONA.lfmm$Pop))
#plot_signif(9,0.13,0.005,2)
anova(lm(RONA.lfmm[,6]~RONA.lfmm$Pop))
#plot_signif(11,0.09,0.005,3)
anova(lm(RONA.lfmm[,7]~RONA.lfmm$Pop))
#plot_signif(13,0.08,0.005,3)

