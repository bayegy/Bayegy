setwd("~/pipelines/Github/CCLasso")

# Basic example
source("R/cclasso.R");
library("gtools")
source("R/SparCC.R");
library(corrplot)
# 1. generate logistic normal variables
n <- 100;
p <- 20;
x <- matrix(rnorm(n * p), nrow = n); 
x.frac <- exp(x) / rowSums(exp((x)));
totCount <- round(runif(n = n,  min = 1000, max = 2000));
x.count <- x.frac * totCount;
# 2. run cclasso 
# using fraction
res_ccl_frac <- cclasso(x = x.frac, counts = F);
# using counts
res_ccl_count <- cclasso(x = x.count, counts = T);
# 3. run SparCC.count and SparCC.frac
res_spa_count <- SparCC.count(x = x.count);
res_spa_frac <- SparCC.frac(x = x.frac);
# 4. get the correlation matrix
{
  cat("CCLasso using fraction data:\n");
  print(round(res_ccl_frac$cor_w, 2));
  cat("CCLasso using count data:\n");
  print(round(res_ccl_count$cor_w, 2));
  cat("SparCC using fraction data:\n");
  print(round(res_spa_frac$cor.w, 2));
  cat("SparCC using count data:\n");
  print(round(res_spa_count$cor.w, 2));
}

class(x.frac)
####################16S
new<-read.table("~/Desktop/Project/VIH/Banocc/Caprisa_454_16s_top100species.txt", header = TRUE)
dim(new)
new$ID<-NULL
summary(new)
colnames(new)
#x.frac.2<-new
#x.frac.3<-as.numeric(as.matrix(x.frac.2))
#x.frac.3[x.frac.3 == 0.000000000] <- 0.000000001
#as.matrix(class(new))
res_ccl_frac <- cclasso(x = new, counts = F);
# using fraction
#res_ccl_frac <- cclasso(x = x.frac.3, counts = F);
print(round(res_ccl_frac$cor_w, 2));
rounded_corr<-round(res_ccl_frac$cor_w, 2)
colnames(rounded_corr) <- colnames(new)
rownames(rounded_corr) <- colnames(new)
rounded_corr

rounded_pval<-round(res_ccl_frac$p_vals, 2)
colnames(rounded_pval) <- colnames(new)
rownames(rounded_pval) <- colnames(new)
rounded_pval

write.table(rounded_corr , file = "Caprisa_454_16s_top100species.cor_w.txt", row.names = T, col.names=T, sep="\t")
write.table(rounded_pval , file = "Caprisa_454_16s_top100species.p_vals.txt", row.names = T, col.names=T, sep="\t")

###############Kraken
new.k<-read.table("~/Desktop/Project/VIH/Banocc/Caprisa_kraken_top100species.txt", header = TRUE)
dim(new.k)
#colnames.k <-colnames(new.k)
new.k$ID<-NULL
summary(new.k)
colnames(new.k)
new.k[new.k == 0] <- 0.000000001
res_ccl_frac.k <- cclasso(x = new.k, counts = F);

print(round(res_ccl_frac.k$cor_w, 2));
rounded_corr.k<-round(res_ccl_frac.k$cor_w, 2)
colnames(rounded_corr.k) <- colnames(new.k)
rownames(rounded_corr.k) <- colnames(new.k)
rounded_corr.k

rounded_pval.k<-round(res_ccl_frac.k$p_vals, 2)
colnames(rounded_pval.k) <- colnames(new.k)
rownames(rounded_pval.k) <- colnames(new.k)
rounded_pval.k

write.table(rounded_corr.k , file = "Caprisa_kraken_top100species.cor_w.txt", row.names = T, col.names=T, sep="\t")
write.table(rounded_pval.k , file = "Caprisa_kraken_top100species.p_vals.txt", row.names = T, col.names=T, sep="\t")


install.packages("corrplot")
library(corrplot)

names<-colnames(new.k)
class(names)

names.short<-substr(names, 1, 200)
names.short

y <- strsplit(names,"s__",fixed=TRUE)
names.short<-lapply(y,FUN=function(x){paste(x[2],"",sep="")})

as.character(names.short)
class(names.short)

colnames(rounded_corr.k)<-names.short
rownames(rounded_corr.k)<-names.short
colnames(rounded_pval.k)<-names.short
rownames(rounded_pval.k)<-names.short

corrplot(rounded_corr.k, type="full", order="original", p.mat = rounded_pval.k, sig.level = 0.01, insig = "blank", tl.cex=0.5, tl.col="black")
################################
new.m<-read.table("~/Desktop/Project/VIH/Banocc/Caprisa_Metaphlan_specieslevel.txt", header = TRUE)
dim(new.m)
new.m$ID<-NULL
summary(new.m)
colnames(new.m)
new.m[new.m == 0] <- 0.000000001
res_ccl_frac.m <- cclasso(x = new.m, counts = F);

print(round(res_ccl_frac.m$cor_w, 2));
rounded_corr.m<-round(res_ccl_frac.m$cor_w, 2)
colnames(rounded_corr.m) <- colnames(new.m)
rownames(rounded_corr.m) <- colnames(new.m)
rounded_corr.m

rounded_pval.m<-round(res_ccl_frac.m$p_vals, 2)
colnames(rounded_pval.m) <- colnames(new.m)
rownames(rounded_pval.m) <- colnames(new.m)
rounded_pval.m

rounded_corr.m[is.na(rounded_corr.m)] <- 0.000000001

write.table(rounded_corr.m , file = "Caprisa_metaphlan_top100species.cor_w.txt", row.names = T, col.names=T, sep="\t")
write.table(rounded_pval.m , file = "Caprisa_metaphlan_top100species.p_vals.txt", row.names = T, col.names=T, sep="\t")

names<-colnames(new.m)
class(names)

y <- strsplit(names,"s__",fixed=TRUE)
names.short<-lapply(y,FUN=function(x){paste(x[2],"",sep="")})

colnames(rounded_corr.m)<-names.short
rownames(rounded_corr.m)<-names.short
colnames(rounded_pval.m)<-names.short
rownames(rounded_pval.m)<-names.short

corrplot(rounded_corr.m, type="full", order="original", p.mat = rounded_pval.m, sig.level = 0.05, insig = "blank", tl.cex=0.3, tl.col="black")
################################

new.s<-read.table("~/Desktop/Project/VIH/Banocc/Caprisa_454_16s_top100species.2.txt", header = TRUE)
dim(new.s)
new.s$ID<-NULL
summary(new.s)
colnames(new.s)
new.s[new.s == 0] <- 0.000000001
res_ccl_frac.s <- cclasso(x = new.s, counts = F);

print(round(res_ccl_frac.s$cor_w, 2));
rounded_corr.s<-round(res_ccl_frac.s$cor_w, 2)
colnames(rounded_corr.s) <- colnames(new.s)
rownames(rounded_corr.s) <- colnames(new.s)
rounded_corr.s

rounded_pval.s<-round(res_ccl_frac.s$p_vals, 2)
colnames(rounded_pval.s) <- colnames(new.s)
rownames(rounded_pval.s) <- colnames(new.s)
rounded_pval.s

rounded_corr.s[is.na(rounded_corr.s)] <- 0.000000001

write.table(rounded_corr.s , file = "Caprisa_16S_top100species.cor_w.txt", row.names = T, col.names=T, sep="\t")
write.table(rounded_pval.s , file = "Caprisa_16S_top100species.p_vals.txt", row.names = T, col.names=T, sep="\t")

names<-colnames(new.s)
class(names)
names.short<-names

#y <- strsplit(names,'\\.',fixed=FALSE)
#head(y)
#n<-length(y[[1]])
#n<-length(y)
#y[n]
#y[1,2]
#class(y)
#names.short<-lapply(y,FUN=function(x){paste(x[8],"",sep="")})
#names.short[[101]]<-'Others'

colnames(rounded_corr.s)<-names.short
rownames(rounded_corr.s)<-names.short
colnames(rounded_pval.s)<-names.short
rownames(rounded_pval.s)<-names.short

corrplot(rounded_corr.s, type="full", order="original", p.mat = rounded_pval.m, sig.level = 0.5, insig = "blank", tl.cex=0.5, tl.col="black")
################################

?corrplot()

################################20180617 to do work 5 prevotella from Andrew

new.s<-read.table("~/Desktop/Project/Andrew/Sparcc/AllPrevotella_GreaterThan1Percent.csv", header = TRUE, sep=",")

#new.s<-read.table("~/Desktop/Project/Andrew/Sparcc/Caprisa_Prevotella_Kraken_NEW.csv", header = TRUE, sep=",")
#new.s<-read.table("~/Desktop/Project/Andrew/Sparcc/Table_S4_Relative_PrevotellaGreaterthan1%25.csv", header = TRUE, sep=",")
#new.s<-read.table("~/Desktop/Project/Andrew/Sparcc/Caprisa_454_Prevotella_NEW.csv", header = TRUE, sep=",")
#new.s<-read.table("~/Desktop/Project/Andrew/Sparcc/Caprisa_Prevotella_Metaphlan_NEW.csv", header = TRUE, sep=",")






dim(new.s)
#new.s$ID<-NULL
summary(new.s)
colnames(new.s)
new.s[new.s == 0] <- 0.000000001
res_ccl_frac.s <- cclasso(x = new.s, counts = F);

print(round(res_ccl_frac.s$cor_w, 2));
rounded_corr.s<-round(res_ccl_frac.s$cor_w, 2)
colnames(rounded_corr.s) <- colnames(new.s)
rownames(rounded_corr.s) <- colnames(new.s)
rounded_corr.s

rounded_pval.s<-round(res_ccl_frac.s$p_vals, 2)
colnames(rounded_pval.s) <- colnames(new.s)
rownames(rounded_pval.s) <- colnames(new.s)
rounded_pval.s

rounded_corr.s[is.na(rounded_corr.s)] <- 0.000000001

#write.table(rounded_corr.s , file = "Caprisa_16S_top100species.cor_w.txt", row.names = T, col.names=T, sep="\t")
#write.table(rounded_pval.s , file = "Caprisa_16S_top100species.p_vals.txt", row.names = T, col.names=T, sep="\t")

names<-colnames(new.s)
class(names)
names.short<-names

#y <- strsplit(names,'\\.',fixed=FALSE)
#head(y)
#n<-length(y[[1]])
#n<-length(y)
#y[n]
#y[1,2]
#class(y)
#names.short<-lapply(y,FUN=function(x){paste(x[8],"",sep="")})
#names.short[[101]]<-'Others'

colnames(rounded_corr.s)<-names.short
rownames(rounded_corr.s)<-names.short
colnames(rounded_pval.s)<-names.short
rownames(rounded_pval.s)<-names.short

corrplot(rounded_corr.s, type="full", p.mat = rounded_pval.s, sig.level = 0.05, insig = "blank", tl.cex=0.5, tl.col="black", order= "hclust", hclust.method= "complete")
#corrplot(rounded_corr.s, type="full", p.mat = rounded_pval.s, sig.level = 0.05, insig = "blank", tl.cex=0.5, tl.col="black", order= "original")



