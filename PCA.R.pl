#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);
use FindBin qw($Bin);
use Getopt::Long;
use lib "$Bin/";
use PATHWAY;
my $cfg = "$Bin/Pathway_cfg.txt";
(-s $cfg) || die"error: can't find config file: $cfg, $!\n";
my ($R0)=get_pathway($cfg,[qw(R)]);
my $T;
GetOptions("R:s"=>\$R0,"T"=>\$T);
die "perl $0 <exp_data> <cutoff(recomend 0.2)> <outdir> [-R R_PATH] [-T(to change input matrix)]\n" unless @ARGV>=3;
my ($exp_data, $cutoff, $outdir) = @ARGV;
my $R ||= $R0;
(-d $outdir) || mkdir($outdir);
$outdir =abs_path($outdir);
if($T){
    my $bname = (split/\//,$exp_data)[-1];
    system"perl $Bin/taxa_table.pl $exp_data $outdir/$bname";
    $exp_data = "$outdir/$bname";
}

my $file =<< "EOF";

	library(gmodels)
	setwd(\"$outdir\")
	dat = read.table(\"$exp_data\", head=T, row.names=1,sep="\\t")
	
b <- matrix(0,nrow = nrow(dat), ncol = ncol(dat))
for(i in 1:ncol(dat)){
        b[,i] = dat[ ,i]/sum(dat[ ,i])
}
colnames(b) <- colnames(dat)
rownames(b) <- rownames(dat)
data <- t(b)
#### PCA 
    data.pca <- fast.prcomp(data,retx=T,scale=F,center=T)
    a <- summary(data.pca)
    b <- a[4]\$importance
	PC1 <-as.numeric(data.pca\$x[,1])
	PC2 <-as.numeric(data.pca\$x[,2])
#	PC3 <-as.numeric(data.pca\$x[,3])
	write.csv(data.pca\$x[,1:2],file="$outdir/PCA.csv")
    pc1 <- as.numeric(sprintf("%.3f",b[2,1]))*100
    pc2 <- as.numeric(sprintf("%.3f",b[2,2]))*100
#    pc3 <- as.numeric(sprintf("%.3f",b[2,3]))*100
	rownames(data.pca\$rotation)=rownames(dat)
	bb<- data.pca\$rotation[abs(data.pca\$rotation[,1])> $cutoff|abs(data.pca\$rotation[,2])>$cutoff,]
	lab=colnames(dat)
	lab = gsub("X","",lab)
	y.length <- ncol(dat)
	legend.ncol<-ceiling(y.length/50)
    cols<- c("pink", "orange", "green", "cyan", "blue", "darkblue", "purple","seagreen","steelblue","red","DarkTurquoise","Sienna","Chocolate","BlueViolet","Magenta", "brown", "gray", "darkred", "darkmagenta","darkcyan")
	pchs <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
#### Plot,pdf
   pdf("$outdir/PCA-2D.pdf",height=7,width=14)
   par(mfrow = c(1,2))
   xlab <- paste("PC1 ", pc1,"%",sep="")
   ylab <- paste("PC2 ", pc2,"%",sep="")
   main <- "PCA-PC1 vs PC2"
   plot(PC1,PC2, xlab = xlab, ylab = ylab,cex = 0.8, main=main, type="n")
   abline(v = 0, h = 0, lty = 3, col="grey")
   points(PC1, PC2, col=cols,cex = 0.8, pch = pchs,bg=pchs)
#   text(PC1, PC2, adj =-0.2 ,labels=lab,pos = NULL, offset = 0.4, vfont = NULL,cex =1.2, col = NULL, font = NULL)
#   arrows(0,0,(bb[,1])*0.12,(bb[,2])*0.12,code=2,length=0.05)
#   text((bb[,1])*0.15,(bb[,2])*0.15,labels=rownames(bb),cex =0.8,font=1)
   barplot(0,axes=F,col="white",border=NA,axisnames=F)
   legend("topleft", legend=lab, col=cols, pch = pchs, cex=0.7,bty="n",pt.bg=cols,ncol =legend.ncol)
   dev.off()
#### Plot,png
   png("$outdir/PCA-2D.png",height=700,width=1400,type="cairo")
   par(mfrow = c(1,2))
   xlab <- paste("PC1 ", pc1,"%",sep="")
   ylab <- paste("PC2 ", pc2,"%",sep="")
   main <- "PCA-PC1 vs PC2"
   plot(PC1,PC2, xlab = xlab, ylab = ylab,cex = 0.8, main=main, type="n")
   abline(v = 0, h = 0, lty = 3, col="grey")
   points(PC1, PC2, col=cols,cex = 0.8, pch = pchs,bg=pchs)
#   text(PC1, PC2, adj =-0.2 ,labels=lab,pos = NULL, offset = 0.4, vfont = NULL,cex =1.2, col = NULL, font = NULL)
#   arrows(0,0,(bb[,1])*0.12,(bb[,2])*0.12,code=2,length=0.05)
#   text((bb[,1])*0.15,(bb[,2])*0.15,labels=rownames(bb),cex =0.8,font=1)
   barplot(0,axes=F,col="white",border=NA,axisnames=F)
  legend("topleft", legend=lab, col=cols, pch = pchs, cex=0.7,bty="n",pt.bg=cols,ncol =legend.ncol)
   dev.off()
EOF

open OUT, ">$outdir/PCA.R" or die $!;
print OUT $file;
close OUT;

`$R -f $outdir/PCA.R`;
