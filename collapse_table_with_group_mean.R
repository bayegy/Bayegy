library(optparse)
library(stringr)
#######arguments
option_list <- list( 
  make_option(c("-i", "--input"),metavar="path", dest="spe",help="Commentless table, table with header at first line. If This one is passed, -f will be ignored",default=NULL),
  make_option(c("-f", "--feature-table"),metavar="path", dest="otu",help="Abundance table with comment at first line. ",default=NULL),
  make_option(c("-t", "--classified-stat-table"),metavar="path", dest="stat",help="path of classified_stat_relative.xls",default=NULL),
  make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
  make_option(c("-r", "--collapse-map"),metavar="logical",dest="cmap", help="Weather to collapse mapping file",default=FALSE),
  make_option(c("-c", "--category"),metavar="string",dest="group", help="Category to calculate mean",default=NULL),
  make_option(c("-e", "--exclude"),metavar="string",dest="ex", help="Specify the numeric variables excluded from calculate mean in mapping file",default="none"),
  make_option(c("-s", "--scale"),metavar = "logical",dest="scale", help="IF TRUE, the sum abundance of each sample will be saled to 1 (relative abundance).",default="F"),
  make_option(c("-o", "--output"),metavar="directory",dest="out", help="Specify the directory of output files. default=./",default="./")
)

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to compute mean of abundance table according to a index in map"))
if(!dir.exists(opt$out)){dir.create(opt$out,recursive = T)}

if(!is.null(opt$spe)){
  spe<-read.table(opt$spe, quote="", header = T, comment.char = "",sep = "\t",check.names = F,stringsAsFactors = F)
#  spe<-read.table("otu_table.Genus.absolute.txt",header = T, comment.char = "",sep = "\t",check.names = F,stringsAsFactors = F)
}else{
  if(!is.null(opt$otu)){
    spe<-read.table(opt$otu,quote="", skip = 1, header = T, comment.char = "",sep = "\t",check.names = F,stringsAsFactors = F)
    }
}

map<-read.table(opt$map, quote="", row.names = 1,na.strings = "",header = T, comment.char = "",sep = "\t",check.names = F,stringsAsFactors = F)
#map<-read.table("mapping_file.txt",row.names = 1,na.strings = "",header = T, comment.char = "",sep = "\t",check.names = F,stringsAsFactors = F)

ifspe<-tryCatch({lspe<-ncol(spe);map<-map[match(colnames(spe)[-c(1,lspe)],rownames(map)),];TRUE},error=function(e){write("####no taxa table detected###",stdout());return(FALSE)})


group<-map[opt$group][,1]
#group<-map["Group1"][,1]


if(!is.null(opt$stat)){
  stat<-read.table(opt$stat, row.names=1, quote="", na.strings = "",header = T, comment.char = "",sep = "\t",check.names = F,stringsAsFactors = F)

  stat<-stat[match(rownames(map),rownames(stat)),]

  stat<-apply(stat,2,function(y){tapply(y,group,mean)})

  stat<-data.frame(Groups=rownames(stat),stat,check.names=F)

  write.table(stat,paste(opt$out,"/",opt$group,"_classified_stat_relative.xls",sep = ""),na = "",sep = "\t",row.names = F,quote = F)

}


if(as.logical(opt$cmap)){
  ex<-str_split(opt$ex,",")[[1]]
  notstr=c()
  for(i in 1:length(map)){
    notstr[i]=is.numeric(map[,i])
  }
  if(ex[1]!="none"){
    notstr<-notstr&(!colnames(map)%in%ex)
  }

  strmap<-apply(map[!notstr],2,function(y){tapply(y,group,function(z){z[1]})})
  nummap<-apply(map[notstr],2,function(y){tapply(y,group,mean)})


  if(!is.null(dim(nummap))){
    outmap<-data.frame(strmap,nummap,check.names = F,check.rows = T)
  }else{
    outmap<-data.frame(strmap,check.names = F)
  }
  outmap<-outmap[!colnames(outmap)=="Description"]
  outmap<-data.frame(id=rownames(outmap),outmap,Description=rownames(outmap))
  colnames(outmap)[1]<-"#SampleID"

  write.table(outmap,paste(opt$out,"/","map_collapsed_with_group_mean.txt",sep = ""),na = "",sep = "\t",row.names = F)

}


if(ifspe){
  outspe<-data.frame(t(apply(spe[,-c(1,lspe)],1,function(y){tapply(y,group,mean)})))
  if(as.logical(opt$scale)){outspe<-apply(outspe,2,function(x){x/sum(x)})}
  outspe<-data.frame(spe[,1],outspe,spe[,lspe])
  colnames(outspe)[c(1,ncol(outspe))]<-colnames(spe)[c(1,lspe)]
  if(!is.null(opt$spe)){
    write.table(outspe,paste(opt$out,"/","abundance_table_collapsed_with_group_mean.txt",sep = ""),na="",sep = "\t",row.names = F,quote = F)
  }else{
    write("#Abundance table collapsed with group mean",paste(opt$out,"/","abundance_table_collapsed_with_group_mean.txt",sep = ""))
    suppressWarnings(write.table(outspe,paste(opt$out,"/","abundance_table_collapsed_with_group_mean.txt",sep = ""),na="",sep = "\t",row.names = F,append = T,quote = F))
  }
}
