#!/usr/bin/env Rscript
library(optparse)
#######arguments
option_list <- list(
  make_option(c("-i", "--input"),metavar="path", dest="table",help="abundance table. Required",default=NULL),
  make_option(c("-m", "--map"),metavar="path",dest="map", help="Specify the path of mapping file. Required",default=NULL),
  make_option(c("-c", "--category"),metavar="string",dest="category", help="Category to compare. Required",default=NULL),
  make_option(c("-n", "--numerator"),metavar="string",dest="numerator", help="group as numerator for fold change",default=NULL),
  make_option(c("-d", "--denominator"),metavar="str", dest="denominator",help="group as denominator for fold change",default=""),
  make_option(c("-a", "--adjustp"),metavar="logical", dest="adjustp",help="use adjusted p value in volcano plot", default=FALSE),
  make_option(c("-P", "--threshold-p"),metavar="float", dest="threshold_p",help="threshold of p value", default=0.05),
  make_option(c("-F", "--threshold-fc"),metavar="float", dest="threshold_fc",help="threshold of log2FoldChange", default=0.05),
  make_option(c("-t", "--top"),metavar="int", dest="top",help="most significant features needed to be labeled, default is 0, meaning all significant features will be labeled",default=0),
  make_option(c("-f", "--filter-freq"),metavar="float", dest="min_freq",help="features absence frequency among samples blow this threshold will be excluded from further analysis",default=0),
  make_option(c("-p", "--prefix"),metavar="str", dest="prefix",help="The prefix of output files, default if null",default=""),
  make_option(c("-o", "--output"),metavar="directory",dest="out_dir", help="Specify the directory of output files. default=./",default="./")
)

opt <- parse_args(OptionParser(option_list=option_list,description = "This script is used to compare the predicted pathway of level 3"))


library(DESeq2)
library(ggplot2)
library(ggrepel)

choose<-function(condition,choice1,choice2){
  if(condition){
    return(choice1)
  }else{
    return(choice2)
  }
}

#'util
or<-function(choice1,choice2){
  if(length(choice1)==1){
      if(is.na(choice1)|choice1==FALSE){
        return(choice2)
      }
  }else{
      if(length(choice1)==0){
        return(choice2)
      }
  }
  return(choice1)
}


plot_volcano<-function(data, out_img="volcano.pdf", axis=c(1,2), top=NA, threshold=c(2, 1.3), threshold_both=TRUE, label="rownames", x_lab=NA, y_lab=NA, bg_colors=c("#E1FFFF05", "#99FF9905", "#FFCC9905"), size="integ"){
    xy<-choose(is.numeric(axis), colnames(data)[axis], axis)
    lbl<-choose(is.numeric(label), choose(label, colnames(data)[label], "rownames"), label)
    sz <- choose(is.numeric(size), colnames(data)[size], size)
    x_index<-which(colnames(data)%in%xy[1])
    y_index<-which(colnames(data)%in%xy[2])
    colnames(data)[x_index]<-"plot_volcano_x"
    colnames(data)[y_index]<-"plot_volcano_y"
    span_x<-0.5*sd(data[,x_index])
    span_y<-0.5*sd(data[,y_index])
    d_max_x<-max(abs(data[,x_index]))+span_x
    d_max_y<-max(abs(data[,y_index]))+span_y
    d_min_x<-min(data[, x_index])
    d_min_y<-min(data[, y_index])
    d_min_x<-choose(d_min_x<0, d_min_x-span_x, d_min_x)
    d_min_y<-choose(d_min_y<0, d_min_y-span_y, d_min_y)
    ifelse(lbl=="rownames", data$plot_volcano_lbl<-rownames(data), colnames(data)[colnames(data)%in%lbl]<-"plot_volcano_lbl")
    if(size=="integ"){
      data$plot_volcano_size <- (abs(data$plot_volcano_x) + abs(data$plot_volcano_y))**2
    }else{
      if(!is.na(size)){
        colnames(data)[colnames(data)%in%sz]<-"plot_volcano_size"
      }
    }
    data <- data[order(abs(data[, y_index]), abs(data[, x_index]), decreasing = T),]
    label_data<-choose(is.na(top),
                       choose(threshold_both,
                              data[(abs(data[x_index])>or(threshold[1], 0))&(abs(data[y_index])>or(threshold[2], 0)), ],
                              data[(abs(data[x_index])>or(threshold[1], 0))|(abs(data[y_index])>or(threshold[2], 0)), ]),
                       head(data, top))

    x_lab<-or(x_lab, xy[1])
    y_lab<-or(y_lab, xy[2])
    xin_lab<-paste(x_lab,"=",threshold[1])
    yin_lab<-paste(y_lab,"=",threshold[2])
    # xin_coord <- c(threshold[1]+span_x, span_y)
    # yin_coord <- c(span_x, threshold[2]+span_y)

    in_lab_data <- data.frame(x=c(threshold[1], d_min_x+span_x*2),
                              y=c(d_min_y+span_y, threshold[2]),
                              label=c(xin_lab, yin_lab))
    data$plot_volcano_dot_colors <- "#222222"
    data$plot_volcano_dot_colors[rownames(data)%in%rownames(label_data)] <- "#DC143C"
    p<-ggplot(data=data,aes(x=plot_volcano_x, y=plot_volcano_y)) +
      geom_hline(yintercept=threshold[2],linetype="dashed") +
      choose(d_min_y>=0, theme(), geom_hline(yintercept=-threshold[2],linetype="dashed")) +
      geom_vline(xintercept=threshold[1],linetype="dashed") +
      choose(d_min_x>=0, theme(), geom_vline(xintercept=-threshold[1],linetype="dashed")) +
      choose(d_min_x>=0, theme(), geom_vline(xintercept = 0)) +
      choose(d_min_y>=0, theme(), geom_hline(yintercept = 0)) +
      choose(is.na(threshold[1]), theme(), geom_rect(xmin=threshold[1], xmax=d_max_x, ymin=d_min_y, ymax=d_max_y, fill=bg_colors[1])) +
      choose(d_min_x>=0, theme(), geom_rect(xmin=-threshold[1], xmax=-d_max_x, ymin=d_min_y, ymax=d_max_y, fill=bg_colors[1])) +
      choose(is.na(threshold[2]), theme(), geom_rect(ymin=threshold[2], ymax=d_max_y, xmin=d_min_x, xmax=d_max_x, fill=bg_colors[2])) +
      choose(d_min_y>=0, theme(), geom_rect(ymin=-threshold[2], ymax=-d_max_y, xmin=d_min_x, xmax=d_max_x, fill=bg_colors[2])) +
      choose(any(is.na(threshold)), theme(), geom_rect(ymin=threshold[2], ymax=d_max_y, xmin=threshold[1], xmax=d_max_x, fill=bg_colors[3])) +
      choose(any(is.na(threshold))|d_min_x>=0, theme(), geom_rect(ymin=threshold[2], ymax=d_max_y, xmin=-threshold[1], xmax=-d_max_x, fill=bg_colors[3])) +
      choose(any(is.na(threshold))|d_min_y>=0, theme(), geom_rect(ymin=-threshold[2], ymax=-d_max_y, xmin=threshold[1], xmax=d_max_x, fill=bg_colors[3])) +
      choose(any(is.na(threshold))|d_min_y>=0|d_min_x>=0, theme(), geom_rect(ymin=-threshold[2], ymax=-d_max_y, xmin=-threshold[1], xmax=-d_max_x, fill=bg_colors[3])) +
      choose(is.na(size), geom_point(size=3, color=data$plot_volcano_dot_colors, alpha=0.3), geom_point(aes(size=plot_volcano_size), color=data$plot_volcano_dot_colors, alpha=0.3)) +
      geom_text_repel(data=label_data, aes(x=plot_volcano_x,y=plot_volcano_y,label=plot_volcano_lbl), color='black', size=3) +
      geom_text(data = in_lab_data, mapping = aes(x=x,y=y,label=label), size=4, color='red') +
      theme_bw() + xlab(x_lab) + ylab(y_lab) + labs(size = size) +
      choose(size=="integ", guides(size=FALSE), theme()) +
      scale_y_continuous(limits=c(d_min_y, d_max_y), expand = c(0, 0)) +
      scale_x_continuous(limits=c(d_min_x, d_max_x), expand = c(0, 0)) +
      theme(panel.grid=element_blank(),
            axis.line = element_line(),
            panel.border =  element_blank(),
            text = element_text(size = 15))
    ggsave(plot = p, out_img, dpi=300, height = 8, width = choose(is.na(size), 8, 10))
}



run_deseq2 <-function(table, map, category, numerator, denominator, adjustp=FALSE, top=0, prefix="", min_freq=0.2, threshold_p=0.05, threshold_fc=2, out_dir="./"){
    if(!dir.exists(out_dir)){dir.create(out_dir,recursive = T)}
    out_dir<-paste(out_dir,"/",prefix,sep="")
    table<-read.table(table, quote="", row.names = 1, na.strings = "", comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
    map<-read.table(map, quote="", row.names = 1, na.strings = "", comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
    table <- table[rownames(map)]
    freq <- apply(table, 1, function(x){return(sum(x!=0, na.rm=TRUE)/length(x))})
    table <- table[freq>min_freq,]
    table <- table + 1
    dds <- DESeqDataSetFromMatrix(countData = table, colData = map, design = as.formula(paste("~", category)))
    dds[[category]] <- relevel(dds[[category]], ref = denominator)
    suppressWarnings(dds <- try(DESeq(dds, quiet = TRUE), silent = TRUE))
    if (inherits(dds, "try-error")) {
        # If the parametric fit failed, try the local.
        suppressWarnings(dds <- try(DESeq(dds, fitType = "local", quiet = TRUE),
        silent = TRUE))
    if (inherits(dds, "try-error")) {
        # If local fails, try the mean
        suppressWarnings(dds <- try(DESeq(dds, fitType = "mean", quiet = TRUE),
            silent = TRUE))
        }
    if (inherits(dds, "try-error")) {
        # If still bad, quit with error.
        return(NULL)
        }
    }
    res <- results(dds, contrast = c(category, numerator, denominator))
    resOrdered <- res[order(res$padj), ]
    df <- data.frame(resOrdered)
    df$`-log10(pvalue)` <- -log(df$pvalue, 10)
    df$`-log10(padj)` <- -log(df$padj, 10)
    write.table(data.frame(Features=rownames(df), df, check.names=FALSE), paste(out_dir, "DESeq2.xls", sep = ""), sep="\t", quote=F, row.names=F)
    df <- na.omit(df)
    p_col <- ifelse(adjustp, "-log10(padj)", "-log10(pvalue)")
    top <- ifelse(top==0, NA, top)
    threshold_p = round(-log(threshold_p, 10), 3)
    plot_volcano(df, out_img=paste(out_dir, "DESeq2.pdf", sep = ""), axis=c("log2FoldChange", p_col), threshold=c(threshold_fc, threshold_p), top=top)
}


run_deseq2(opt$table, opt$map, opt$category, opt$numerator, opt$denominator,
  as.logical(opt$adjustp), as.integer(opt$top), opt$prefix, as.double(opt$min_freq),
  threshold_p=as.double(opt$threshold_p), threshold_fc=as.double(opt$threshold_fc),
  out_dir=opt$out_dir)


