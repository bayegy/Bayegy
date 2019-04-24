library(rjson)
library(RColorBrewer)

get_colors=function(category=NULL, mapping_file=NULL){
    colors_plan_path=Sys.getenv("COLORS_PLAN_PATH")
    if(colors_plan_path!=""&!is.null(category)){
        colors_plan=fromJSON(file = colors_plan_path)
        return(colors_plan[category][[1]])
    }else{
        pallet<-c(rainbow(3),brewer.pal(12,"Set3"))
        if(!is.null(mapping_file)&!is.null(category)){
            map<-read.table(mapping_file,quote="",row.names = 1,na.strings = "",comment.char="",check.names=F,stringsAsFactors=F, header = TRUE, sep = "\t")
            number_of_group<-length(unique(map[category][,1]))
            return(pallet[1:number_of_group])
        }else{
            return(pallet)
        }
    }
}



