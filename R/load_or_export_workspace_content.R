#' Generates description for a flowFrame
#'
#' @param ff A flowframe
#' @return A flowframe with description consistent with pData(ff@parameters)
#' @export

generate_description<-function(ff){
    res=ff
    res.desc=pData(res@parameters)
    new.names=sprintf("$P%sS",1:nrow(res.desc))
    new.values=setNames(res.desc[,"desc"],new.names)
    new.values[is.na(new.values)]=" "
    new.values=as.list(new.values)
    res@description[names(res@description)%in%names(new.values)]=NULL
    res@description=c(res@description,new.values)

    new.names=sprintf("$P%sR",1:nrow(res.desc))
    new.values=setNames(ceiling(res.desc[,"maxRange"]-res.desc[,"minRange"]),new.names)
    new.values[is.na(new.values)]=" "
    new.values=as.list(new.values)
    res@description[names(res@description)%in%names(new.values)]=NULL
    res@description=c(res@description,new.values)

    new.names=sprintf("$P%sE",1:nrow(res.desc))
    new.values=setNames(rep("0,0",length(new.names)),new.names)
    res@description[names(res@description)%in%names(new.values)]=NULL
    res@description=c(res@description,new.values)

    new.names=sprintf("$P%sB",1:nrow(res.desc))
    new.values=setNames(rep("32",length(new.names)),new.names)
    res@description[names(res@description)%in%names(new.values)]=NULL
    res@description=c(res@description,new.values)

    new.names=sprintf("$P%sG",1:nrow(res.desc))
    new.values=setNames(rep("1",length(new.names)),new.names)
    res@description[names(res@description)%in%names(new.values)]=NULL
    res@description=c(res@description,new.values)
    res
}

#' Wrapper to flowCore:::write.flowSet to preserve channel annotation
#'
#' @param fs A flowSet
#' @param directory Where to save the flowFrames of the flowSet
#' @export

en.write_flowset<-function(fs,directory){
    dir.create(directory)
    sapply(fs@phenoData$name,function(x){
        res=fs[[x]]

        res=generate_description(res)

        write.FCS(res,filename=file.path(directory,x))
    })
}
