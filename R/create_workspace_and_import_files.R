#' Reads a flowset and converts to a list with 1) a big expression matrix with targets as colnames 2) a named vector of matrix rownames to samples 3) a named vector of colnames to channels
#' @export

en.read_flowset<-function(files,filter=TRUE,compensate=FALSE,...){
    if(length(files)==1&!filter){
        fs_src=as(read.FCS(files,...),"flowSet")
    } else if(length(files)==1&filter){
        fs_src=as(read.FCS(filter_fcs(files),...),"flowSet")
    } else if(filter){
        fs_src=flowCore::read.flowSet(filter_fcs(files),...)
    } else {
        fs_src=flowCore::read.flowSet(files,...)
    }

    if(compensate){
        w=unlist(fsApply(fs_src,function(x){"SPILL"%in%names(x@description)}))
        if(any(!w)){
            error(paste("Compensation matrix missing from files",names(fs_src[w])))
        } else {
            fs_src=fsApply(fs_src,function(x){
                compensate(x,x@description$`SPILL`)
            })
        }
    }
    
    ## Flowset to matrix
    xp_src=as.matrix(do.call(rbind,fsApply(fs_src,exprs,simplify=F)))
    rownames(xp_src)=as.character(1:nrow(xp_src))

    samples=unlist(fsApply(fs_src,nrow,simplify=F))
    samples=rep(names(samples),times=samples)
    events_to_samples=setNames(samples,as.character(1:nrow(xp_src)))

    targets=fs_src[[1]]@parameters$desc
    targets[is.na(targets)]=fs_src[[1]]@parameters$name[is.na(targets)]

    colnames_to_channels=setNames(colnames(xp_src),targets)
    colnames(xp_src)=targets

    return(list(fs_src=fs_src,xp_src=xp_src,channels.code=colnames_to_channels,events.code=events_to_samples))
}

#' Use regular expressions to only filter files ending in .fcs from a vector of files
#'
#' @param files Files to filter
#' @return A character vector of .fcs files
#' @export

filter_fcs=function(files){
    files[grep(pattern="^.*.(fcs)$",x=files,ignore.case=T)]
}

#' Create a 'names_template.csv' file from fcs files in the /fcs/config folder
#' @details The template can be copied, edited to select relevant variables, and saved to be reused in future functions such as tSNE, ONESENSE etc...
#' @return NULL
#' @examples create_namesfile_template()
#' @export

create_namesfile_template <- function(fs=fs_src,path="./config/names_template.csv"){
    require(flowCore)
    require(digest)
    checksums=flowCore::fsApply(fs,function(ff){
        digest::digest(pData(ff@parameters[,c("name","desc")])) #Make sure all name and desc are equal across the flowSet
    })
    n.checksums=table(checksums[,1])
    if(length(n.checksums)>1){
        inconsistent_checksums=names(n.checksums)[n.checksums<max(n.checksums)]
        inconsistent_files=rownames(checksums[checksums[,1]%in%inconsistent_checksums,])
        stop("The following files' channels contained inconsistent Metals and/or targets:",inconsistent_files,sep="\n")
    } else {
        res=pData(fs[[1]]@parameters[,c("desc","name")])
        colnames(res)=c("Target","Reporter")
        res=cbind(res,Select=rep("y",nrow(res)),stringsAsFactors=F)
        res[sapply(res[,"Target"],function(x){x=="BC"|is.na(x)}),"Select"]=""
        write.table(res,file=path,sep=",",row.names=FALSE)
    }
}
