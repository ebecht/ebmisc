#' Aggregates matrix or data frame
#'
#' @param df A matrix or data frame
#' @param groups A vector of groups (discrete values)
#' @param fun A function that aggregates a vector into objects of length 1
#' @param margin If 1, aggregates rows, if 2 aggregates columns. Defaults to 1
#' @param ... passed to fun
#' @return A data.frame with aggregated rows or columns
#' @export

aggregate_df=function(df,groups,fun=mean,margin=1,...){
    if(length(groups)!=dim(df)[margin]){
        stop("Size of 'groups' vector is different that the one of the specified data margin")
    }
    if(is.data.frame(df)){
        if(margin==2){
            df=as.data.frame(t(df))
        } else {
            df=as.data.frame(df)
        }
        df=split(df,groups)
    } else if(is.matrix(df)){
        if(margin==2){
            df=t(df)
        }
        df=split_matrix(df,groups,byrow=TRUE)
    }
    
    res=do.call(
        rbind,
        lapply(
            df,
            function(x){
                apply(x,2,fun,...)
            }
        )
    )

    if(margin==2){
        return(t(res))
    } else {
        return(res)
    }
}

#' Wrapper to locator to properly plot segments on the fly
#' @export

en.locator<-function(){
  input=TRUE
  x=vector()
  y=vector()
  while(!is.null(input)){
    input=locator(1)
    x=c(x,input$x)
    y=c(y,input$y)

    if(length(x)>1){
      segments(x0=x[length(x)-1],x1=x[length(x)],y0=y[length(y)-1],y1=y[length(y)],lty=2)
    }
  }
  segments(x0=x[1],x1=x[length(x)],y0=y[1],y1=y[length(y)],lty=2)
  list(x=x,y=y)
}

#' Remove self intersection in polygons
#' @param poly a polygon (list with two components x and y which are equal-length numerical vectors)
#' @return A polygon without overlapping edges and new vertices corresponding to non-inner points of intersection
#' @export

polygon.clean<-function(poly){
  require(rgeos)
  require(sp)
  coords=cbind(poly$x,poly$y)
  coords=rbind(coords,coords[1,])
  s = SpatialLines(list(Lines(list(Line(coords)),ID=1)))
  s_outer = gUnaryUnion(gPolygonize(gNode(s)))
  x=s_outer@polygons[[1]]@Polygons[[1]]@coords[,"x"]
  y=s_outer@polygons[[1]]@Polygons[[1]]@coords[,"y"]
  return(list(x=x[-length(x)],y=y[-length(y)]))
}

#' Makes a color transparent
#' @export
#' @param colors A vector of colors as in `?col2rgb`
#' @param alpha transparency value (0=fully transparent, 255=fully opaque)
#'
## credit http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color

makeTransparent = function(colors, alpha=255) {
  sapply(colors,function(col){
    col=col2rgb(col)
    rgb(red=col[1,], green=col[2,], blue=col[3,], alpha=alpha, maxColorValue=255)
  })
}


#' Allows to bind the output of load() to a variable
#' @export
#' @param what character, filename passed to load

en.load<-function(what){
    tmp=new.env()
    load(what,envir=tmp)
    if(length(tmp)==1){
        return(get(ls(envir=tmp),envir=tmp))
    } else {
        return(as.list(tmp))
    }
}

#' FlowSOM wrapper to call on matrices
#' @export

en.FlowSOM<-function(mat,scale=FALSE,scaled.center=TRUE,scaled.scale=TRUE,nClus=NULL,maxMeta,...){
    fsom=list(data=mat)
    if(scale){
        fsom$data=scale(fsom$data,,scaled.center,scaled.scale)
        fsom$scaled.center <- attr(fsom$data, "scaled:center")
        attr(fsom$data, "scaled:center") <- NULL
        fsom$scaled.scale <- attr(fsom$data, "scaled:scale")
        attr(fsom$data, "scaled:scale") <- NULL
    }
    fsom=BuildSOM(fsom)
    fsom=BuildMST(fsom)
    
    if (is.null(nClus)) {
        cl <- MetaClustering(fsom$map$codes, 
            "metaClustering_consensus", maxMeta)
    }
    else {
        cl <- metaClustering_consensus(fsom$map$codes,nClus)
    }
    list(fsom,cl,cl[fsom$map$mapping[,1]])
}

#' writeClipboard for Linux
#' @export
#' @author Credit to Freecube on http://stackoverflow.com/questions/10959521/how-to-write-to-clipboard-on-ubuntu-linux-in-r
writeClipboard <- function(x, sep="\t", row.names=FALSE, col.names=TRUE){
     con <- pipe("xclip -selection clipboard -i", open="w")
     write.table(x, con, sep=sep, row.names=row.names, col.names=col.names)
     close(con)
}

#' Map a vector to colors (for categorical variables)
#' @export
toColors_discrete <- function(vector,palette=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")){
    n=sort(unique(vector))
    setNames(colorRampPalette(palette)(length(n)),n)[as.character(vector)]
}

#' Map a vector to colors (for continuous variables); old version
#' @export
toColors_continuous_fromRank<-function(vector,palette=rev(c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695"))){
    M=min(200,length(unique(vector)))
    levels=as.numeric(cut(vector,breaks=M,include.lowest=TRUE))
    colorRampPalette(palette)(M)[levels]
}

#' Map a vector to colors (for continuous variables); new version
#' @export
toColors_continuous<-function(vector,palette=rev(c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695"))){
    M=min(200,length(unique(vector)))
    levels=as.numeric(cut(vector,breaks=seq(min(vector),max(vector),length.out=M),include.lowest=TRUE))
    colorRampPalette(palette)(M)[levels]
}

#' Plot a 2D heatmap of frequencies using rectangular bins (a la FlowJo)
#' @export
freqplot=function(x,y,breaks=200,na.rm=TRUE,palette = rev(c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695")), add_white=TRUE,...){

    w=is.na(x)|is.na(y)|is.nan(x)|is.nan(y)|is.infinite(x)|is.infinite(y)
    if(any(w)){
        if(na.rm){
            x=x[!w]
            y=y[!w]
        } else {
            stop("NA values found and na.rm is FALSE")
        }
    }
    
    w.x=length(unique(x))>1
    w.y=length(unique(y))>1
    
    if(w.x){
        breaks.x=seq(min(x),max(x),length.out=breaks)
        labels.x=breaks.x[-length(breaks.x)]
        X=cut(x,breaks=breaks.x,labels=labels.x,include.lowest=TRUE)
    } else {
        X=x
    }
    if(w.y){
        breaks.y=seq(min(y),max(y),length.out=breaks)
        labels.y=breaks.y[-length(breaks.y)]
        Y=cut(y,breaks=breaks.y,labels=labels.y,include.lowest=TRUE)
    } else {
        Y=y
    }

    tab=log10(1+table(X,Y))

    if(length(x)<1|length(y)<1){
        plot.new()
        return(tab)
    }
    
    if(w.x&w.y){
        if(add_white){
            null_color = "white"
        } else {
            null_color = NULL
        }
        image(tab,col=c(null_color,colorRampPalette(palette)(100)),x=breaks.x,y=breaks.y,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",...)
        ticks=seq(0,1,by=0.25)
        axis(side=1,at=quantile(breaks.x,ticks),labels=signif(quantile(breaks.x,ticks),2),line=0.5)
        axis(side=2,at=quantile(breaks.y,ticks),labels=signif(quantile(breaks.y,ticks),2),line=0.5)
    } else {
        if(!w.x){
            X=runif(length(x))
        } else {
            X=x
        }
        if(!w.y){
            Y=runif(length(y))
        } else {
            Y=y
        }
        freqplot(X,Y,breaks=breaks,na.rm=na.rm,...)
    }
    tab
}

#' Convert line to user coordinates
#' @export
#' @author Credit to jbaums https://stackoverflow.com/questions/30765866/get-margin-line-locations-in-log-space/30835971#30835971
line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
  y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
  switch(side,
         `1` = grconvertY(-line * y_off, 'npc', 'user'),
         `2` = grconvertX(-line * x_off, 'npc', 'user'),
         `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
         `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
         stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

#' Bayes classifier - training
#' @export
prior_free_bayes=function(data,classes){
    data=split(as.data.frame(data),classes)
    lapply(data,function(x){
        cbind(colMeans(x),apply(x,2,sd))
    })
}

#' Bayes classifier - prediction
#' @export
#' @note To incorporate a class-frequency prior, multiply each score by the class frequency on the training data
prior_free_bayes_call=function(data,model){
    do.call(cbind,lapply(model,function(x){
        rowSums(do.call(cbind,lapply(1:nrow(x),function(y){
            dnorm(data[,y],x[y,1],x[y,2],log=TRUE)
        })))
    }))
}

#' Auto-logicle adapted from flowCore
#' @export
eb.autolgcl=function(x){
    x
    t=max(x)
    m=4.5
    q=0.05
    r=.Machine$double.eps + quantile(x,q)
    w=max((m-log10(t/abs(r)))/2,0.1)
    a=0
    logicleTransform(w=w,t=t,m=m,a=a) ##Just use summary() to retrive the parameters
}

#' 2D to color map
#' @export
toColors_2d=function(mat,rank=TRUE){
    if(ncol(mat)!=2){
        stop("Input a two-columns matrix")
    }
    if(rank){
        mat=apply(mat,2,function(x){order(order(x))})
    }
    mat=apply(mat,2,function(x){
        x=x-min(x)
        x=x/max(x)
    })
    o1=mat[,1]
    o2=mat[,2]
    rgb(1-o1,((o1+o2)-min(o1+o2))/(max(o1+o2)-min(o1+o2)),1-o2,maxColorValue=1)
}

#' Fast splitting of matrix to list (avoids conversion do data.frame)
#' @export
split_matrix=function(mat,vector,byrow=TRUE){
    if(byrow&nrow(mat)!=length(vector)){
        stop("if byrow=TRUE, vector's length should have length nrow(mat)")
    } else if(!byrow&ncol(mat)!=length(vector)) {
        !byrow&ncol(mat)!=length(vector)
        stop("if byrow=FALSE, vector's length should have length ncol(mat)")
    }
    
    if(byrow){
        levels=split(1:nrow(mat),vector)
        res=lapply(levels,function(x){
            mat[x,,drop=FALSE]
        })
    } else {
        levels=split(1:ncol(mat),vector)
        res=lapply(levels,function(x){
            mat[,x,drop=FALSE]
        })
    }
    res
}

#' faust:::tsGates wrapper
#' @export
tsGates=function(xVec, modePrior = 0L, maxElements = 10^4){
    require(faust)
    if(length(xVec) > maxElements){
        xVec = xVec[sample(1:length(xVec), maxElements)]
    }
    xVec = sort(xVec)
    faust:::tsGates(xVec, modePrior)
}
