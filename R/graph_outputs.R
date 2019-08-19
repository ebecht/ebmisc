                                        #Outputs, for each flowframe, for each channel, the range it should use. (so its a named list containing matrices with columns = channels and rows = min, max)
                                        #Then we ca use this in our fs.plot loop to look at the required scale.

#' Computes data ranges for each channel of each flowFrame of a flowSet
#'
#' @param fs A flowSet
#' @param global_across_channels Boolean
#' @param na.rm Boolean : wether to ignore missing values or change output to NA
#' @details See \code{\link{color_biplots_by_channels}} for details of how the boolean arguments affect the output
#' @export

data_ranges=function(fs,global_across_samples,global_across_channels,na.rm=T){
    if(global_across_samples){
        global.xp=do.call(rbind,fsApply(fs,exprs,simplify=F))
        if(global_across_channels){
            return(
                setNames(
                    rep(
                        list(
                            matrix(
                                ncol=ncol(fs[[1]]),
                                nrow=2,
                                data=rep(range(global.xp,na.rm=na.rm),ncol(fs[[1]])),
                                dimnames=list(c("min","max"),colnames(fs[[1]]))
                            )
                        ),
                        length(fs)
                    ),
                    fs@phenoData$name
                )
            )
        }
        if(!global_across_channels){
            res=apply(global.xp,2,range,na.rm=na.rm)
            rownames(res)=c("min","max")
            return(
                setNames(
                    rep(
                        list(res),
                        length(fs)
                    ),
                    fs@phenoData$name
                )
            )
        }
    }
    if(!global_across_samples){
        if(global_across_channels){
            return(
                fsApply(
                    fs,
                    function(ff){
                        matrix(
                            ncol=ncol(ff),
                            nrow=2,
                            data=rep(range(exprs(ff),na.rm=na.rm),ncol(ff)),
                            dimnames=list(c("min","max"),colnames(ff))
                        )
                    },
                    simplify=F
                )
            )
        }
        if(!global_across_channels){
            return(
                fsApply(fs,function(ff){
                    res=apply(exprs(ff),2,range,na.rm=na.rm)
                    rownames(res)=c("min","max")
                    res
                },
                simplify=F
                )
            )
        }
    }
}

#' Colors points of a biplot (2d-tSNE, 2d-PCA...) by the intensity of channels for each flowframe in the flowset
#'
#' @param matrix A matrix
#' @param file_name String of the location of the output file (should end in .pdf)
#' @param x_axis The column name of the x_axis. Usually the same name you used previously in fs.append
#' @param y_axis The column name of the y_axis. Usually the same name you used previously in fs.append
#' @param global_across_channels Boolean specificying whether the color key should be calibrated across all channels or individually per channel.
#' @param palette A vector of colors that'll be passed to colorRampPalette
#' @param resolution The resolution of the files that'll be imported in the pdf. Default is 72, increase for higher resolution
#' @param data_transformation_reverse The colors will be linearly scaled across the range of the data. If the data was transformed, you may however want the labels which appear in your color-keys to reflect the raw intensities. In this case, this should be the inverse of the function you used to transform your data
#' @param raster.width Width of each embedded raster. Defaults to 480
#' @param raster.height height of each embedded raster. Defaults to 480
#' @param ... passed to plot (suggested: pch=16, cex=0.5 or less)
#' @return NULL
#' @note Since pdf files are vectorized, they can get really big if a lot of data point are plotted. This function thus used bitmap images that are stored in a temporary directory (tmpDir()) and then import them in a single pdf. If you're interested in using the bitmap images, you can fetch them in tmpDir()
#' @export

color_biplot_by_channels <- function(
                                     matrix,
                                     x_axis,
                                     y_axis,
                                     global_across_channels=T,
                                     palette=c("blue","green","red"),
                                     resolution=72,
                                     data_transformation_reverse=identity,
                                     file_name="biplot.pdf",
                                     raster.height=480,
                                     raster.width=480,
                                     ... #pass to plot for e.g. tSNE biplot
                                     )
{
    sapply(c("png","raster","grid"),function(package){
        if(!require(package,character.only=T)){
            install.packages(pkgs=package)
            library(package,character.only=T)
        }
    })

    regular_channels=setdiff(colnames(matrix),c(x_axis,y_axis))

    if(global_across_channels){
        data_range=range(matrix[,regular_channels],na.rm=T)
        data_range=matrix(rep(data_range,length(regular_channels)),ncol=length(regular_channels),nrow=2,byrow=F,dimnames=list(c("min","max"),regular_channels))
    } else {
        data_range=apply(matrix[,regular_channels],na.rm=T,2,range,na.rm=T)
        rownames(data_range)=c("min","max")
    }

    x = matrix[,x_axis]
    y = matrix[,y_axis]
    xp = matrix[,regular_channels,drop=FALSE]

    if(any(!(is.na(x)&is.na(y)))){
        rasters=sapply(regular_channels,function(pname,xp,data.range,x,y)
        {
            color.scale=unique(colorRampPalette(palette)(1000))
            n=length(color.scale)

            breaks=unique(seq(data.range["min",pname],data.range["max",pname],length.out=n+1))
            if(length(unique(breaks))>1){
                points.colors=as.character(cut(xp[,pname],breaks=breaks,labels=color.scale))
            } else {
                points.colors=rep("lightsteelblue",length(xp[,pname]))
            }
            mainplot=paste(tmpDir(),"/mainplot_",pname,".png",sep="")
            png(mainplot,res=resolution,height=raster.height*resolution/72,width=raster.width*resolution/72)
            par("bty"="l")
            plot(
                x,
                y,
                col=points.colors,
                xlab=x_axis,
                ylab=y_axis,
                main=pname,
                ...
            )
            dev.off()

            colorscale=paste(tmpDir(),"/colorscale_",pname,".png",sep="")
            png(colorscale,res=resolution,height=raster.height/2*resolution/72,width=raster.width*resolution/72)
            plot.new()
            par("mar"=c(2,1,2,1))


            xlims=par("usr")[1:2]
            x_coords=seq(xlims[1],xlims[2],length.out=n+1)
            ylims=par("usr")[3:4]

            labels=signif(data_transformation_reverse(breaks),2)
            labels=labels[round(seq(1,length(labels),length.out=5))]
            labels.x_coords=seq(x_coords[1],x_coords[length(x_coords)],length.out=5)

            rect(border=NA,ybottom=ylims[1],ytop=ylims[2],xleft=x_coords[-length(x_coords)],xright=x_coords[-1],col=color.scale)
            text(xpd=T,y=ylims[1],pos=1,labels=labels,x=labels.x_coords)
            text(xpd=T,y=ylims[2],pos=3,labels=paste(pname,"intensity"),x=mean(xlims))
            dev.off()

            return(list(main.file=mainplot,scale.file=colorscale))
        },simplify=F,xp=xp,data.range=data_range,x=x,y=y)

        file=file_name

        pdf(file)
        if(global_across_channels){
            plot.new()
            grid.raster(readPNG(rasters[[1]]$scale.file,native=T))
        }
        sapply(rasters,function(x){
            par("mar"=c(0,0,0,0))
            grid.newpage()
            label=sub(".png","",sub("mainplot_","",tail(strsplit(x$main.file,"/")[[1]],1),fixed=TRUE),fixed=TRUE)
            
            if(!global_across_channels){
                grid.raster(readPNG(x$main.file,native=T),y=0.6,height=0.8)
                grid.raster(readPNG(x$scale.file,native=T),y=0.1,height=0.2)
            }
            if(global_across_channels){
                grid.raster(readPNG(x$main.file,native=T))
            }
            grid.text(x=unit(1,"npc"),y=unit(1,"npc"),label=label,just=c(1,1),gp=gpar(col="white",cex=0.1))
            return(NULL)
        })
        dev.off()
    }
    NULL
}

#' Wrapper to image where the coordinate system is more intuitive
#'
#' @param matrixToPlot A matrix with two dimensions
#' @param labCol if specified, the text to annotate rows
#' @param labRow if specified, the text to annotate colums
#' @param na.col color for missing values
#' @param cex.text cex value for both rows and columns
#' @param cex.row cex value for rows only. overrides cex.text
#' @param cex.col cex value for columns only. overrides cex.text
#' @param ... passed to image
#' @export

imageWrapper=function(matrixToPlot,labCol=colnames(matrixToPlot),labRow=rownames(matrixToPlot),...,cex.text=1,cex.row=1,cex.col=1,na.col=par("bg")){
    nrow=nrow(matrixToPlot)
    ncol=ncol(matrixToPlot)

    ## Force evaluation now
    labCol
    labRow 
    
    matrixToPlot=t(matrixToPlot[nrow(matrixToPlot):1,,drop=F])

    NAs=which(is.na(matrixToPlot),arr.ind=T)

    image(0:nrow(matrixToPlot),0:ncol(matrixToPlot),matrixToPlot,...,xaxt="n",yaxt="n",ylab="",xlab="",bty="n")

                                        #    segments(x0=0,x1=ncol,y0=0:nrow,col="white")
                                        #    segments(x0=0:ncol,y0=0,y1=nrow,col="white")

    if(missing(cex.col)){
        cex.col=cex.text
    }
    if(missing(cex.row)){
        cex.row=cex.text
    }

    text(x=(1:ncol)-0.5,y=nrow,pos=4,srt=90,labels=labCol,xpd=T,adj=c(0,0.5),offset=0,cex=cex.col)
    text(x=0,y=(nrow:1)-0.5,pos=2,labels=labRow,xpd=T,cex=cex.row,offset=0,adj=c(1,0.5))
    
    if(length(NAs)>0){
        rect(ybottom=NAs[,"col"]-1,ytop=NAs[,"col"],xleft=NAs[,"row"]-1,xright=NAs[,"row"],col=na.col,border=F)
    }
}

#'Computes breaks for a color scale based on a matrix
#'
#' @param Mat a matrix
#' @param n How many breaks
#' @param sym Boolean. Whether breaks should be symetric (centered at 0)
#' @return Numerical vector of length n
#' @export

fastBreaks<-function(Mat,n,sym=T){
    if(sym){
        absMax=max(abs(Mat),na.rm=T)
        return(seq(-absMax,absMax,length.out=n))
    }
    if(!sym){
        seq(min(Mat,na.rm=T),max(Mat,na.rm=T),length.out=n)
    }
}

#' Colors a biplot according to a vector with discrete values
#'
#' @param matrix a two columns matrix
#' @param discrete_vector a vector of size nrow(matrix)
#' @param ... passed to plot
#' @export

color_biplot_by_discrete<-function(matrix,discrete_vector,...,bty="l",pch=16,cex=0.5,colors=NULL,alpha=100){
    levels=unique(discrete_vector)
    if(missing(colors)){
        colors=setNames(makeTransparent(c("black",rainbow(length(levels)-1)),alpha=alpha),levels)
    }
    plot(matrix,bty=bty,pch=pch,cex=cex,col=colors[as.character(discrete_vector)],...)
}

#' Computes AUC from a vector of sorted values
#' @export
aucFromRules=function(sortedLabels,levels){
    sortedLabels=na.omit(sortedLabels)
    sortedLabels=sortedLabels%in%levels

    totalTrues=sum(sortedLabels)
    totalFalses=length(sortedLabels)-totalTrues

    tpr=c(0,cumsum(!sortedLabels)/totalFalses)
    fpr=c(0,cumsum(sortedLabels)/totalTrues)

    heights=(tpr[2:length(tpr)]+tpr[1:(length(tpr)-1)])/2
    widths=fpr[2:length(fpr)]-fpr[1:(length(fpr)-1)]
    sum(heights*widths)
}
