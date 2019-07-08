#' A wrapper to flowFrame constructor that also updates the "desc" parameters based on fs_src
#'
#' @param matrix A numeric matrix
#' @return a flowFrame whose expression correspond to the matrix provided. If some channels are common with fs_src, the function updates the flowframe with the corresponding description (ie protein targets)
#' @details If replace is true, each element ff will have n events. If it is false, ff will have min(n,nrow(ff))
#' @export
matrix_to_flowframe <- function(matrix){
    res=flowFrame(matrix)

    res@parameters$desc=as.character(res@parameters$desc)
    res@parameters$name=as.character(res@parameters$name)

    shared=intersect(as.character(fs_src[[1]]@parameters$desc),colnames(matrix))

    tester=res@parameters$desc%in%shared

    res@parameters$name[tester]=as.character(channels.code[as.character(res@parameters$name)[tester]])

    res@parameters$desc[!tester]=res@parameters$name[!tester]

    res
}

#' Takes from input a channelsname.csv file and returns columns to select as a vector
#'
#' @param channels.namefile A csv file
#' @param  to_config TRUE if file is in ./config/ subfolder, FALSE if it's just a file path
#' @export
#' @return a character vector of column names to be selected

namefile_to_colnames<-function(channels.namefile,to_config=TRUE)
{
    if(to_config){
        channels.namefile=to_config(channels.namefile)
    }
    s=read.csv(channels.namefile,stringsAsFactors=F,na.strings=c("","NA"))
    s=subset(s,Select=="y")
    s=s$Target
    s
}

#' From a biplot let the user interactively draw polygons to create a "Gate" column in the expression data
#'
#' @param matrix A matrix
#' @param ... passed to plot
#' @param alpha Transparency value (between 0 and 255, 0 is fully transparent)
#' @export
#' @return A named vector of length nrow(matrix) and names rownames(matrix). Ungated events are set to NA

gate_from_biplot<-function(matrix,x_axis,y_axis,...,bty="l",pch=16,cex=0.5,alpha=100,sample=NULL)
{
  require("sp")
  xp=matrix[,c(x_axis,y_axis)]
  if(!is.null(sample)){
    s=sort(sample(1:nrow(xp),sample))
  } else {
    s=1:nrow(xp)
  }

  gate_updated=rep(0,nrow(xp))
  color_biplot_by_discrete(xp[s,],gate_updated[s],bty=bty,pch=pch,cex=cex,alpha=alpha,...)

  input.message=" n : new gate, c : redo last, s : stop gating. "
  cat("\nPlease use the mouse pointer to draw")
  polygons=list()
  i=0
  u="n"
  while(u!="s"){
    if(!u%in%c("n","s","c")){
      u=readline(paste("Incorrect input.",input.message,sep=""))
    }
    if(u=="n"){
      gate=gate_updated
      i=i+1
      col=setNames(makeTransparent(c("black",rainbow(i)),alpha=alpha),0:i)

      new.pol=en.locator()

      new.pol=polygon.clean(new.pol)
      polygons=c(polygons,list(new.pol))
      gate_updated=update_gate(xp,polygons[[i]],gate,i)

      color_biplot_by_discrete(xp[s,],gate_updated[s],bty=bty,pch=pch,cex=cex,colors=col,...)

    }
    if(u=="c"){
      gate_updated=gate
      color_biplot_by_discrete(xp[s,],gate_updated[s],bty=bty,pch=pch,cex=cex,colors=col,...)
      new.pol=en.locator()
      new.pol=polygon.clean(new.pol)
      polygons[[i]]=new.pol
      gate_updated=update_gate(xp,polygons[[i]],gate_updated,i)
      color_biplot_by_discrete(xp[s,],gate_updated[s],bty=bty,pch=pch,cex=cex,colors=col,...)
    }
    u=readline(paste("Input ?",input.message,"\n",sep=""))
  }

  sapply(1:length(polygons),function(i,polygons){
    poly=polygons[[i]]
    coords = cbind(poly$x, poly$y)
    coords=rbind(coords,coords[1,])
    s = SpatialLines(list(Lines(list(Line(coords)),ID=1)))
    text(coordinates(gCentroid(s))[,1],coordinates(gCentroid(s))[,2],labels=as.character(i),xpd=T,font=2)
  },polygons=polygons)
  gate=gate_updated

  gate[gate==0]=NA
  gate[apply(xp,1,function(x)any(is.na(x)))]=NA

  setNames(gate,rownames(matrix))
}

#' Updates a gate vector
#'
#' @param xp A two colums matrix
#' @param polygon.list A list of lists with two components x and y of equal lenghts and numeric values
#' @param gate_vector a vector of length nrow(xp) with integer values
#' @param value The number that will be assigned to gate_vector, corresponding to points that lie in the polygon
#' @return The updated gate_vector
#' @export

update_gate=function(xp,polygon,gate_vector=rep(0,nrow(xp)),value=1){
  gate_vector[point.in.polygon(xp[,1],xp[,2],polygon$x,polygon$y)!=0]=value
  gate_vector
}
