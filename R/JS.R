#' @export
js.div<-function(xp2d,events.code,n=40,n.x=NULL,n.y=NULL){
    if(missing(n.x)){
        n.x=n
    }
    if(missing(n.y)){
        n.y=n
    }

    w=!apply(xp2d,1,function(x)any(is.na(x)))
    xp2d=xp2d[w,]
    events.code=events.code[w]
    lab.x=as.character(cut(xp2d[,1],n.x,labels=1:n.x))
    lab.y=as.character(cut(xp2d[,2],n.y,labels=1:n.y))
    lab=paste(lab.x,lab.y,sep="-")

    tab=table(lab,events.code)

    ##Normalizing the data to avoid non-zero entries. Credit http://mathoverflow.net/questions/72668/how-to-compute-kl-divergence-when-pmf-contains-0s
    tab=apply(tab,2,function(x){
        x=(x+1)/(sum(x)+length(x))
    })

    res=matrix(nrow=ncol(tab),ncol=ncol(tab),dimnames=list(colnames(tab),colnames(tab)),data=NA)
    diag(res)=0
    for(i in 1:(ncol(tab)-1)){
        for(j in (i+1):ncol(tab)){
            x=tab[,i]
            y=tab[,j]
            m=1/2*(x+y)
            kl=(sum(x*log2(x/m))+sum(y*log2(y/m)))/2 ## log2 will make JS bounded between 0/1
            res[i,j]=kl
            res[j,i]=kl
        }
    }
    res
}
