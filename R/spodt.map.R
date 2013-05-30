spodt.map <-
function(object, x=NULL, y=NULL, map.col=FALSE, col.sgmts="black")   #ancien carte
{
    
    if (map.col)
    {
        col.sgmts <- "white"  #white
        col.fond <- "gray50" #gray50
        color<-colors()[c(116:124, 361:489, 491:657)]#rainbow(n)#colors()[c(1:151, 361:489, 491:657)]#[c(116:151, 361:489, 491:657)]    
    }
    else
    {
        col.fond <- "white"
    }    
    if (!is.null(x)  |  !is.null(y))
    {
        minX <- min(x) - (max(x) - min(x))/50
        maxX <- max(x) + (max(x) - min(x))/50
        minY <- min(y) - (max(y) - min(y))/50
        maxY <- max(y) + (max(y) - min(y))/50
        
        X <- c(minX, minX, maxX, maxX)
        Y <- c(minY, maxY, maxY, minY)
        
        plot(c(minX,maxX), c(minY,maxY), col="white", xlab="x", ylab="y", asp=1)

        rect(minX, minY, maxX, maxY, col=col.fond)
    }    

    recursif.carte(object@racine, map.col)

    if (nrow(object@sgmts.grf) != 0)
    {                                                                   
        segments(object@sgmts.grf[,1], object@sgmts.grf[,2], object@sgmts.grf[,3], object@sgmts.grf[,4], col=col.fond, lwd=3)   #lwd3     col=col.fond
    }    
    if (map.col  &  (!is.null(x)  |  !is.null(y)))
    {
        points(x,y, pch=16, cex=0.7, col=color[(object@partition*7)%%332+1])    #col=color[(object@partition*7)%%332+1])       pch=16, cex=0.7
    }
}
