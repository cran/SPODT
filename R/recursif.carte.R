recursif.carte <-
function(noeud, map.col)
{
    color<-colors()[c(116:151, 361:489, 491:657)]#[c(116:151, 361:489, 491:657)]  

    if (class(noeud) != "f.spodt")
    {
        if (class(noeud) == "sp.spodt")
        {
            X1 <- noeud@int[1]
            X2 <- noeud@int[2]
            Y1 <- noeud@coeff[1] * X1 + noeud@coeff[2]
            Y2 <- noeud@coeff[1] * X2 + noeud@coeff[2]

            segments(X1, Y1, X2, Y2, col="gray75")    # col="gray75"  pas lwd
        }
        recursif.carte(noeud@fg, map.col)
        recursif.carte(noeud@fd, map.col)
    }
    else
    {
        if (map.col)
        {
            text(noeud@G[1], noeud@G[2], noeud@id, cex=0.7, col=color[(noeud@id*7)%%332+1])
        }
        else
        {
            text(noeud@G[1], noeud@G[2], noeud@id, cex=0.7)
        }
    }
}
