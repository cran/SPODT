spodt <-
function(formula, cofactors ,data,
                  weight=FALSE, graft=0,
                  level.max=5, min.parent=10, min.child=5, rtwo.min=0.001)
{   
#extracting formula
     Call <- match.call()
        indx <- match(c("formula", "cofactors", "data"),
                      names(Call), nomatch = 0L)
        if (indx[1] == 0L) stop("a 'formula' with the location ID and coordinates is required")
        if (indx[2] == 0L) stop("a 'formula' with the cofactors is required\n for single spatial analysis (with no cofactor) the right hand side should be ~1")
        dataset_coords <- model.frame(Call$formula, data=data)
        dataset_fact <- model.frame(Call$cofactors, data=data)

    dataset<-cbind(dataset_coords[,2], dataset_coords[,3],dataset_coords[,4],dataset_coords[,1], dataset_fact)
    colnames(dataset)[1:4] <- c("loc", "x", "y", "z")
    ponderer<-weight
    greffer<-graft
    nv.max<-level.max
    min.pere<-min.parent
    min.fils<-min.child
    var.exp.min<-rtwo.min
    n <- nrow(dataset)
    rownames(dataset) <- 1:n
    ###############"
    #extracting numeric and factor cofactors
    mv<-data.frame(dataset[,-c(1:4)])
    colnames(mv)<-colnames(dataset)[-c(1:4)]

 L <- dim(mv)[2]
 ind_num <- 0
 ind_fact <- 0
 
 if(L ==0){
    vqt<-NULL
    vql<-NULL
 }
 if(L != 0){
    for (i in 1:L){
        if (class(mv[,i])=="numeric" | class(mv[,i])=="integer"){
            ind_num <- c(ind_num, i)
        }
        else if(class(mv[,i])=="factor"){
            ind_fact <- c(ind_fact, i)
        }
        else{
            print(paste("class of ", colnames(mv)[i], " must be numeric or factor", sep=""))
        }
    }
    ind_num  <- ind_num[-1]
    ind_fact <- ind_fact[-1]
    vqt <- colnames(mv)[ind_num]
    vql <- colnames(mv)[ind_fact]
 }


#######

    data.sp <- prep.data.sp(dataset)

    perm <- calculer.pentes(dataset, data.sp)
    bord <- prep.bord(dataset)
    
    if (length(vqt) > 1)
    {
        ordre.vqt <- apply(dataset[,vqt], MARGIN=2, order)
    }
    else
    {
        if (length(vqt) == 1)
        {
            ordre.vqt <- order(dataset[,vqt])
        }
        else
        {
            ordre.vqt <- NULL
        }
    }


    res <- creer.noeud(dataset, data.sp, perm, vql, vqt, ordre.vqt, bord,
                       n, 1, 1,
                       ponderer, greffer,
                       nv.max, min.pere, min.fils, var.exp.min)

    if (class(res$noeud) == "f.spodt")
    {
        warning("Root is a leaf")
    }

    arbre <- new("spodt")
    arbre@racine <- res$noeud
    arbre@partition <- res$part

    res <- classe.adj(res$bord)

    arbre@adj <- res$adj    


    if (greffer)
    {
        if (is.numeric(greffer))
        {
            var.exp.min <- greffer
        }
        
        bord <- transformer.bord(res$bord)

        
        res <- realiser.greffe(dataset, bord, arbre@partition, ponderer, var.exp.min, arbre@adj)

        if (!res$grf)
        {
            arbre@R2 <- R2.global(dataset$z, arbre@partition)
            return(arbre)
        }
                                                             
        arbre@adj <- unlist(res$adj)
        arbre@R2 <- R2.global(dataset$z, arbre@partition)

        cl.grf <- matrix(as.integer(res$cl.grf ), ncol=3)
        dimnames(cl.grf) <- list(1:nrow(cl.grf), c("c1", "c2", "union(c1,c2)"))
        arbre@cl.grf <- cl.grf

        sgmts.grf <- matrix(res$sgmts.grf, ncol=4)
        dimnames(sgmts.grf) <- list(1:nrow(sgmts.grf), c("X1", "Y1", "X2", "Y2"))
        arbre@sgmts.grf <- sgmts.grf

        arbre@partition <- unlist(res$part)
        
        arbre@R2 <- R2.global(dataset$z, arbre@partition)
    }

    return(arbre)
}
