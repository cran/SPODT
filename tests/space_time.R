#assess the SpODT algorithm analysing time covariate


require(SPODT)
require(tree)
data(dataCOV)
dataset <- dataCOV

#mapping the data
plot(dataset$x, dataset$y, cex=dataset$z)


#SpODT tuning parameter

gr<-0 #graft parameter 0.2
rtw<-0.1 #rtwo.min  0.01
parm<-10 #min.parent 10
childm<-5 #min.child 
lmx<-3 #level.max 5

# SpODT classification

    dataset_res <- spodt(z~i+x+y,~V1,data=dataset, graft=gr, min.ch=childm, min.parent=parm, level.max=lmx, rtwo.min=rtw, weight=TRUE)

#mapping the SpODT classification

    spodt.map(dataset_res, dataset$x, dataset$y,map.col=TRUE)
    points(dataset$x, dataset$y, cex=dataset$z)
   title(main = "SpODT classification")


#SpODT tree

    spodt.tree(dataset_res)
    title(main = "SpODT classification")


## SpODT classification test
#
#    dataset_tst <- test.spodt(z~i+x+y,~V1,data=dataset, dataset_res@R2,  'rnorm', c(length(dataset$i), mean(dataset$z), var(dataset$z)), min.ch=childm, nb.sim=9, weight=TRUE, graft=gr, level.max=lmx, min.parent=parm, rtwo.min=rtw)
#

#CART classification with the tree package


    dataset_cart <- tree(z ~ x + y, data=dataset, split = c("deviance"),mincut = 5, minsize = 10, mindev = 0.01)#4

    plot(dataset$x, dataset$y, cex=dataset$z, xlab="x",ylab="y")
    partition.tree(dataset_cart, ordvars=c("x","y"), add=TRUE)
    title(main = "CART classification")

