#assess the SpODT algorithm analysing the number of malaria episodes per child at each household, from November to December 2009, Bandiagara, Mali. Copyright: Pr Ogobara Doumbo, MRTC, Bamako, Mali. email: okd@icermali.org


require(SPODT)
require(tree)
data(dataMALARIA)
dataset <- dataMALARIA

#mapping the data

 plot(dataset$x, dataset$y, cex=5*log(0.5+dataset$z))
    title(main = "Malaria episodes, Bandiagara, Mali")

#SpODT tuning parameter

gr<-0.1 #graft parameter
rtw<-0.01 #rtwo.min
parm<-25 #min.parent
childm<-2 #min.child
lmx<-7 #level.max

# SpODT classification
system.time(
    dataset_res <- spodt(z~loc+x+y,~1,data=dataset, graft=gr, min.ch=childm, min.parent=parm, level.max=lmx, rtwo.min=rtw, weight=TRUE)
)
#mapping the SpODT classification

    spodt.map(dataset_res, dataset$x, dataset$y,map.col=TRUE)
    points(dataset$x, dataset$y, cex=5*log(0.5+dataset$z))
   title(main = "SpODT classification of Malaria episodes, Bandiagara, Mali")


#SpODT tree

    spodt.tree(dataset_res)
    title(main = "SpODT classification of Malaria episodes, Bandiagara, Mali")


## SpODT classification test
#system.time(
#    dataset_tst <- test.spodt(z~loc+x+y,~1,data=dataset, dataset_res@R2,  'rpois', c(length(dataset$loc), mean(dataset$z)), min.ch=childm, nb.sim=9, weight=TRUE, graft=gr, level.max=lmx, min.parent=parm, rtwo.min=rtw)
#)
#
#CART classification with the tree package


    dataset_cart <- tree(z ~ x + y, data=dataset, split = c("deviance"),mincut = 5, minsize = 10, mindev = 0.01)#4


    plot(dataset$x, dataset$y, cex=5*log(0.5+dataset$z), xlab="x",ylab="y")
    partition.tree(dataset_cart, ordvars=c("x","y"), add=TRUE)
    title(main = "CART classification of Malaria episodes, Bandiagara, Mali")

