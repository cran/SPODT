#assess the SpODT algorithm to detect a rotated square shape situation: Clustered data with higher values inside the rotated square shape


require(SPODT)
require(tree)

data(dataSQUARE0) #no cluster
data(dataSQUARE0_5) #cluster with moderate values
data(dataSQUARE1_5) #cluster with high values
data(dataSQUARE2) #cluster with very high values

beta<-c(0,0.5,1.5,2)

dataset <- vector('list', 4)
dataset[[1]]$dat<-dataSQUARE0
dataset[[2]]$dat<-dataSQUARE0_5
dataset[[3]]$dat<-dataSQUARE1_5
dataset[[4]]$dat<-dataSQUARE2

#mapping the data)

par(mfrow=c(2,2), mar=c(3,3,1,0.1), mgp=c(2,0.7,0))
for (i in 1:4) {
    plot(dataset[[i]]$dat$x, dataset[[i]]$dat$y, cex=dataset[[i]]$dat$z)
    title(main = bquote(paste(beta,"=",.(beta[i]))), line=0.5)
}

#SpODT tuning parameter

gr<-0.2 #graft parameter
rtw<-0.001 #rtwo.min
parm<-10 #min.parent
childm<-5 #min.child
lmx<-5 #level.max
n<-300 #number of lines

# SpODT classification
for (i in 1:4) {
    dataset[[i]]$res <- spodt(z~i+x+y,~1,data=dataset[[i]]$dat, graft=gr, min.ch=childm, min.parent=parm, level.max=lmx, rtwo.min=rtw, weight=TRUE)
    cat(i, '\n')
}

#mapping the SpODT classification
par(mfrow=c(2,2), mar=c(3,3,1,0.1), mgp=c(2,0.7,0))
for (i in 1:4) {
    spodt.map(dataset[[i]]$res, dataset[[i]]$dat$x, dataset[[i]]$dat$y,map.col=TRUE)
    points(dataset[[i]]$dat$x, dataset[[i]]$dat$y, cex=dataset[[i]]$dat$z)
    title(main = bquote(paste(beta,"=",.(beta[i]))), line=0.5)
}

#SpODT tree
for (i in 2:4) {
    spodt.tree(dataset[[i]]$res)
    title(main = bquote(paste(beta,"=",.(beta[i]))), line=0.5)
}

## SpODT classification test
#for (i in 1:4) {
#    dataset[[i]]$tst <- test.spodt(z~i+x+y,~1,data=dataset[[i]]$dat,dataset[[i]]$res@R2,  'rnorm', c(n, mean(dataset[[i]]$dat$z), var(dataset[[i]]$dat$z)), min.ch=childm, nb.sim=9, weight=TRUE, graft=gr, level.max=lmx, min.parent=parm, rtwo.min=rtw)
#    cat(i, '\n')
#}
#
#par(mfrow=c(2,2), mar=c(3,3,1,0.1), mgp=c(2,0.7,0))
#for (i in 1:4) {
#    hist(dataset[[i]]$tst$R2.sim, freq=TRUE, xlim=0:1, xlab='R2', main=bquote(paste(beta,"=",.(beta[i]))))
#    abline(v=dataset[[i]]$res@R2, col=2)
#}
#
#
#CART classification with the tree package

for (i in 1:4) {
    dataset[[i]]$cart <- tree(z ~ x + y, data=dataset[[i]]$dat, split = c("deviance"),mincut = 5, minsize = 10, mindev = 0.01)#4
}
par(mfrow=c(2,2), mar=c(3,3,1,0.1), mgp=c(2,0.7,0))
for (i in 1:4) {
    plot(dataset[[i]]$dat$x, dataset[[i]]$dat$y, cex=dataset[[i]]$dat$z, xlab="x",ylab="y")
    partition.tree(dataset[[i]]$cart, ordvars=c("x","y"), add=TRUE)
    title(main = bquote(paste(beta,"=",.(beta[i]))), line=0.5)
}
