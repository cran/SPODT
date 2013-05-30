test.spodt <-
function(R2.obs, data, rdist, par.rdist, nb.sim,
           qt.fact=NULL, ql.fact=NULL, weight, graft,
           level.max, min.parent, min.child, rtwo.min)
{
vqt<-qt.fact
vql<-ql.fact
  tloi<-get(rdist)
  if (rdist=="rbinom" | rdist=="runif" | rdist=="rnorm")
  {
    
    data.sim <- data.frame(replicate(nb.sim, tloi(par.rdist[1],par.rdist[2],par.rdist[3])))
  }
  else if  (rdist=="rpois") 
   {
  
    data.sim <- data.frame(replicate(nb.sim, tloi(par.rdist[1],par.rdist[2])))
   }
   else if  (rdist=="rnbinom")
   {
    data.sim <- data.frame(replicate(nb.sim, tloi(par.rdist[1],par.rdist[2],par.rdist[3],par.rdist[4])))
   }
   else
   {
   warning("wrong distribution")
    }
    R2.sim <- apply(data.sim, MARGIN=2, simuler.spodt, data=data,
                    vqt=vqt, vql=vql, weight=weight, graft=graft,
                    level.max=level.max, min.parent=min.parent, min.child=min.child, rtwo.min=rtwo.min)

   hist(R2.sim, xlim=c(0,1))
   abline(v=R2.obs, col="red")
   
   return(list(R2.sim=R2.sim, p=length(R2.sim[which(R2.sim > R2.obs)])/(length(R2.sim))))
}
