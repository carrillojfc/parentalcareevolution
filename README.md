# parentalcareevolution

install.packages("ape")
install.packages("phangorn")
install.packages("phytools")
install.packages("geiger")
install.packages("picante")
install.packages("phyloch")
install.packages("ips")
install.packages("treeio")
install.packages("ggtree")

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(picante)
library(ips)
library(treeio)
library(ggtree)
library(lmtest)

devtools::install_github("YuLab-SMU/ggtree")
library(BiocManager)
BiocManager::install("ggtree")

#### loading

sqData<-read.csv("mode.csv",row.names=1)
dim(sqData)

sqTree<-read.nexus("tree228.txt")
print(sqTree,printlen=2)
sqTree$tip.label

plotTree(sqTree,type="phylogram",lwd=1,fsize=0.8,ftype="i")

## check name matching
chk<-name.check(sqTree,sqData)
summary(chk)

## drop tips of tree that are missing from data matrix
sqTree.pruned<-drop.tip(sqTree,chk$tree_not_data)
## drop rows of matrix that are missing from tree
sqData.pruned<-sqData[!(rownames(sqData)%in%
                          chk$data_not_tree),,drop=FALSE]

## extract discrete trait
modes<-setNames(as.factor(sqData.pruned[,"mode"]),
               rownames(sqData.pruned))
head(modes)
str(modes)

## fit ER model to squamate toe data using fitDiscrete
fitER<-fitDiscrete(sqTree.pruned,modes,model="ER")
print(fitER,digits=3)
## plot fitted ER model
plot(fitER,mar=rep(0,4),signif=5)
text(1,1,"ER")

## fit SYM model
fitSYM<-fitDiscrete(sqTree.pruned,modes,model="SYM")
print(fitSYM,digits=3)
## graph fitted SYM model
plot(fitSYM,show.zeros=TRUE,mar=rep(0,4),signif=5)
text(1,1,"SYM")

## fit ARD model
fitARD<-fitDiscrete(sqTree.pruned,modes,model="ARD")
print(fitARD,digits=3)
## plot fitted model
plot(fitARD,show.zeros=TRUE,mar=rep(0,4),signif=3)
text(1,1,"ARD")

## create design matrix for bi-directional
## ordered model

ordered.model<-matrix(c(
  0,1,0,
  2,0,3,
  0,4,0),3,3,byrow=TRUE,
  dimnames=list(0:2,0:2))
ordered.model

## create design matrix for directional ordered
## model
directional.model<-matrix(c(
  0,0,0,
  1,0,0,
  0,2,0),3,3,byrow=TRUE,
  dimnames=list(0:2,0:2))
directional.model

## fit bi-directional ordered model
fitOrdered<-fitDiscrete(sqTree.pruned,modes,
                        model=ordered.model,surpressWarnings=TRUE)
print(fitOrdered,digits=3)

## fit directional (loss only) ordered model
fitDirectional<-fitDiscrete(sqTree.pruned,modes,
                            model=directional.model,surpressWarnings=TRUE)
print(fitDirectional,digits=3)

## split plot area into panels
par(mfrow=c(3,2))
## plot ordered and directional models
plot(fitOrdered,show.zeros=FALSE,signif=3,
mar=c(0.1,1.1,0.1,0.1))
#mtext("ordered",line=-2,adj=0,cex=1.5)
text(1,1,"ORD")

  plot(fitDirectional,show.zeros=TRUE,signif=3,
mar=c(0.1,1.1,0.1,0.1))
#mtext("direct",line=-2,adj=0,cex=1.5)
text(1,1,"DIR")
dev.off()
## comapring models
library(lmtest)
## likelihood-ratio test comparing ER & SYM
lrtest(fitER,fitSYM)
## likelihood-ratio test comparing ER & ARD
lrtest(fitER,fitARD)
## likelihood-ratio test comparing SYM & ARD
lrtest(fitSYM,fitARD)
## compare directional and ordered
lrtest(fitDirectional,fitOrdered)
## compare direction and ARD
lrtest(fitDirectional,fitARD)
## compare ordered and ARD
lrtest(fitOrdered,fitARD)

## accumulate AIC scores of all five models into
## a vector
aic<-setNames(c(AIC(fitER),AIC(fitDirectional),
                AIC(fitOrdered),AIC(fitSYM),AIC(fitARD)),
              c("ER","Directional","Ordered","SYM","ARD"))
aic ## lower aic best support data

aic.w(aic)
results<-round(data.frame(
  k=c(fitER$opt$k,fitDirectional$opt$k,
      fitOrdered$opt$k,fitSYM$opt$k,fitARD$opt$k),
  logL=c(logLik(fitER),logLik(fitDirectional),
         logLik(fitOrdered),logLik(fitSYM),logLik(fitARD)),
  AIC=aic,Akaike.w=as.vector(aic.w(aic))),3)
results

## ancestral recontruction

## extract feeding mode as a vector

modes


## set colors for plotting
cols<-setNames(c("lightblue","red","blue"),levels(modes))
## plot the tree & data
plotTree.datamatrix(sqTree.pruned,as.data.frame(modes),
                    colors=list(cols),header=FALSE,fsize=0.8
                    )
## add legend
legend("topright",legend=levels(modes),pch=22,
       pt.cex=1.5,pt.bg=cols,bty="n",cex=0.8)



## ancestral reconstruction - with ORD as best model

library(corHMM)

lepto.data<-data.frame(Genus_sp=names(modes),
                     cp.mode=as.numeric(modes)-1)

head(lepto.data)
print(lepto.data)

## generate 1,000 stochastic character maps in which
## the transition rate is sampled from its posterior
## distribution
mtrees<-make.simmap(sqTree.pruned,modes,model=ordered.model,
                    nsim=1000,Q="mcmc",vQ=0.01,
                    prior=list(use.empirical=TRUE),samplefreq=10)
mtrees


## set plot margins
par(mar=c(5.1,4.1,2.1,2.1))
## create a plot of the posterior density from stochastic
## mapping 
plot(d<-density(sapply(mtrees,function(x) x$Q[1,2]),
                bw=0.005),bty="n",main="",xlab="q",xlim=c(0,0.5),
     ylab="Posterior density from MCMC",las=1,
     cex.axis=0.8)
polygon(d,col=make.transparent("blue",0.25))
## add line indicating ML solution for the same parameter
abline(v=fit.marginal$solution[1,2])
text(x=fit.marginal$solution[1,2],y=max(d$y),"MLE(q)",
     pos=4)

## create a 10 x 10 grid of plot cells
par(mfrow=c(10,10))
## graph 100 stochastic map trees, sampled evenly from
## our set of 1,000
null<-sapply(mtrees[seq(10,1000,by=10)],
             plot,colors=cols,lwd=1,ftype="off")

## compute posterior probabilities at nodes
pd<-summary(mtrees)
pd

dev.off()
## create a plot showing PP at all nodes of the tree
plot(pd,colors=cols,fsize=0.9,ftype="i",lwd=1.5,
     offset=0.2,ylim=c(2,Ntip(sqTree.pruned)-2),
     cex=c(0.4,0.3))

## add a legend
legend("bottomleft",legend=levels(modes),pch=22,
       pt.cex=1.5,pt.bg=cols,bty="n",cex=0.8)



#########

##coevolution
##

datamatr<-read.csv("mode-coevo.csv",row.names=1,
                   stringsAsFactors=TRUE)
datamatr
str(datamatr)
sqTree ## complete tree 

## check name matching

chk1<-name.check(sqTree,datamatr)
summary(chk1)


## drop tips of tree that are missing from data matrix
sqTree.pruned1<-drop.tip(sqTree,chk1$tree_not_data)
plotTree(sqTree.pruned1,type="phylogram",lwd=1,fsize=0.8,ftype="i")

## drop rows of matrix that are missing from tree
sqData.pruned1<-datamatr[!(rownames(datamatr)%in%
                             chk1$data_not_tree),,drop=FALSE]
sqData.pruned1
### testing best model 

## extractmodes as a vector
rp.mode<-setNames(sqData.pruned1[,2],rownames(sqData.pruned1))
co.mode<-setNames(sqData.pruned1[,3],rownames(sqData.pruned1))
pc.mode<-setNames(sqData.pruned1[,1],rownames(sqData.pruned1))
rp.modebin <- ifelse(rp.mode == "Aquatic", 0, 1)
co.modebin <- ifelse(co.mode == "Uncover", 0, 1)
pc.modebin <- ifelse(pc.mode == "Absent", 0, 1)
rp.modebinname<-setNames(rp.modebin,rownames(sqData.pruned1))
co.modebinname<-setNames(co.modebin,rownames(sqData.pruned1))
pc.modebinname<-setNames(pc.modebin,rownames(sqData.pruned1))


### reproductive mode best model
## fit ER model
fitERrp<-fitMk(sqTree.pruned1,rp.mode,model="ER")
fitERrp
## fit ARD model
fitARDrp<-fitMk(sqTree.pruned1,rp.mode,model="ARD")
fitARDrp

## fit bite->suction model
fit01rp<-fitMk(sqTree.pruned1,rp.mode,
             model=matrix(c(0,1,0,0),2,2,byrow=TRUE))
## fit suction->bite model
fit10rp<-fitMk(sqTree.pruned1,rp.mode,
             model=matrix(c(0,0,1,0),2,2,byrow=TRUE))
## extract AIC values for each model
aic<-c(AIC(fitERrp),AIC(fitARDrp),
       AIC(fit01rp),AIC(fit10rp))
## print summary table # suction->bite Gain # suction->bite loss
data.frame(model=c("ER","ARD",
                   "Gain","Loss"),
           logL=c(logLik(fitERrp),logLik(fitARDrp), #logLik(fitARDrp),
                  logLik(fit01rp),logLik(fit10rp)),
           AIC=aic,delta.AIC=aic-min(aic))

## likelihood-ratio test comparing ER & ARD
lrtest(fitERrp,fitARDrp)

###covering mode best model
## fit ER model
fitERco<-fitMk(sqTree.pruned1,co.mode,model="ER")
fitERco
## fit ARD model
fitARDco<-fitMk(sqTree.pruned1,co.mode,model="ARD")
fitARDco

## fit bite->suction model
fit01co<-fitMk(sqTree.pruned1,co.mode,
               model=matrix(c(0,1,0,0),2,2,byrow=TRUE))
## fit suction->bite model
fit10co<-fitMk(sqTree.pruned1,co.mode,
               model=matrix(c(0,0,1,0),2,2,byrow=TRUE))
## extract AIC values for each model
aic<-c(AIC(fitERco),AIC(fitARDco),
       AIC(fit01co),AIC(fit10co))
## print summary table 
data.frame(model=c("ER","ARD",
                   "Gain","Loss"),
           logL=c(logLik(fitERco),logLik(fitARDco), #logLik(fitARDco),
                  logLik(fit01co),logLik(fit10co)),
           AIC=aic,delta.AIC=aic-min(aic))

## likelihood-ratio test comparing ER & ARD
lrtest(fitERco,fitARDco)

## parental care mode best model
## fit ER model
fitERpc<-fitMk(sqTree.pruned1,pc.mode,model="ER")
fitERpc
## fit ARD model
fitARDpc<-fitMk(sqTree.pruned1,pc.mode,model="ARD")
fitARDpc

## directional
fit01pc<-fitMk(sqTree.pruned1,pc.mode,
               model=matrix(c(0,1,0,0),2,2,byrow=TRUE))
## fit suction->bite model
fit10pc<-fitMk(sqTree.pruned1,pc.mode,
               model=matrix(c(0,0,1,0),2,2,byrow=TRUE))
## extract AIC values for each model
aic<-c(AIC(fitERpc),AIC(fitARDpc),
       AIC(fit01pc),AIC(fit10pc))
## print summary table 
data.frame(model=c("ER","ARD",
                   "Gain","Loss"),
           logL=c(logLik(fitERco),logLik(fitARDco), #logLik(fitARDco),
                  logLik(fit01co),logLik(fit10co)),
           AIC=aic,delta.AIC=aic-min(aic))

## likelihood-ratio test comparing ER & ARD
lrtest(fitERpc,fitARDpc)

#### setting names
leptos.tree<-sqTree.pruned1
print(leptos.tree,printlen=3)
leptos.data<-sqData.pruned1
head(leptos.data)
##
#### plot with rev level of nest

leptos.data$nest <- factor(leptos.data$nest, levels = rev(levels(leptos.data$nest)))

object <- plotTree.datamatrix(
  leptos.tree,
  leptos.data,
  fsize = 0.9,
  yexp = 1,
  header = FALSE,
  xexp = 1.45,
  palettes = c("Blues", "YlOrRd", "PuBuGn")
)



## plot the tree with adjacent data matrix
#dev.off()
object<-plotTree.datamatrix(leptos.tree,leptos.data,
                            fsize=0.9,yexp=1,header=FALSE,xexp=1.45,
                            palettes=c("Blues","YlOrRd","PuBuGn"))



## Leg 1: Parental care
leg1 <- legend(
  x = "topright",
  legend = names(object$colors$pcare),
  cex = 0.9, pch = 22, pt.bg = object$colors$pcare,
  pt.cex = 1.5, bty = "n", title = "Parental care"
)

## Leg 2: Reproductive mode (debajo de la anterior)
leg2 <- legend(
  x = leg1$rect$left-2,
  y = leg1$rect$top - leg1$rect$h - 0.5,  # ajusta el espacio entre leyendas
  legend = names(object$colors$repmode),
  cex = 0.9, pch = 22, pt.bg = object$colors$repmode,
  pt.cex = 1.5, bty = "n", title = "Reproductive mode"
)

## Leg 3: Nest covering (debajo de la anterior)
leg3 <- legend(
  x = leg1$rect$left,
  y = leg2$rect$top - leg2$rect$h - 0.5,
  legend = names(object$colors$nest),
  cex = 0.9, pch = 22, pt.bg = object$colors$nest,
  pt.cex = 1.5, bty = "n", title = "Nest covering"
)


## add a legend for trait 1
leg<-legend(x="topright",names(object$colors$pcare),
            cex=0.9,pch=22,pt.bg=object$colors$pcare,
            pt.cex=1.5,bty="n",title="Parenal care")
## add a second legend for trait 2
leg<-legend(x=leg$rect$left-1.5,y=leg$rect$top-leg$rect$h,
            names(object$colors$repmode),cex=0.9,
            pch=22,pt.bg=object$colors$repmode,pt.cex=1.5,
            bty="n",title="Reproductive mode")

## add a 3 legend for trait 3
leg<-legend(x=leg$rect$left,y=leg$rect$top-leg$rect$h,
            names(object$colors$nest),cex=0.9,
            pch=22,pt.bg=object$colors$nest,pt.cex=1.5,
            bty="n",title="Nest covering")


## parental care vs reproductive modes

paternal_care<-setNames(leptos.data[,1],
                        rownames(leptos.data))

rep_mode<-setNames(leptos.data[,2],
                   rownames(leptos.data))
cover_mode<-setNames(leptos.data[,3],
                     rownames(leptos.data))


# repmode

parentalCare.fit<-fitPagel(leptos.tree,paternal_care,
                           rep_mode)
print(parentalCare.fit)


plot(parentalCare.fit,signif=4,cex.main=1,
     cex.sub=0.8,cex.traits=0.7,cex.rates=0.7,
     lwd=1)

#covering

parentalCare.fit1<-fitPagel(leptos.tree,paternal_care,
                            cover_mode)
print(parentalCare.fit1)

plot(parentalCare.fit1,signif=4,cex.main=1,
     cex.sub=0.8,cex.traits=0.7,cex.rates=0.7,
     lwd=1)
dev.off()

#### philo signal #####
## bimodal
#### pruning

sqData3<-read.csv("modebinariop01.csv",row.names=1,
                  stringsAsFactors=TRUE)
dim(sqData3)

## check name matching
chk3<-name.check(sqTree,sqData3)
summary(chk3)

## drop tips of tree that are missing from data matrix
sqTree.pruned3<-drop.tip(sqTree,chk3$tree_not_data)
plotTree(sqTree.pruned3,type="phylogram",lwd=1,fsize=0.8,ftype="i")


## drop rows of matrix that are missing from tree
sqData.pruned3<-sqData3[!(rownames(sqData3)%in%
                            chk3$data_not_tree),,drop=FALSE]
str(sqData.pruned3)

feed.mode3<-setNames(sqData.pruned3[,1],rownames(sqData.pruned3))


phylosig(sqTree.pruned3,feed.mode3,method="lambda",test=TRUE)


### multimodal
#### pruning

sqData
sqTree
sqTree.pruned
sqData.pruned
str(sqData.pruned)
str(sqTree.pruned)


mode3<-setNames(sqData.pruned[,1],rownames(sqData.pruned))


phylosig(sqTree.pruned,mode3,method="lambda",test=TRUE)
