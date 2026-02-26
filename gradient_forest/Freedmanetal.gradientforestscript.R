
library(gradientForest)
library(data.table)
library(tidyverse)

Gchimp=fread("newestGchimp.csv",header=T) #SNPdata with NA's
Gchimp=as.data.frame(Gchimp)
Gchimp=Gchimp[,colSums(is.na(Gchimp))==0] #SNP data without NA's
Gchimp=Gchimp[-c(2:3,7,11),]
Echimp=read.csv("nechimp2nll.csv")
Echimp=Echimp[-c(2:3,7,11),]

preds <- colnames(Echimp)
specs <- colnames(Gchimp)

nSites <- dim(Gchimp)[1]
nSpecs <- dim(Gchimp)[2]
# set depth of conditional permutation
lev <- floor(log2(nSites*0.368/2))
lev

chimpforest=gradientForest(cbind(Echimp,Gchimp), predictor.vars=preds, response.vars=specs,ntree=100, transform = NULL, compact=F,nbin=100, maxLevel=1,trace=T)

load(file="allnllchimpforest.rda")



for (i in 1:200){

  chimpforestR=gradientForest(cbind(Echimp,Gchimp), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=F,nbin=100, maxLevel=1,trace=T)
  
  realtotal=chimpforestR$species.pos.rsq
  realaverage=sum(chimpforestR$result)/realtotal
  
  write.table(realtotal,file="realtotal.csv",sep=",",col.names=F,row.names=F,append=TRUE)
  write.table(realaverage,file="realaverage.csv",sep=",",col.names=F,row.names=F,append=TRUE)
}

setwd("~/Desktop/R/Chimps/RealSNPResults")

for (i in 1:200){
  
  chimpforestR=gradientForest(cbind(Echimp,Gchimp), predictor.vars=preds, response.vars=specs, ntree=10, transform = NULL, compact=F,nbin=100, maxLevel=1,trace=T)
  
  resulttotal=chimpforestR$result
  
  write.table(resulttotal,file=paste0("resulttotal",i,".csv"),sep=",",col.names=F,row.names=T,append=FALSE)
}
#crosstry=crossing(Echimp,Gchimp)
#realreal=cbind(Echimp,Gchimp)
#realreal=as_tibble(realreal)
#crosstry2=anti_join(crosstry,realreal)

setwd("~/Desktop/R/Chimps/RandomSNPResults")

for (i in 1:200){
  
  predsR=Echimp[sample(nrow(Echimp)),]
  
  chimpforestR=gradientForest(cbind(predsR,Gchimp), predictor.vars=colnames(predsR), response.vars=colnames(Gchimp), ntree=10, transform = NULL, compact=F,nbin=100, maxLevel=1,trace=T)
  
  resulttotal=chimpforestR$result
  
  write.table(resulttotal,file=paste0("randresulttotal",i,".csv"),sep=",",col.names=F,row.names=T,append=FALSE)
}

for (i in 1:200){
  
  predsR=Echimp[sample(nrow(Echimp)),]
  
  #predsR=crosstry2[sample(nrow(crosstry2),nrow(Echimp)),]
  #predsR=as.data.frame(predsR)

  chimpforestR=gradientForest(cbind(predsR,Gchimp), predictor.vars=colnames(predsR), response.vars=colnames(Gchimp), ntree=10, transform = NULL, compact=F,nbin=100, maxLevel=1,trace=T)
  
  randtotal=chimpforestR$species.pos.rsq
  randaverage=sum(chimpforestR$result)/randtotal
  
  write.table(randtotal,file="randtotal.csv",sep=",",col.names=F,row.names=F,append=TRUE)
  write.table(randaverage,file="randaverage.csv",sep=",",col.names=F,row.names=F,append=TRUE)
}

randa=read.csv("randaverage.csv",sep=",",header=FALSE)
randt=read.csv("randtotal.csv",sep=",",header=FALSE)

hist(randa$V1,xlim=c(0,1.0),main="Average R-squared of Random Gradient Forests",xlab="Average R-squared of SNPs")
abline(v=realaverage,col="red")
hist(randt$V1,xlim=c(0,1000),main="Total Correlated SNPs of Random Gradient Forests",xlab="Total SNPs of R-squared > 1")
abline(v=realtotal,col="red")


chimpgrid=read.table("chimp100k.csv",header=T,sep=",")
#chimpforest=chimpforesta

#predictoroverallimportance #allnll
plot(chimpforest,plot.type="O",	col = c("darkgreen","darkgreen","indianred4","lightskyblue3","lightskyblue3","darkgreen","darkgreen","lightskyblue3","darkgreen","indianred4","lightskyblue3","grey","grey","indianred4","grey","lightskyblue3"),
     par.args=list(mfrow = c(1,1), mar = c(2,8,1,1)+0.3, cex = 2))
legend("bottomright", inset = 0.08, legend=c("Temperature", "Precipitation", "Vegetation","Other"),
       col=c("indianred4", "lightskyblue3","darkgreen","grey"), pch=15, cex=1,pt.cex = 3)

#allwll
plot(chimpforest,plot.type="O",	col = c("darkgreen","darkgreen","indianred4","lightskyblue3","lightskyblue3","darkgreen","darkgreen","lightskyblue3","darkgreen","indianred4","lightskyblue3","grey","grey","grey","grey","indianred4","grey","lightskyblue3"),
     par.args=list(mfrow = c(1,1), mar = c(2,8,1,1)+0.3, cex = 2))
legend("bottomright", inset = 0.08, legend=c("Temperature", "Precipitation", "Vegetation","Other"),
       col=c("indianred4", "lightskyblue3","darkgreen","grey"), pch=15, cex=1,pt.cex = 3)

#elliotiwll
plot(chimpforest,plot.type="O",	col = c("darkgreen","grey","darkgreen","indianred4","darkgreen","lightskyblue3","indianred4","indianred4","lightskyblue3","darkgreen","grey","grey","lightskyblue3","grey","darkgreen","lightskyblue3","grey","lightskyblue3"),
     par.args=list(mfrow = c(1,1), mar = c(2,8,1,1)+0.3, cex = 2))
legend("bottomright", inset = 0.08, legend=c("Temperature", "Precipitation", "Vegetation","Other"),
       col=c("indianred4", "lightskyblue3","darkgreen","grey"), pch=15, cex=1,pt.cex = 3)

#elliotinll
plot(chimpforest,plot.type="O",	col = c("darkgreen","darkgreen","indianred4","lightskyblue3","indianred4","darkgreen","indianred4","lightskyblue3","darkgreen","grey","grey","lightskyblue3","grey","darkgreen","lightskyblue3","lightskyblue3"),
     par.args=list(mfrow = c(1,1), mar = c(2,8,1,1)+0.3, cex = 2))
legend("bottomright", inset = 0.08, legend=c("Temperature", "Precipitation", "Vegetation","Other"),
       col=c("indianred4", "lightskyblue3","darkgreen","grey"), pch=15, cex=1,pt.cex = 3)


predictors=names(importance(chimpforest))

#splitsdensityplots
plot(chimpforest, plot.type="S", imp.vars=allpredictors, leg.posn="topright", cex.legend=0.4, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1)))

#speciescumulativeplot
plot(chimpforest, plot.type="Cumulative.Importance", imp.vars=allpredictors, show.overall=T, legend=T,common.scale=T,leg.posn="topleft", leg.nspecies=5, cex.lab=0.7, cex.legend=0.4, cex.axis=0.6, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1),omi=c(0,0.3,0,0)))

#predictorcumulative
plot(chimpforest, plot.type="C", imp.vars=allpredictors, show.species=F, common.scale=T, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(2.5,1.0,0.1,0.5), omi=c(0,0.3,0,0)))

#R2
plot(chimpforest, plot.type="P", show.names=F, horizontal=F, cex.axis=1, cex.labels=0.7, line=2.5)

#transform grid and environmental predictors
predictors <- names(importance(chimpforest)[1:16])
tchimpgrid=cbind(chimpgrid[,c("long","lat")], predict(chimpforest,chimpgrid[,predictors]))
Trns_site <- predict(chimpforest)

#pcs
PCs=prcomp(tchimpgrid[,3:5])
sgn <- sign(PCs$rotation["bio19",])
PCs$rotation <- sweep(PCs$rotation,2,sgn,"*")
PCs$x <- sweep(PCs$x,2,sgn,"*")
# set up a colour palette for the mapping
a1 <- PCs$x[,1]
a2 <- PCs$x[,2]
a3 <- PCs$x[,3]
r <- a1+a2
g <- -a2
b <- a3+a2-a1
r <- (r-min(r)) / (max(r)-min(r)) * 255
g <- (g-min(g)) / (max(g)-min(g)) * 255
b <- (b-min(b)) / (max(b)-min(b)) * 255

nvs <- dim(PCs$rotation)[1] # number of variables
vec <- c("bio19","srtmstd","bio3") 
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
# choose a scaling factor to plot the vectors over the grid
scal <- 40
xrng <- range(PCs$x[,1], PCs$rotation[,1]/scal)*1.1
yrng <- range(PCs$x[,2], PCs$rotation[,2]/scal)*1.1
plot((PCs$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=7, col=rgb(r,g,b, max = 255), asp=1)
# plot the other predictors with "+"
points(PCs$rotation[! vind,1:2]/scal, pch="+")  
# plot the chosen predictors as arrows
arrows(rep(0,lv), rep(0,lv), PCs$rotation[,1]/scal, PCs$rotation[,2]/scal, length = 0.01)
jit <- 0.00015
text(PCs$rotation[vec,1]/scal+jit*sign(PCs$rotation[vec,1]), PCs$rotation[vec,2]/scal+jit*sign(PCs$rotation[vec,2]), labels = vec)

# first predict the PCs for the transformed site data
PCsites <- predict(PCs,Trns_site[,predictors])
# plot all the sites as points on the biplot
points(PCsites[,1:2])
# calc & plot the weighted mean locations of species (from gf$Y)
#SpsWtd <- sweep(chimpforest$Y,2,apply(chimpforest$Y,2,min),"-")
#SpsWtdPCs <- (t(SpsWtd) %*% (PCsites[,1:2]))/colSums(SpsWtd)
#points(SpsWtdPCs, col="red", pch="+")
# interactively label some of the species, if desired
#identify(SpsWtdPCs, labels = as.character(rownames(SpsWtdPCs)), col="blue")


#plot these in geographic space

chimp.pred <- predict(chimpforest, chimpgrid[,predictors])
plot(tchimpgrid[,c("long","lat")],pch=15,cex=0.5,asp=1,col=rgb(r,g,b, max=255))
points(Echimp[,1:2])

#export map for use in QGIS
greencols=rgb(r,g,b,max=255)
greencols2=col2rgb(greencols)
greencols3=t(greencols2)
gradients=cbind(tchimpgrid[c("long","lat")],greencols3)
write.csv(gradients,file="gradients4arcgis.csv")

#stats on average runs and associations

randa=read.csv("randaverage.csv",sep=",",header=FALSE)
randt=read.csv("randtotal.csv",sep=",",header=FALSE)
reala=read.csv("realaverage.csv",sep=",",header=FALSE)
realt=read.csv("realtotal.csv",sep=",",header=FALSE)
colnames(randa)=c("randa")
colnames(randt)=c("randt")
colnames(reala)=c("reala")
colnames(realt)=c("realt")


#realaverage=0.15513504
#realtotal=588

r1 = rgb(255,0,0, max = 255, alpha = 95, names = "lt.pink")

hist(randa$randa,xlim=c(0.1,0.25),main="Average R-squared of Observed Versus Random Gradient Forests",xlab="Average R-squared of SNPs")
hist(reala$reala,xlim=c(0.1,0.25),col=r1,add=TRUE)
#abline(v=realaverage,col="red")
hist(randt$randt,xlim=c(100,1000),main="Total Correlated SNPs of Observed Versus Random Gradient Forests",xlab="Total SNPs of R-squared > 1")
hist(realt$realt,xlim=c(100,1000),col=r1,add=TRUE)
#abline(v=realtotal,col="red")

rravg=cbind(reala,randa)
rravg=as_tibble(rravg)

rrtot=cbind(realt,randt)
rrtot=as_tibble(rrtot)

ttesttot=t.test(rrtot$realt,rrtot$randt)
ttestavg=t.test(rravg$reala,rravg$randa)

ttesttot
ttestavg

save(chimpforest,file="elliotinllforest.rda")
