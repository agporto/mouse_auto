library(geomorph)
library(tools)
library(stringdist)
library( abind )

load_data = function(x){
  lmlist = list( )
  for ( i in 1:length(x) ) {
    data = read.csv(x[i],skip = 3, header=F,sep = ",", stringsAsFactors = F )
    lm = as.matrix(data[,2:4])
    lmlist[[ i ]] = lm
      }
  return(lmlist)
}

closest = function(string, stringVector){
  
  stringVector[amatch(string, stringVector, maxDist=Inf)]
  
}

bd = '/home/ap/Dropbox/Year_2020/For_Murat/mouse_auto'
set.seed(1234)
if ( ! dir.exists( bd ) )
  stop("set base directory to point to the data folder of the cloned repository.")

#original = Sys.glob( paste(bd,'/original_mouse/*.fcsv',sep='') )

multireg_balb = Sys.glob( paste(bd,'/predictions_multi_BALB/*.fcsv',sep='') )
multireg_aj = Sys.glob( paste(bd,'/predictions_multi_AJ/*.fcsv',sep='') )
pairreg_balb = Sys.glob( paste(bd,'/predictions_pair_balb/*.fcsv',sep='') )
pairreg_aj = Sys.glob( paste(bd,'/predictions_pair_aj/*.fcsv',sep='') )


##CHOOSE HERE WHICH DATASET TO ANALYZE
pseudo_data = load_data(pairreg_aj)

mydims = c( nrow( pseudo_data[[1]] ), ncol( pseudo_data[[1]]), length(pseudo_data) )
names = file_path_sans_ext(basename(pairreg_aj))

#load LM data
load(paste(bd,"/LMs.RData",sep="/")) #LMs manually annotated  from the individual skull.
old_names = dimnames(LM)[[3]]
pseudoLM = array( unlist( pseudo_data ), dim = mydims )
dimnames(pseudoLM)[[3]] = names
new_names=NULL
for (i in 1:length(old_names)){new_names[[i]] = closest(old_names[i],names)}
dimnames(LM)[[3]] = new_names

exclude.lm = which(!as.character(1:51) %in% rownames(LM[,,1])) 

names = new_names

GPA=gpagen(LM[,,names])
GPA.pca=plotTangentSpace(GPA$coords)

pseudoGPA=gpagen(pseudoLM[-exclude.lm,,names])
pseudoGPA.pca=plotTangentSpace(pseudoGPA$coords)
joint=abind(LM,pseudoLM[-exclude.lm,,names])
joint.GPA=gpagen(joint)

par(mfrow=c(1,2))
plot(GPA$Csize, pseudoGPA$Csize,xlab="Manual LMs Centroid Size",ylab="Pseudo LM Centroid Size",pch=20)
plot(GPA.pca$pc.scores[,1],pseudoGPA.pca$pc.scores[,1],xlab="Manual LMs PC1",ylab="Pseudo LM PC1",pch=20)
for (i in 1:5) print(cor.test(GPA.pca$pc.scores[,i],pseudoGPA.pca$pc.scores[,i]))

#plots and tables 

par(mfrow=c(1,1))

plot(joint.GPA$consensus[,1:2],pch=32,col="black",cex=2)
LM_labels=rownames(joint.GPA$consensus[, 1:2])
for (i in 1:51) points(joint.GPA$coords[, 1:2, i], pch=20, col="red")
for (i in 52:102) points(joint.GPA$coords[, 1:2, i], pch=20, col="cyan")
#points(joint.GPA$consensus[, 1:2], pch=3, col="black", cex=2, lwd=4)
pseudo_Mean=apply(joint.GPA$coords[, , 52:102], c(1, 2), mean)
GPA_Mean=apply(joint.GPA$coords[, , 1:51], c(1, 2), mean)
points(GPA_Mean[, 1:2], col="black", pch=3, cex=2, lwd=2)
text((GPA_Mean[, 1:2]+c(0, 0.001)), labels = LM_labels, pos=3, cex=1.5)
points(pseudo_Mean[, 1:2], col="black", pch=3, cex=2, lwd=2)

groups=as.factor(c(rep("GPA",51),rep("Pseudo",51)))
proc.anova=procD.lm(joint.GPA$coords~groups,iter = 1000,RRPP = T)

#Figure 5.
par(mfrow=c(1,2))
plot(round((GPA.pca$sdev^2)/sum(GPA.pca$sdev^2),3),pch=20,cex=2,cex.axis=2,main="GPA PCA",
     ylab="Variance Explained",xlab="Principal Components")
plot(round((pseudoGPA.pca$sdev^2)/sum(pseudoGPA.pca$sdev^2),3),pch=20,cex=2,cex.axis=2,main="Pseudo LM PCA",
     ylab="Variance Explained",xlab="Principal Components")


par(mfrow=c(1,1))
plot(GPA.pca$pc.scores[,1],pseudoGPA.pca$pc.scores[,1],pch=20, xlab="GPA PC1",ylab="Pseudo LM PC1")
for (i in 1:1) {
  print(cor.test(GPA.pca$pc.scores[,i],pseudoGPA.pca$pc.scores[,i]))}

#further analysis regarding the PD of landmarks

proc.Dist=function(aligned){
  return(sqrt(sum((aligned - joint.GPA$consensus)^2)))
}

#to assess the magnitude in the shape variation we will first look at the procrustes distance
#as a quick way measure of measuring total variation in a landmark configuration
PD=NULL
for (i in 1:dim(joint)[3]) PD[i]= proc.Dist(joint.GPA$coords[,,i])
tapply(PD, groups, sd)/tapply(PD, groups, mean)*100
t.test(PD[groups=="GPA"], PD[groups=="Pseudo"], paired=T)
wilcox.test(PD[groups=="GPA"], PD[groups=="Pseudo"], paired=T, alternative="greater") 

#we can also look at each individual landmark
#but we have to keep in mind that these are based on already aligned coordinates

LM_dist=array(dim=c(45,102))
for (i in 1:102) for (j in 1:45) LM_dist[j, i] = sqrt(sum((joint.GPA$coords[j,,i] - joint.GPA$consensus[j,])^2))

for (i in 1:45) print(wilcox.test(LM_dist[i, groups=="GPA"], LM_dist[i, groups=="Pseudo"], paired=T, alternative="greater") )

summary(proc.anova)