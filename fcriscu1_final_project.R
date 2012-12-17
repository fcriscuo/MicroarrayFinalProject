# Final Project
# AS.410.671 Microarrays & Analysis
# Fred Criscuolo   fcriscu1@jhu.edu
#
# An analysis pipeline for expression data submitted to NCBI GEO data repository
# GEO accession: GDS3096 Inflammatory breast cancer: stroma
# An analysis of the stromata surrounding tumors from patients with
# inflammatory breast cancer (IBC) and invasive non-IBC
# http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3096
# Ref: PMID 17999412
# 47 samples
#
# setup

source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
biocLite("impute")
biocLite("pcaMethods")
biocLite("genefilter")
biocLite("multtest")
biocLite("MASS")
biocLite("ssize")

library(ssize)
library(gdata)
library(GEOquery)
library(impute)
library(pcaMethods)
library(genefilter)
library(multtest)
library(MASS)

install.packages("calibrate")
library(calibrate)

install.packages("outliers")
library(outliers)

#SVM package for feature selection
#re: http://www.uccor.edu.ar/paginas/seminarios/Software/SVM_RFE_R_implementation.pdf
install.packages("e1071")
library(e1071)
source(file="SVM-RFE\\code\\SVM-RFE.r")

library(ggplot2)
# function definitions
# function to change column names based on sample type
new.column.name<-function(id, df) {
  loc <- grep(id, Columns(df)$sample)
  if( length(loc) < 1 ) {id} else {
    des <- Columns(df)$disease.state[loc]
    if (regexpr("non-inflammatory",des) < 1) { paste (id,"_IBC",sep="") }
    else {paste (id,"_nonIBC",sep="")}
    
  }
}


setwd("C:/JHU Program/Courses/AS.410.671 Microarrays & Analysis/Final Project")
# use BioConductor GEOQuery pakage to load the data from NCBI
gds <- getGEO("GDS3096")
exp.data <- Table(gds)

#read in the annotation file and keep these data in sync with the expression data
annon.file.name <-  "GDS3096_breast_stroma/GPL96.annot"
annon.data <- read.delim(annon.file.name,header=T, row.names=1,skip=27,sep="\t")
#edit out the annotations for the Affy controls
annon.data <- annon.data[1:nrow(exp.data),]
gene.info<- data.frame(title=annon.data$Gene.title,symbol=annon.data$Gene.symbol)
rownames(gene.info) <- rownames(annon.data)

# rename the data columns to indicate IBC or non-IBC
new.names <- c()
ibc.sample.count <-0
for (name in names(exp.data)) { 
  new.name <- new.column.name(name,gds)
  if(regexpr("non",new.name) < 1){ibc.sample.count <- ibc.sample.count +1}
  new.names <- c(new.names,new.name)}
names(exp.data) <- new.names
cat("\nItem 01: New data column names:\n")
cat (names(exp.data))
nonibc.sample.count <- ncol(exp.data) - ibc.sample.count -2
cat(sprintf("nonIBC sample count = %d; IBC sample count = %d", nonibc.sample.count, ibc.sample.count))

gsm.matrix <- data.matrix(exp.data[,(3:ncol(exp.data))])
#Descriptive satistics for initial expression data
sample.count<- nrow(gsm.matrix)
exp.df <- exp.data[,(3:ncol(exp.data))]
exp.data.mean <- rowMeans(gsm.matrix,na.rm=TRUE)
exp.data.sd<- apply(exp.df,1,sd,na.rm=TRUE)

mean.title<-paste("Histogram of mean log2 expression values \nfor",sample.count,"samples")
sd.title<-paste("Histogram of standard deviations \nfor log2 expression values for\n",sample.count,"samples")
#par(mfrow=c(1,2))
png("des_mean_plot.png")
hist(exp.data.mean, col="orange",
     xlab="Mean log2 expression value for IBC/nonIBC samples",
     ylab="Frequency",
     main=mean.title)
dev.off()
png("des_sd_plot.png")
hist(exp.data.sd, col="blue",
     xlab="StdDev log2 expression value for IBC/nonIBC samples",
     ylab="Frequency",
     main=sd.title)
dev.off()

#fold change
fc <- log2(3.0)
power<-0.8
sig.level <- 0.001
all.size <- ssize( sd=exp.data.sd, delta=fc, sig.level=sig.level, power=power)
ssize.plot( all.size, lwd=2, col="magenta", xlim=c(1,20), marks=c(1,5,10,15,20))
xmax<- par("usr")[2]-1
ymin<-par("usr")[3]+0.05

# Task 1- Identify and remove outliers

cat(sprintf("\nDimensions for edited data matrix: %d rows, %d columns",nrow(gsm.matrix), 
            ncol(gsm.matrix)))
# generate correlation matrix
gsm.corr.matrix <- cor(gsm.matrix, method="pearson", use="pairwise.complete.obs")
gsm.heatmap <- heatmap(gsm.corr.matrix, col=heat.colors(256), scale="column",
                       xlab="GSM data", ylab="GSM Data", 
                       main="GSM expresson")

#cluster plot
t.gsm.matrix <- t(gsm.matrix)
dat.dist <-dist(t.gsm.matrix,method="euclidean")
dat.clust <- hclust(dat.dist, method="single")
plot(dat.clust,labels=names(t.gsm.matrix), cex=0.75, xlab="Breast Stroma Samples",
     main="Cluster Dendrogram for breast stroma samples: \nnonIBC and IBC")

#Plot the coefficient of variance vs the mean
gsm.col.mean <- apply(log2(gsm.matrix),2, mean) # column (2) means
gsm.col.var <- apply(log2(gsm.matrix),2, var)  # column variances
gsm.col.cv <- gsm.col.var / gsm.col.mean

plot(gsm.col.mean, gsm.col.cv, xlab="Mean of log2(Intensity)",
     ylab="Coefficient of Variation log2(Intensity)",
     main="CV vs Mean Plot for IBC Intensity Measurements",
     col=c(rep("red",nonibc.sample.count),rep("blue",ibc.sample.count)),
     pch=c(rep(16,nonibc.sample.count),rep(17,ibc.sample.count)))
legend("topright", c("nonIBC","IBC"), pch=c(16,17), col=c("red","blue"),bg="lightyellow")
# use textxy() function from calibrate package to label points (cleaner than text() function)
textxy(gsm.col.mean, gsm.col.cv,names(exp.data)[3:ncol(exp.data)])


#Avg correlation plot
gsm.corr.row.mean <- apply(log2(gsm.corr.matrix),1, mean) #  correlation matrix row (1) means

par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(gsm.corr.row.mean)), range(gsm.corr.row.mean), type="n", xlab="",ylab="Average r",
     main="Avg correlation of tumor/normal samples",
     axes=F)
points(gsm.corr.row.mean,
       bg="red",
       ,col=c(rep("red",nonibc.sample.count),rep("blue",ibc.sample.count)), 
       pch=c(rep(16,nonibc.sample.count),rep(17,ibc.sample.count)),
       cex=1.25)
axis(1, at=c(1:length(gsm.corr.row.mean)), labels=dimnames(gsm.matrix)[[2]],
     las=2, cex.lab=0.4, cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")
legend("bottomleft", c("nonIBC","IBC"), pch=c(16,17), col=c("red","blue"),bg="lightyellow")

#use outlier() function from outliers package to identify largest outlier
x <- gsm.corr.row.mean <=  outlier(gsm.corr.row.mean)
y <- gsm.corr.row.mean[x]
cat(sprintf("Outlier function identified %s as the largest outlier ",names(y)))

#Based on the results from the above outlier detection methods 
# remove IBC smaple GSM136326_IBC from the data used for futher analysis
gsm.matrix.ed <- gsm.matrix[, -(grep("GSM136326_IBC",colnames(gsm.matrix)))]
cat(sprintf("\nNumber of columns after removing GSM136326_IBC:  %d",ncol(gsm.matrix.ed)))

########################### Task 02 ##################################################
#Identify and remove genes that have low expression values
#
# standardize the rownames 
rownames(gsm.matrix.ed) <- rownames(annon.data)
exp.df <- data.frame(exp=rowMeans(gsm.matrix.ed))
#plot histogram of mean gene eexpression values
ggplot(exp.df, aes(x=exp)) +geom_histogram(binwidth=0.25,colour="black",fill="white") +
  geom_vline(aes(xintercept=median(exp)), color="red", size=1)+
  geom_vline(aes(xintercept=3.5), color="blue", linetype="dashed", size=1)+
  xlab("Mean Expression Value per Gene") + ylab("Count") +
  theme(panel.background = element_rect(fill='lightyellow')) +
  ggtitle("Distribution of Mean Gene Expression Values for GDS3096")

# based on distribution of mean expression values, remove genes whose mean expression value
# is <= 3.5
cutoff <- 3.5
exp.flag <- lapply(rowMeans(gsm.matrix.ed), function(x){x>cutoff})
gsm.exp.matrix <- gsm.matrix.ed[as.logical(exp.flag),]
# remove the same entries from the annotation data to keep that data frame in sync
gene.info <- gene.info[as.logical(exp.flag),]
#rownames are no longer consecutive - keep track of them
exp.rownames <- rownames(gsm.exp.matrix)
rownames(gsm.exp.matrix) <- rownames(gene.info)


##########################    Task 03 #################################################
# feature selection using two factor test
#

# Use Welch's t-test as a 2 factor selection analysis

ibc.type <- lapply(colnames(gsm.exp.matrix), function(x) {
if(regexpr("nonIBC",x) <1) {"IBC"} else { "nonIBC"}
})



pval <- c()
#calculate individual p-values
for (i in seq(nrow(gsm.exp.matrix))){
  tt <- t.test(gsm.exp.matrix[i,ibc.type=="nonIBC"],gsm.exp.matrix[i,ibc.type=="IBC"], alternative="two.sided")
  pval <- c(pval, tt$p.value)
}

cat(sprintf("The number of probesets whose p-value < %f is %d",
            threshold, sum(pval<threshold)))

#add p-value data to gene.info data fame
gene.info$pvalue <- pval
#test application of Bonferroni correction
bc <- threshold / nrow(gsm.exp.matrix)

cat(sprintf("Bonferroni correction: the number of probesets with a p-value < %10.8f = %d",
          bc, sum(pval <bc)))


# Output  the probesets that meet the threhold criterion to a csv file
gi <- gene.info[order(gene.info$pvalue),]
sink("pvalue_probes.csv")
for (r in seq(sum(pval<threshold))) {
  
  cat(sprintf("\n %s ,  %s ,  %8.6f ",
              rownames(gi)[r], gi$symbol[r],
              gi$pvalue[r] ))
  
}
sink()

#create a data frame to represent pass/fail for p-value criterion
pval.criterion <-lapply(as.list(pval),function(p){ if (p < threshold ) TRUE else FALSE})
pval.df <- as.data.frame(do.call(rbind,pval.criterion),row.name=exp.rownames)
names(pval.df) <- c("PvalCriterion")
rownames(pval.df) <- rownames(gene.info)

#plot histogram of gene p-values
p.df <-data.frame(pvalue=pval)
bin.width<-0.01
exp.genes.per.bin <- nrow(gsm.exp.matrix) / (1.0 / bin.width)
annotation.text.1 <- paste("Blue line = ",threshold, " cutoff", sep="")
annotation.text.2 <- paste("Red line = 0.05 cutoff",sep="")
annotation.text.3 <- paste("Green line = the expected number of probesets\nper bin by chance (",
                           exp.genes.per.bin,")",sep="")
title.text <-"Distribution of Student t-test p-Values by probeset for GDS3096 data set\nIBC vs. non-IBC"

ggplot(p.df, aes(x=pvalue)) +geom_histogram(binwidth=bin.width,colour="black",fill="white") +
  xlab("p-value by Gene") + ylab("Count") +
  geom_vline(aes(xintercept=threshold), color="blue", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=0.05), color="red", linetype="dashed", size=1)+
  geom_hline(aes(yintercept=exp.genes.per.bin), colour="green",size=1) +
  theme(panel.background = element_rect(fill='lightblue')) +
  annotate("text",x=0.40,y=400, label=annotation.text.1, colour="blue")+
  annotate("text",x=0.40,y=350, label=annotation.text.2, colour="red")+
  annotate("text",x=0.50,y=300, label=annotation.text.3, colour="darkgreen")+
  ggtitle(title.text)

# determine probesets that demonstrate a abs 1.5x fold change between IBC and nonIBC samples
fold.change.criterion <- log2(1.5)
ibc.exp.matrix <- gsm.exp.matrix[,ibc.type=="IBC"]
ibc.exp.mean <- apply(ibc.exp.matrix,1,mean, na.rm=TRUE)
nonibc.exp.matrix <- gsm.exp.matrix[,ibc.type=="nonIBC"]
nonibc.exp.mean <- apply(nonibc.exp.matrix,1,mean, na.rm=TRUE)
fold.change <- ibc.exp.mean - nonibc.exp.mean
fold.criterion <- lapply (as.list(fold.change), function(x){ if (abs(x)> fold.change.criterion) TRUE else FALSE })
fold.df <- as.data.frame(do.call(rbind,fold.criterion))
rownames(fold.df) <- rownames(gene.info)
names(fold.df) <- c("FoldCriterion")
cat(sprintf("The nummber of probesets that have a 1.5x fold change between IBC and non IBC is %d",
    sum(fold.df$FoldCriterion == TRUE)))
#merge the two boolean dataframes to determine the probesets that meet both criteria
# selected.df <- merge(fold.df, pval.df, by="row.names",sort=TRUE)
# add fold change to gene.info data frame
gene.info$fold.change <- fold.change

#plot fold change histogram
fc.title.text <-"Distribution of log2 fold change by probeset for GDS3096 data set\nIBC vs. non-IBC"
fc.annotation <- "Red lines indicate \nfold change significance cutoff"

ggplot(gene.info, aes(x=fold.change)) +geom_histogram(binwidth=bin.width,colour="black",fill="white") +
  xlab("log2 Fold Change by Probeset") + ylab("Probeset Count") +
  theme(panel.background = element_rect(fill='lightyellow')) +
  annotate("text",x=-1.2,y=700, label=fc.annotation, colour="red")+
  geom_vline(aes(xintercept=fold.change.criterion), color="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=-fold.change.criterion), color="red", linetype="dashed", size=1)+
  ggtitle(fc.title.text)

gene.info$FoldCriterion <- fold.df$FoldCriterion
gene.info$PvalCriterion <- pval.df$PvalCriterion



selected.probesets <- subset(gene.info, (FoldCriterion & PvalCriterion ))

# create report for selected probesets
selected.probesets<-selected.probesets[order(selected.probesets$pvalue),]
for (r in seq(nrow(selected.probesets))) {
  cat(sprintf("\nProbe: %s \t Symbol: %s \tp-value  %8.6f \t Fold change: %6.6f",
              rownames(selected.probesets)[r], selected.probesets$symbol[r],
              selected.probesets$pvalue[r],selected.probesets$fold.change[r]  ))
}




#select the genes that meet the cutoff criterion
# keep track of their original index value and p-value
top.genes <-order(pval)[1:gene.count]
top.exp.df <- data.frame(gene.index=top.genes,
                         exp.matrix=gsm.exp.matrix[top.genes,],
                         p.value=pval[top.genes])

top.exp.matrix <-gsm.exp.matrix[top.genes,]
rownames(top.exp.matrix) <-rownames(gsm.exp.matrix[top.genes,])


#Feature Selection using SVM-RFE algorithm implementation
# provided by http://www.uccor.edu.ar/paginas/seminarios/Software/SVM_RFE_R_implementation.pdf
#scale samples with mean zero and standard deviation one
gsm.exp.scaled <- t(top.exp.matrix)  # transform expression matrix
for(i in 1:nrow(gsm.exp.scaled))gsm.exp.scaled[i,]<- 
   (gsm.exp.scaled[i,]-mean(gsm.exp.scaled[i,]))/sd(gsm.exp.scaled[i,])
#scale genes (i.e.) with mean zero and standard deviation one
for(i in 1:ncol(gsm.exp.scaled))gsm.exp.scaled[,i] <-
  (gsm.exp.scaled[,i]-mean(gsm.exp.scaled[,i]))/sd(gsm.exp.scaled[,i])
gsm.exp.scaled<-  2*atan(gsm.exp.scaled/2)

# Feature Ranking with SVM-RFE
fact <- factor(unlist(ibc.type))
topfeatureRankedList = svmrfeFeatureRanking(gsm.exp.scaled,
                                         fact)
# computation of the feature rankings is compute intensive
# Output the results to a file for reuse if necessary
df <- data.frame(rank=topfeatureRankedList)
write.table(df,file="top_svmrfe_rankings.txt")


#plot concordance between the two ordered lists 

plot.df <- data.frame(pval.order=seq(1:gene.count), 
                      feat.rank=topfeatureRankedList[1:gene.count])
plot2.title.text <-paste("Feature Selection Concordance: t-test vs SVM-RFE rankings",
                         "\nTop ",gene.count," rankings",sep="")
ggplot(plot.df, aes(x=pval.order,y=feat.rank)) +
  geom_point(colour="red") +
  xlab("T-test p-value Order") + ylab("SVM-RFE rankings") +
  theme(panel.background = element_rect(fill='lightgreen')) +
  ggtitle(plot2.title.text)

#Select the probesets that passed the Welch t-test criterion using the order
# calculated by SVM-RFE refinement
svmrfe.top.exp.matrix <- top.exp.matrix[topfeatureRankedList,]
rownames(svmrfe.top.exp.matrix) <- rownames(top.exp.matrix[topfeatureRankedList,])
top.gene.info <- gene.info[rownames(svmrfe.top.exp.matrix),]
#plot histogram of p-values for selected probesets
pv.title.text<- paste("Histogram of Welch's t-test p-values for top ",nrow(top.gene.info)," probesets",
                                      sep="")
pv.annot01<-"p-value selection criterion: < 0.001"
ggplot(top.gene.info, aes(x=pvalue)) +geom_histogram(binwidth=0.00005,colour="blue",fill="lightblue") +
  xlab("Welch's t-test p-value") + ylab("Probeset Count") +
  theme(panel.background = element_rect(fill='lightyellow')) +
  annotate("text",x=0.0005,y=5.5, label=pv.annot01, colour="red")+
  ggtitle(pv.title.text)

####################################Task 4###########################################
#
# Sample classification
#
# PCA analysis of samples for top probesets
#
top.probe.pca <- prcomp(t(svmrfe.top.exp.matrix), cor=F)
pca.loadings <- top.probe.pca$x[,1:3]

mycols <- c("deeppink", "blue2")
plot.pch <-c(1,2)
png("pca_plots.png")
par(mfrow=c(3,1))

#n.b. the range for y axis was expanded beyond the range of the values to provide room for the legend

plot(range(pca.loadings[,1]), range(c(-5,5)),xlab="Principal Component 1",col.lab="darkblue",
     ylab="Principal Comonent 2", 
     main="PCA plot for GDS3096 data IBC and nonIBC samples\np1 vs. p2",col.main="darkblue")
points(pca.loadings[,1][as.numeric(ibc.type=="nonIBC")==1], (pca.loadings[,2][as.numeric(ibc.type=="nonIBC")==1]),
       col=mycols[1], pch=plot.pch[1], cex=1.5)
points(pca.loadings[,1][as.numeric(ibc.type=="IBC")==1], (pca.loadings[,2][as.numeric(ibc.type=="IBC")==1]),
       col=mycols[2], pch=plot.pch[2], cex=1.5)
legend("bottomright", ,c("nonIBC", "IBC"), col=mycols, text.col=mycols,pch=plot.pch)



plot(range(pca.loadings[,1]), range(pca.loadings[,3]),xlab="Principal Component 1",col.lab="darkblue",
     ylab="Principal Comonent 3", 
     main="PCA plot for GDS3096 data IBC and nonIBC samples\np1 vs. p3",col.main="darkblue")
points(pca.loadings[,1][as.numeric(ibc.type=="nonIBC")==1], (pca.loadings[,3][as.numeric(ibc.type=="nonIBC")==1]),
       col=mycols[1], pch=plot.pch[1], cex=1.5)
points(pca.loadings[,1][as.numeric(ibc.type=="IBC")==1], (pca.loadings[,3][as.numeric(ibc.type=="IBC")==1]),
       col=mycols[2], pch=plot.pch[2], cex=1.5)
legend("bottomright", ,c("nonIBC", "IBC"), col=mycols, text.col=mycols,pch=plot.pch)




plot(range(pca.loadings[,2]), range(pca.loadings[,3]),xlab="Principal Component 2",col.lab="darkblue",
     ylab="Principal Comonent 3", 
     main="PCA plot for GDS3096 data IBC and nonIBC samples\np2 vs. p3",col.main="darkblue")
points(pca.loadings[,2][as.numeric(ibc.type=="nonIBC")==1], (pca.loadings[,3][as.numeric(ibc.type=="nonIBC")==1]),
       col=mycols[1], pch=plot.pch[1], cex=1.5)
points(pca.loadings[,2][as.numeric(ibc.type=="IBC")==1], (pca.loadings[,3][as.numeric(ibc.type=="IBC")==1]),
       col=mycols[2], pch=plot.pch[2], cex=1.5)
legend("bottomright", ,c("nonIBC", "IBC"), col=mycols, text.col=mycols,pch=plot.pch)

dev.off()

#Scree plot to determine the degree of variability 
ibc.pca.var <- round(top.probe.pca$sdev^2 / sum(top.probe.pca$sdev^2)*100,2)
plot(c(1:length(ibc.pca.var)), ibc.pca.var, type="b", xlab="# of components",col.lab="darkblue",
     ylab="% variance", 
     main="Scree plot of IBC sample variance", col="blue",col.main="darkblue")

#MDS analysis
top.probe.dist <- dist(t(svmrfe.top.exp.matrix))
top.probe.loc <- cmdscale(top.probe.dist)
top.probe.mds <- isoMDS(top.probe.dist)
png("mds_plots.png")
par(mfrow=c(3,1))
# plot classical MDS
plot(top.probe.loc, type="n")
points(top.probe.loc[,1][as.numeric(ibc.type=="nonIBC")==1], 
       (top.probe.loc[,2][as.numeric(ibc.type=="nonIBC")==1]),
       col=mycols[1], pch=plot.pch[1], cex=1.5)
points(top.probe.loc[,1][as.numeric(ibc.type=="IBC")==1], 
       (top.probe.loc[,2][as.numeric(ibc.type=="IBC")==1]),
       col=mycols[2],pch=plot.pch[2], cex=1.5)
title(main="Classical (metric) MDS plot of IBC data",col.main="darkblue",col.lab="darkblue")
legend("bottomright",c("nonIBC", "IBC"), col=mycols, text.col=mycols,pch=plot.pch)

#Plot Kruskal's non-metric MDS
plot(top.probe.mds$points, type="n")
points(top.probe.mds$points[,1][as.numeric(ibc.type=="nonIBC")==1], 
       (top.probe.mds$points[,2][as.numeric(ibc.type=="nonIBC")==1]),
       col=mycols[1], pch=plot.pch[1], cex=1.5)
points(top.probe.mds$points[,1][as.numeric(ibc.type=="IBC")==1], 
       (top.probe.mds$points[,2][as.numeric(ibc.type=="IBC")==1]),
       col=mycols[2], pch=plot.pch[2], cex=1.5)
title(main="Kruskal (non-metric) MDS plot of IBC data",col.main="darkblue", col.lab="darkblue")
legend("bottomright",c("nonIBC", "IBC"), col=mycols, text.col=mycols,pch=plot.pch)
dev.off()

ibc.t <-t(svmrfe.top.exp.matrix)
ibc.scale <- scale(ibc.t, center=TRUE, scale=TRUE)
#weighted graph Laplacian
# n.b. k.specClust2 function copied from Lecture 8 R code
k.speClust2 <- function (X, qnt=NULL) {
  dist2full <- function(dis) {
    n <- attr(dis, "Size")
    full <- matrix(0, n, n)
    full[lower.tri(full)] <- dis
    full + t(full)
  }
  dat.dis <- dist(t(X),"euc")^2
  if(!is.null(qnt)) {eps <- as.numeric(quantile(dat.dis,qnt))}
  if(is.null(qnt)) {eps <- min(dat.dis[dat.dis!=0])}
  kernal <- exp(-1 * dat.dis/(eps))
  K1 <- dist2full(kernal)
  diag(K1) <- 0
  D = matrix(0,ncol=ncol(K1),nrow=ncol(K1))
  tmpe <- apply(K1,1,sum)
  tmpe[tmpe>0] <- 1/sqrt(tmpe[tmpe>0])
  tmpe[tmpe<0] <- 0
  diag(D) <- tmpe
  L <- D%*% K1 %*% D
  X <- svd(L)$u
  Y <- X / sqrt(apply(X^2,1,sum))
}
qnt<-NULL
phi <- k.speClust2(t(ibc.scale),qnt)  
if( !is.null(qnt)) {
  plot.title <-paste("Weighted Graph Laplacian plot of IBC Data\nepsilon="
                     ,qnt,sep="")
} else {
  plot.title <-"Weighted Graph Laplacian plot of IBC Data"
}

plot(range(phi[,1]),range(phi[,2]),xlab="phi1",ylab="phi2",
     main=plot.title,col.main="darkblue")
points(phi[,1][as.numeric(ibc.type=="nonIBC")==1],phi[,2][as.numeric(ibc.type=="nonIBC")==1],col="red",
       pch=plot.pch[1],cex=1.5)
points(phi[,1][svmrfe.top.exp.matrix],phi[,2][as.numeric(ibc.type=="IBC")==1],col="blue",
       pch=plot.pch[2],cex=1.5)
legend("topright",c("nonIBC", "IBC"),col=c("red", "blue"),pch=plot.pch,cex=.7,horiz=F)

# Task 05 - Cluster Analysis
#
#
# perform hierarchical clustering in two dimenstions (i.e. genes and samples)
# plot a heatmap; y-axis=genes, x-axis=samples

#calculate a distance matrix using the manhattan method
sample.dist.man <- dist(t(svmrfe.top.exp.matrix), method="manhattan")

#generate a hierachical cluster using median
hc.med.sample <-  hclust(sample.dist.man, method="median")
hc.med.sample$labels <- lapply(ibc.type, function(x){
  if (x =="IBC") {"I"}else {"N"}
})
                         
plot(hc.med.sample, main="Cluster Dendogram by IBC status \n(Manhattan distance method)",
     col.main="darkblue",
     xlab="Clusters by IBC status",ylab="Distance")
legend("bottomleft",c("N=nonIBC", "I=IBC") ,col="blue")

# color vector obtained from lecture notes
hm.rg <- c("#FF0000","#CC0000","#990000","#660000","#330000","#000000",
           "#000000","#0A3300","#146600","#1F9900","#29CC00","#33FF00")
tmp.matrix<-svmrfe.top.exp.matrix
colnames(tmp.matrix)<-ibc.type
hm.title <- paste("heatmap for top",gene.count,"probesets and IBC status")
heatmap(tmp.matrix,col=hm.rg,
        main=hm.title,col.main="darkblue"
        xlab="IBC status ", 
        ylab="Topr ranked probesets")

#perform k-means clustering using top 2 eigenvectors from PCA analysis
top.probe.eigen <- top.probe.pca$rotation[,1:2]
top.probe.kmeans <- kmeans(top.probe.eigen,centers=2, iter.max=30)
print(top.probe.kmeans$cluster)

# Plot a 2-D scatter plot of the sample classification labels
# embedded within the first two (2) eigenvectors
pca.class <- lapply(rownames(top.probe.pca$x), function(x){
  if( regexpr("non",x)>0) "nonIBC"else "IBC"
})
sample.classification <- lapply(top.probe.kmeans$cluster,function(x){
  if (x==1){"nonIBC"}else{"IBC"}
})
plot.df <- data.frame(pc1=top.probe.pca$rotation[,1], 
                      pc2=top.probe.pca$rotation[,2], 
                      classification=factor(unlist(sample.classification)))

ggplot(plot.df, aes(x=as.numeric(pc1), y=as.numeric(pc2), colour=classification)) +
  geom_point(shape=1)+
  xlab("Top Probesests Eigenvector 1") + ylab("Top Probesets EigenVector 2") +
  theme(panel.background = element_rect(fill='lightgreen')) +
  ggtitle("k-means classification using \nPCA Eigenvector 1 and Eigenvector 2")



#LDA analysis
lda.pred <- function(lda.matrix, start.index, stop.index,ibc.train.size=2 ){
  
  cat(sprintf("\nLDA processing %d samples %d arrays start=%d stop=%d",
              nrow(lda.matrix), ncol(lda.matrix), start.index, stop.index))  
  n.ibc <- ncol(lda.matrix[,as.numeric(ibc.type=="IBC")==1]) 
  n.nonibc <- ncol(lda.matrix[,as.numeric(ibc.type=="nonIBC")==1])
  
  t.exp.matrix <- t(lda.matrix)
  
  nonibc.train.size=3*ibc.train.size
  ibc.test.size <- n.ibc - ibc.train.size
  nonibc.test.size <- n.nonibc - nonibc.train.size
  
  training.map <- c(rep(TRUE,nonibc.train.size),rep(FALSE,nonibc.test.size),
                    rep(TRUE,ibc.train.size),rep(FALSE,ibc.test.size))
  
  testing.map <- !(training.map)
  
  train.class.names <- c( rep("nonIBC", nonibc.train.size), rep("IBC", ibc.train.size))
  test.class.names <- c( rep("nonIBC", nonibc.test.size), rep("IBC", ibc.test.size))
  
  class <- factor(train.class.names)
  training.data <- t.exp.matrix[training.map,]
  testing.data <- t.exp.matrix[testing.map,]
  train.lda <- lda(training.data[, start.index:stop.index],grouping=class)
  # predict using test data
  
  dat.pred <- predict(train.lda, testing.data[,start.index:stop.index],dimen=2)
  
  
  return (table(dat.pred$class,test.class.names))
}


small.pred <-lda.pred(svmrfe.top.exp.matrix,1,2,4)
top11.pred <-lda.pred(svmrfe.top.exp.matrix,1,11,4)

top.pred <-lda.pred(svmrfe.top.exp.matrix,1,nrow(svmrfe.top.exp.matrix),4)
all.pred<-lda.pred(gsm.exp.matrix,1,nrow(svmrfe.top.exp.matrix),4)





lda.matrix <- svmrfe.top.exp.matrix
ibc.matrix<- lda.matrix[,as.numeric(ibc.type=="IBC")==1]
n.ibc <- ncol(ibc.matrix)
nonibc.matrix<- lda.matrix[,as.numeric(ibc.type=="nonIBC")==1]
n.nonibc <- ncol(nonibc.matrix)
t.exp.matrix <- t(lda.matrix)

ibc.train.size=2
nonibc.train.size=3*ibc.train.size
ibc.test.size <- n.ibc - ibc.train.size
nonibc.test.size <- n.nonibc - nonibc.train.size

training.map <- c(rep(TRUE,nonibc.train.size),rep(FALSE,nonibc.test.size),
                  rep(TRUE,ibc.train.size),rep(FALSE,ibc.test.size))
testing.map <- !(training.map)

train.class.names <- c( rep("nonIBC", nonibc.train.size), rep("IBC", ibc.train.size))
test.class.names <- c( rep("nonIBC", nonibc.test.size), rep("IBC", ibc.test.size))
class <- factor(train.class.names)
training.data <- t.exp.matrix[training.map,]
testing.data <- t.exp.matrix[testing.map,]
train.lda <- lda(t.exp.matrix[training.map, 1:2],grouping=class)
# predict using test data

dat.pred <- predict(train.lda, t.exp.matrix[testing.map,1:2])
table(dat.pred$class,test.class.names)












