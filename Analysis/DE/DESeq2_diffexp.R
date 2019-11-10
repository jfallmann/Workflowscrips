#https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
library(DESeq2)
args <- commandArgs(trailingOnly = TRUE)

anname<-args[1]
inname<-args[2]
outcsv<-args[3]
outdir<-args[4]

anno <- as.matrix(read.table(anname,row.names=1))
colnames(anno) <- c("condition")
anno <- as.data.frame(anno)
head(anno)

countData <- as.matrix(read.table(inname,header=T,row.names=1))
head(countData)

#Check if names are consistent
all(rownames(anno) == colnames(countData))

#Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = anno,
                              design= ~ condition)

#filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)

res <- results(dds,contrast=c("condition",levels(anno$condition)), name=paste(levels(anno$condition),collapse='_vs_'))

#sort and output
resOrdered <- res[order(res$log2FoldChange),]
#mcols(res,use.names=TRUE)
#DataFrame with 6 rows and 2 columns
#                       type

#baseMean       intermediate
#log2FoldChange      results
#lfcSE               results
#stat                results
#pvalue              results
#padj                results
#description

#baseMean                        the base mean over all rows
#log2FoldChange log2 fold change (MAP): condition treated vs untreated
#lfcSE           standard error: condition treated vs untreated
#stat            Wald statistic: condition treated vs untreated
#pvalue          Wald test p-value: condition treated vs untreated
#padj            BH adjusted p-values
#write the table to a csv file

pdf(paste(inname,"DESeq2","plot.pdf",sep="_"))
plotMA(res, ylim=c(-3,3))
dev.off()
write.table(as.data.frame(resOrdered), file=paste(outcsv), sep="\t")

###
#Now we want to transform the raw discretely distributed counts so that we can do clustering. (Note: when you expect a large treatment effect you should actually set blind=FALSE (see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).
rld<- rlogTransformation(dds, blind=TRUE)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)

#We also write the normalized counts to file
write.table(as.data.frame(assay(rld)), file=paste(inname,"DESeq2_rld.txt",sep="_"), sep="\t", col.names=NA)
write.table(as.data.frame(assay(vsd)), file=paste(inname,"DESeq2_vsd.txt",sep="_"), sep="\t", col.names=NA)

#Here we choose blind so that the initial conditions setting does not influence the outcome, ie we want to see if the conditions cluster based purely on the individual datasets, in an unbiased way. According to the documentation, the rlogTransformation method that converts counts to log2 values is apparently better than the old varienceStabilisation method when the data size factors vary by large amounts.

par(mai=ifelse(1:4 <= 2, par('mai'), 0))
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord]<150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c('blue', 'black')
pdf(paste(inname,"DESeq2_VST","and_log2.pdf",sep="_"))
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type='l', lty=1, col=vstcol, xlab='n', ylab='f(n)')
legend('bottomright', legend = c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.off()
#The x axis is the square root of variance over the mean for all samples, so this will naturally included variance due to the treatment. The goal here is to flattern the curve so that there is consistent variance across the read counts, and that is what we got.
##############################
#
 #library('vsn')
 #par(mfrow=c(1,3))
 #notAllZero <- (rowSums(counts(dds))>0)
 #pdf(paste(inname,"SdPlot.pdf",sep="_"))
 #meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim = c(0,2.5))
 #meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,2.5))
 #meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,2.5))
 #dev.off()
#This interesting plot shows the standard deviation across all samples against the mean counts using three different methods of transformation. With this data you can see that the shifted logarithm method (left) seems to do pretty badly at for low count genes, with both regularized log (center) and DESeqs variance stabilisation (right) doing a much better job across the entire dynamic range of counts.
##############################
#
#
#
##############################
library('RColorBrewer')
library('gplots')
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
pdf(paste(inname,"DESeq2","heatmap1.pdf",sep="_"))
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
Rowv = FALSE, Colv = FALSE, scale='none',
dendrogram='none', trace='none', margin=c(10,6))
dev.off()
pdf(paste(inname,"DESeq2","heatmap2.pdf",sep="_"))
heatmap.2(assay(rld)[select,], col = hmcol,
Rowv = FALSE, Colv = FALSE, scale='none',
dendrogram='none', trace='none', margin=c(10, 6))
dev.off()
pdf(paste(inname,"DESeq2","heatmap3.pdf",sep="_"))
heatmap.2(assay(vsd)[select,], col = hmcol,
Rowv = FALSE, Colv = FALSE, scale='none',
dendrogram='none', trace='none', margin=c(10, 6))
dev.off()
#The above shows heatmaps for 30 most highly expressed genes (not necessarily the biggest fold change). The data is of raw counts (left), regularized log transformation (center) and from variance stabilizing transformation (right) and you can clearly see the effect of the transformation has by shrinking the variance so that we don’t get the squish effect shown in the left hand graph.
##############################
#Now we calculate sample to sample distances so we can make a dendrogram to look at the clustering of samples.
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
#updated in latest vignette (See comment by Michael Love)
#this line was incorrect
#heatmap.2(mat, trace='none', col = rev(hmcol), margin=c(16, 16))
#From the Apr 2015 vignette
hc <- hclust(distsRL)
pdf(paste(inname,"DESeq2","heatmap_samplebysample.pdf",sep="_"))
heatmap.2(mat, Rowv=as.dendrogram(hc),
symm=TRUE, trace='none',
col = rev(hmcol), margin=c(13, 13))
dev.off()

##############################
pdf(paste(inname,"DESeq2","PCA.pdf",sep="_"))
print(plotPCA(rld, intgroup=c('condition')))
dev.off()

##############################
#Previously we talked about the cooks distance treatment of outliers, in that a gene is thrown away if one of its samples is deemed to be an outlier. You may not want this to happen so DESeq2 we can take a different approach by replacing the outlier value with one estimated value as predicted by the distribution using the trimmed mean approach. DESeq2 recomends you only do this if you have several replicates per treatment, and indeed it automatically uses this feature if you have 7 or more replicates in your datatable.
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < .1,
cleaned = results(ddsClean)$padj < .1)
addmargins(tab)
write.table(as.data.frame(tab),file=paste(inname,'deseq2','initial'), sep="\t")
resClean <- results(ddsClean)
write.table(as.data.frame(resClean),file=paste(inname,'deseq2','cook_dist_cleaned'), sep="\t")

##############################
#Dispersion plot shows how the estimates are shrunk from the gene wise values (black dots) toward the fitted estimates, with the final values used in testing being the blue dots.
pdf(paste(inname,"DESeq2","PCA.pdf",sep="_"))
plotDispEsts(dds)
dev.off()
##############################
##Now independent filtering to remove any tests that have little chance of pass to reduce the number of tests we have to perform, thus reducing the effects of multiple testing error. (false discovery). You can see how many genes are rejected based on different values of alpha (FDR)
##filtering threshold
#attr(res,'filterThreshold')
#plot(attr(res,'filterNumRej'),type='b', ylab='number of rejections')
#dev.copy(png,'deseq2_filtering_treshold.png')
#dev.off()

##############################
#Plot of the maximum Cook’s distance per gene over the rank of the Wald statistics for the condition.

W <- res$stat
maxCooks <- apply(assays(dds)[['cooks']],1,max)
idx <- !is.na(W)
pdf(paste(inname,"DESeq2","Cooksdist.pdf",sep="_"))
plot(rank(W[idx]), maxCooks[idx], xlab='rank of Wald statistic',ylab='maximum Cook\'s distance per gene', ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))
dev.off()
#Here more about independent filtering. What it shows in genes with very low counts are unlikely to have a significant p-value due to excessive dispersion at the left side of the dynamic range of counts. The y-axis here is -log10, so bigger numbers are smaller p-values (better).

##############################
pdf(paste(inname,"DESeq2","indep_filt.pdf",sep="_"))
plot(res$baseMean+1, -log10(res$pvalue),
     log='x', xlab='mean of normalized counts',
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
dev.off()
#All those dots on the left hand side the graph represent failed tests due to very low count values, thus we can really just get rid of them to reduce our chance of making a type I error.
#And again, you can see that only a few small (or no) p-values are discarded by the filtering.
##############################
#The graph on the left ranks the p-values from smallest to biggest (x-axis) and plots them. The black line is the actual p-value numbers (remember only about 23 genes had a p-value lower than 0.05). The red line has a slope that represents the number of tests divided by the false discovery rate (0.1). The idea here is the FDR is controlled at the 0.1% value for all tests that occur to the left of the right-most intersection of the black and red line.
use <- res$baseMean > attr(res,'filterThreshold')
table(use)
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c('do not pass'='khaki', 'pass'='powderblue')
pdf(paste(inname,"DESeq2","filterbar.pdf",sep="_"))
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = '', ylab='frequency')
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend('topright', fill=rev(colori), legend=rev(names(colori)))
dev.off()
#resFilt <- res[use & !is.na(res$pvalue),]
#orderInPlot <- order(resFilt$pvalue)
#showInPlot <- (resFilt$pvalue[orderInPlot] <= 0.08)
#alpha <- 0.1
#plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
#     pch='.', xlab = expression(rank(p[i])), ylab=expression(p[i]))
#abline(a=0, b=alpha/length(resFilt$pvalue), col='red3', lwd=2)
