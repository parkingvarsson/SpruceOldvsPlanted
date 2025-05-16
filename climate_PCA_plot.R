library(corrplot)
library(psych)

cols<-rep(c("forestgreen","tan2","tan2"),15)
pca<-prcomp(samples[,c(8,9,10,12,13,14,17)],scale=T)

plot(pca$x[,c(1,2)],pch=19,col=cols,cex=1.8,xlab="Climate PC1",ylab="Climate PC2",las=1)
abline(v=0,lty=3)
abline(h=0,lty=3)
sf<-3
arrows(x0 = 0, x1 = sf*pca$rotation[,1], 
       y0 = 0, y1 = sf*pca$rotation[,2], 
       col = "blue", 
       length = 0.1, 
       lwd = 1.5,
       angle = 30)

text(x = sf*pca$rotation[,1], y = sf*pca$rotation[,2], 
     labels = c("AIT","CMT","CONT","DegD0","DegD5","MTC","PCQ"), 
     cex = 1.2,
     font = 2,
     col = "black", 
     pos = c(2, 4, 2, 1, 3, 1))

corrplot.mixed(cor(samples[,c(7:22)]), upper = "ellipse", lower = "number",
         title = "Environmental variables",
         tl.pos = "n", mar = c(2, 1, 3, 1)) 

pairs(samples[,c(7:22)],cor=TRUE)

samples.new<-samples[,c(7:22)]

names(samples.new)<-c("annual PET","AIT","CMI","CONT","embQ","DegD0","DegD5","MTC","MTW","MT10","PETc","PETd","PETs","PETwaq","PETweq","TIdx")


pairs.panels(samples.new,
             smooth = FALSE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = FALSE,    # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             pch = 21,           # pch symbol
             lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             hist.col = 4,       # Histograms color
             stars = TRUE,       # If TRUE, adds significance level with stars
             ci = TRUE)          # If TRUE, adds confidence intervals
