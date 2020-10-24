library(corrplot)
library(Hmisc)
library(correlation)
library(ggcorrplot)
library(ggplot2)
my_data <- read.csv("mets_corr.csv")
cormat <- correlation(
  my_data,
  method = "pearson",
  p_adjust = "fdr",
  ci = 0.95)
pvalue_matrix<- cormat[,c(1,2,8)]
correlationcoeff_matrix<-cormat[,c(1,2,3)]
data<- as.matrix(pvalue_matrix)
data2 <- as.matrix(correlationcoeff_matrix)
png(height=1200,pointsize=13, width=1500,file="metscorr-excess.png")

corrplot(data2, method ="square", type = "upper", order = "hclust",tl.cex = 2,cl.cex=2, tl.offset=0.6, tl.col = "black", tl.srt = 90, p.mat = data,sig.level = .05, insig="blank",diag=FALSE)
dev.off()