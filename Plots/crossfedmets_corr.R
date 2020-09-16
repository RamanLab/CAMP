library(corrplot)
library(Hmisc)
setwd("D:/Paper_Draft/GitHub/Plots")
my_data <- read.csv("mets_corr.csv")
res3 <- rcorr(as.matrix(my_data))
res2 <- cor(as.matrix(my_data))
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
png(height=1200,pointsize=12, width=1500,file="metscorr.png")
flattenCorrMatrix <- function(cormat, pmat) {
ut <- upper.tri(cormat)+ data.frame(row = rownames(cormat)[row(cormat)[ut]],column = rownames(cormat)[col(cormat)[ut]],cor=(cormat)[ut],p = pmat[ut])
}
corrplot(res2, method ="square", type = "upper", order = "hclust",tl.cex = 2, tl.offset=0.6, tl.col = "black", tl.srt = 90, p.mat = res3$P, sig.level = .05, insig="blank",diag=FALSE)
dev.off()