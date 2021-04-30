## Question 1 ----
library(rospca)
?dataGen

# Q.1a ----
N <- 100 # number of observations
P <- 9   # number of variables
var <- 3 # number of variables for each group
corr <- c(0.7, 0.9, 0.8) # inner correlation for each group
sd <- c(20, 5, 10)   # standard deviation for each group

simulate1 <- dataGen(m = 1, n = N, p = P, bLength = var, a = corr, SD = sd)
data1 <- data.frame(simulate1[[1]])

# Q.1b ----
# check the standard deviations
apply(data1, 2, sd)

# check the inner correlations
library(multicon)
?inner.outer
# create a list indicating the items belonging to each group
list1 <- list(data1[,1:3], data1[,4:6], data1[,7:9]) 
inner.outer(list1)

# check the correlations
cor(data1[, 1:3])
cor(data1[, 4:6])
cor(data1[, 7:9])
# correlation plots
library(corrplot)
corrplot(cor(data1), type="upper", order="hclust")
# par(mfrow=c(2,2))
# corrplot(cor(data[, 1:3]), type="upper", order="hclust")
# corrplot(cor(data[, 4:6]), type="upper", order="hclust")
# corrplot(cor(data[, 7:9]), type="upper", order="hclust")

# Q.1c ----
simulate2 <- dataGen(m = 1, n = N, p = P, bLength = var, 
                     a = c(0, 0, 0), SD = sd)
data2 <- data.frame(simulate2[[1]])
head(data2)

apply(data2, 2, sd)
list2 <- list(data2[,1:3], data2[,4:6], data2[,7:9])
inner.outer(list2)

# Q.1d ----
# (1) on the basis of the sample covariance matrix (not-scaled)
var(data1)
eig.cov1 <- eigen(var(data1))
eig.cov1$values
eig.cov1$vectors
# change the divisor (var() function uses (N-1) as divisor)
sqrt(((N-1)/N) * eig.cov1$values)
# (2) on the basis of the sample correlation matrix (scaled)
cor(data1)
eig.cor1 <- eigen(cor(data1))
eig.cor1$values
eig.cor1$vectors
sqrt(((N-1)/N) * eig.cor1$values)

# (3) use princomp (and not-scaled) -> compared with (1)
pca.cov1 <- princomp(data1)
summary(pca.cov1)    # equivalent to: sqrt(((N-1)/N) * eig.cov1$values)
loadings(pca.cov1)    # equivalent to：eig.cov1$vectors
# (4) use princomp (and scaled) -> compared with (2)
pca.cor1 <- princomp(data1, cor = TRUE)
summary(pca.cor1)   # close to: sqrt(((N-1)/N) * eig.cor1$values)
loadings(pca.cor1)   # equivalent to: eig.cor1$vectors

# (5) use prcomp (and not-scaled)
pr.cov1 <- prcomp(data1)
summary(pr.cov1)   # close to: summary(pca.cov1)
pr.cov1$rotation   # equivalent to: loadings(pca.cov1) 
# (6) use prcomp (and scaled)
pr.cor1 <- prcomp(data1, scale = TRUE)
summary(pr.cor1)   # 同：summary(pca.cor1)  
pr.cor1$rotation   # 同：loadings(pca.cor1) 

# (7) compare difference between scaling or not on biplot
biplot(pr.cov1, scale = 0)   # dataset1, not scaled
biplot(pr.cor1, scale = 0)   # dataset1, scaled
# dataset2
pr.cov2 <- prcomp(data2)
pr.cor2 <- prcomp(data2, scale = TRUE)
biplot(pr.cov2, scale = 0)   # dataset2, not scaled
biplot(pr.cor2, scale = 0)   # dataset2, scaled

# PCA on data1 ----
pr.out1 <- prcomp(data1, scale = TRUE)
summary(pr.out1)
# 'rotation' matrix - the PC loadings; each column is a PC loading vector
pr.out1$rotation   # loadings
# expect 4 columns, because min(n-1, p)=9

# mean and sd of original variables（NOT PC's） 
pr.out1$center   # apply(data1, 2, mean)
pr.out1$scale   # apply(data1, 2, sd)
# variance of pc
pr.out1$sdev
# 'x' - all of the pc
# the matrix X -> its columns as the PC score vectors
# the k-th column is the k-th principal component score vector
dim(pr.out1$x)
pr.out1$x[,1:2]   # x scores for PC1 PC2

# plot the first 2 PC
# scale=0: arrows are scaled to represent the loadings
biplot(pr.out1, scale = 0)

# PVE - proportion of variance explained by each PC
# the variance explained by each PC: obtained by squaring SE
pr.var1 <- pr.out1$sdev^2
pr.var1
pve1 <- pr.var1/sum(pr.var1)
pve1
# the 1st PC explains 62.0% of the variance in the data
# the 2nd PC explains 24.7% of ...

# plot the PVE explained by each PC, and the cumulative PVE
par(mfrow = c(1, 2))
plot(pve1, xlab = 'Principal Component', ylab = 'Proportion of Variance Explained',
     ylim = c(0, 1), type = 'b')
plot(cumsum(pve1), xlab = 'Principal Component', 
     ylab = 'Cumulative Proportion of Variance Explained',
     ylim = c(0, 1), type = 'b')


# PCA on data2 ----
pr.out2 <- prcomp(data2, scale = TRUE)
biplot(pr.out2, scale = 0)

# PVE
pr.var2 <- pr.out2$sdev^2
pr.var2
pve2 <- pr.var2/sum(pr.var2)
pve2

par(mfrow = c(1, 2))
plot(pve2, xlab = 'Principal Component', ylab = 'Proportion of Variance Explained',
     ylim = c(0, 1), type = 'b')
plot(cumsum(pve2), xlab = 'Principal Component', 
     ylab = 'Cumulative Proportion of Variance Explained',
     ylim = c(0, 1), type = 'b')
