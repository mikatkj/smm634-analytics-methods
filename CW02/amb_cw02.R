## load libraries ----
library(fpc)
library(mvtnorm)
library(dplyr)
library(ggplot2)
?clusterboot

## simulate data ----
generateGaussianData <- function(n, center, sigma, label) {
        data = rmvnorm(n, mean = center, sigma = sigma)
        data = data.frame(data)
        names(data) = c('dim1', 'dim2', 'dim3', 'dim4')
        data = data %>% mutate(class=factor(label))   
        data
}
symMatrix <- function(nrow, diag, off.diag) {
        n <- nrow
        d <- diag
        od <- off.diag
        matrix <- matrix(rep(d, n*n), nrow = n)
        matrix[upper.tri(matrix)] <- od
        matrix[lower.tri(matrix)] <- od
        matrix
}
# cluster: 1/2/3 - normal distribution
# parameters for multivariate normal distribution
n1 <- 1000
n2 <- 500
n3 <- 400
center1 <- c(0,0,0,0)
center2 <- c(3,2,3,2)
center3 <- c(-2,2,-2,2)
sigma1 <- symMatrix(nrow = 4, diag = 0.2, off.diag = 0)
sigma2 <- symMatrix(nrow = 4, diag = 0.5, off.diag = 0.4)
sigma3 <- symMatrix(nrow = 4, diag = 0.5, off.diag = -0.4)
# 2 clusters
data1 <- generateGaussianData(n1, center1, sigma1, 1)
data2 <- generateGaussianData(n2, center2, sigma2, 2)
data3 <- generateGaussianData(n3, center3, sigma3, 3)
# cluster: 4 - uniform distribution
# parameters for noise
n4 <- 20
dim1 <- runif(n = n4, min = -4, max = 5)
dim2 <- runif(n = n4, min = -2, max = 5)
dim3 <- dim1
dim4 <- dim2
data4 <- data.frame(dim1, dim2, dim3, dim4) %>% mutate(class = factor('noise'))
# generate one simulation data set
data <- bind_rows(data1, data2, data3, data4)

# produce "noise" dimension 5 and 6
dim5 <- rnorm(nrow(data), mean = 0, sd = 1)   # using normal~N(0,1)
dim6 <- rf(nrow(data), df1 = 10, df2 = 20)   # using F distribution
hist(dim5)
hist(dim6)
summary(dim5)
summary(dim6)
data <- data %>% mutate(dim5 = dim5)
data <- data %>% mutate(dim6 = dim6)
# move 'class' column to the last
data[,8] <- data[,5]
data[,5] <- list(NULL)
colnames(data)[7] <- 'class'

# have a look at data set
# plot(data[,1:4])
str(data)
# dimension 1 vs 2
data %>% ggplot(aes(x=dim1, y=dim2, shape=class)) + 
        geom_point() + 
        coord_fixed() +   # the same scale of x and y
        scale_shape_manual(values = c(1,2,3,4))   # symbols for points
# dimension 3 vs 4
data %>% ggplot(aes(x=dim3, y=dim4, shape=class)) + 
        geom_point() + 
        coord_fixed() +   # the same scale of x and y
        scale_shape_manual(values = c(1,2,3,4))   # symbols for points

# save to csv
write.csv(data, file = 'cw02_data.csv')
# import data
data <- read.csv("cw02_data.csv", row.names = 1, stringsAsFactors = TRUE)
str(data)
# MDS: first two dimensions
data.sc <- scale(data[,-7])
dist <- dist(data.sc)
mds <- cmdscale(dist)
mds
# plot mds
mds <- data.frame(mds)
mds <- mds %>% mutate(class = data$class)
mds %>% ggplot(aes(x=X1, y=X2, shape=class)) + 
        geom_point() + 
        labs(x = 'PCA1', y = 'PCA2') +
        coord_fixed() +   # the same scale of x and y
        scale_shape_manual(values = c(1,2,3,4))   # symbols for points

# clustering bootstrap
# 1) 4-means ----
cb.4means <- clusterboot(data[,-7], B = 100, bootmethod = c('boot','subset','noise','jitter'),
                         subtuning = 960,
                         noisetuning = c(0.05,2),
                         jittertuning = 0.1,
                         clustermethod = kmeansCBI,
                         krange = 4, seed = 111, count = FALSE)
cb.4means$bootmean  

# 2) 3-means ----
cb.3means <- clusterboot(data[,-7], B = 100, bootmethod = c('boot','subset','noise','jitter'),
                         subtuning = 960,
                         noisetuning = c(0.05,2),
                         jittertuning = 0.1,
                         clustermethod = kmeansCBI,
                         krange = 3, seed = 111, count = FALSE)

cb.3means$bootmean

# 3) hierarchical ----
# agglomerative hierarchical clustering with noise component
# noisecut: All clusters of size <= noisecut in the disthclustCBI/hclustCBI-partition are joined and declared as noise/outliers
cb.hier <- clusterboot(data[,-7],
                       B = 100, bootmethod = c('boot','subset','noise','jitter'),
                       subtuning = 960,
                       noisetuning = c(0.05,2),
                       jittertuning = 0.1,
                       clustermethod = hclustCBI,
                       k=3, cut='number', method='average', seed=111, count = FALSE)
cb.hier$bootmean

# 4) model-based clustering ----
?mclustBIC
cb.mclust <- clusterboot(data[,-7], B = 100, bootmethod = c('boot','subset','noise','jitter'),
                         subtuning = 960,
                         noisetuning = c(0.05,2),
                         jittertuning = 0.1,
                         clustermethod = noisemclustCBI,
                         G = NULL,
                         modelNames = c("VVE","VVV"), 
                         nnk = 5,  # clustering was first applied to data after 5 nearest neighbor denoising
                         summary.out = TRUE,
                         noisemethod = TRUE,   # If TRUE, the last cluster is regarded as "noise cluster"
                         seed = 111, count = FALSE)
cb.mclust$result$mclustsummary
