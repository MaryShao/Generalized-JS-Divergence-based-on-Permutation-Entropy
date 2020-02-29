install.packages('philentropy')
library('philentropy')

# 3 groups of sample, each size is 40
N_in_group = 40
N_group = 3
N_total = N_in_group * N_group
x=matrix(0, ncol = 2, nrow = N_total)
label = rep(c(1,2,3), rep(N_in_group,N_group))
# 4 graphs
layout(matrix(seq(4), 2,2))

# use polar coordinate to generate sample
# theta ~ UNIF(0, 2pi)
r = rep(c(1, 3, 6), rep(N_in_group,N_group))
theta = runif(N_total)*2*pi
x[,1] = (r  ) * cos(theta) + rnorm(N_total, sd = .2)
x[,2] = (r  ) * sin(theta) + rnorm(N_total, sd = .2)

plot(x[,1], x[,2], col = rainbow(3)[label], main = "Origin"
     , xlab="First dimension", ylab="Second dimension")

X = x
XtX = t(X) %*% X
XtX
#XtX<-read.csv('matrix.csv',header=F)
res = eigen(XtX)


V = res$vectors
V
D = diag(res$values)
D

# verify eigen decop
# sum(abs(XtX %*% V - V %*% (D)))
Y = X%*% V
Y
plot(Y[,1], Y[,2], col = rainbow(3)[label], main = "Traditional PCA" 
     , xlab="First component", ylab="Second component")


# Kernel PCA
# Polynomial Kernel
# k(x,y) = t(x) %*% y + 1
k1 = function (x,y) { (x[1] * y[1] + x[2] * y[2] + 1)^2 }
K = matrix(0, ncol = N_total, nrow = N_total)
for (i in 1:N_total) {
  for (j in 1:N_total) {
    K[i,j] =k1(X[i,], X[j,])
  }}
ones = 1/N_total* matrix(1, N_total, N_total)
K_norm = K - ones %*% K - K %*% ones + ones %*% K %*% ones
res = eigen(K_norm)

V = res$vectors
D = diag(res$values)

Y = K %*% V
plot(Y[,1], Y[,2], col = rainbow(3)[label], main = "Kernel PCA (Poly)"
     , xlab="First component", ylab="Second component")



# Gaussian Kernel
# k(x,y) = exp(-sum((x-y)^2)))
k2 = function (x,y) { dnorm(norm(matrix(x-y), type="F"))}
K = matrix(0, ncol = N_total, nrow = N_total)
for (i in 1:N_total) {
  for (j in 1:N_total) { 
    K[i,j] = k2(X[i,], X[j,])
  }}
ones = 1/N_total* matrix(1, N_total, N_total)
K_norm = K - ones %*% K - K %*% ones + ones %*% K %*% ones
res = eigen(K_norm)

V = res$vectors
D = diag(res$values)

Y = K %*% V
plot(Y[,1], Y[,2], col = rainbow(3)[label], main = "Kernel PCA (Gaussian)"
     , xlab="First component", ylab="Second component")


# Kernel PCA
# JSD Kernel
# k(x,y) = JSD

K = matrix(0, ncol = N_total, nrow = N_total)
for (i in 1:N_total) {
  for (j in 1:N_total) { 
    K[i,j] = JSD(rbind(abs(X[i,]), abs(X[j,])))
  }}

ones = 1/N_total* matrix(1, N_total, N_total)
K_norm = K - ones %*% K - K %*% ones + ones %*% K %*% ones
res = eigen(K_norm)

V = res$vectors
D = diag(res$values)

Y = K %*% V
plot(Y[,1], Y[,2], col = rainbow(3)[label], main = "Kernel PCA (JSD)"
     , xlab="First component", ylab="Second component")

###################
K = matrix(0, ncol = N_total, nrow = N_total)
for (i in 1:N_total) {
  for (j in 1:N_total) { 
    M<-rbind(X[i,], X[j,])
    max(M)
    X1<-(X[i,]-min(M))/(max(M)-min(M))
    X2<-(X[j,]-min(M))/(max(M)-min(M))
    K[i,j] = JSD(rbind(X1,X2))
  }}

ones = 1/N_total* matrix(1, N_total, N_total)
K_norm = K - ones %*% K - K %*% ones + ones %*% K %*% ones
res = eigen(K_norm)

V = res$vectors
D = diag(res$values)

Y = K %*% V
plot(Y[,1], Y[,2], col = rainbow(3)[label], main = "Kernel PCA (JSD_unit)"
     , xlab="First component", ylab="Second component")
