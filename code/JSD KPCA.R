
install.packages('philentropy')
library('philentropy')
XtX<-read.csv('matrix.csv',header=F)
M<-cbind(XtX[-1,-1])
M
res = eigen(100*M,symmetric = TRUE)




c1<-c(0,0.494347239, 0.678173478, 0.676409525, 0.679537971, 0.692390148,  0.67437596)
c2<-c(0.494347239, 0, 0.141903064, 0.180075889, 0.235436309, 0.540386165, 0.232293403)
c3<-c(0.678173478,0.141903064,0,0.029462989,0.104824512, 0.464524475, 0.105609224)
c4<-c(0.676409525, 0.180075889, 0.029462989,0, 0.036365132, 0.346765217, 0.036585096)
c5<-c(0.679537971, 0.235436309, 0.104824512, 0.036365132, 0, 0.251882519, 0.000545387)
c6<-c(0.692390148, 0.540386165,0.464524475, 0.346765217, 0.251882519,0, 0.251803192)
c7<-c(0.67437596, 0.232293403, 0.105609224, 0.036585096,0.000545387, 0.251803192,0)
K<- cbind(c1,c2,c3,c4,c5,c6,c7)
as.matrix(K)


# Kernel PCA
# Polynomial Kernel
# k(x,y) = t(x) %*% y + 1
#k1 = function (x,y) { (x[1] * y[1] + x[2] * y[2] + 1)^2 }
#K = matrix(0, ncol = N_total, nrow = N_total)
#for (i in 1:N_total) {
#  for (j in 1:N_total) { 
#    K[i,j] = JSD(rbind(X[i,], X[j,]))
#  }}
N_total=7
ones = 1/N_total* matrix(1, N_total, N_total)
K_norm = K - ones %*% K - K %*% ones + ones %*% K %*% ones
res = eigen(K_norm,symmetric = TRUE)
res

V = res$vectors
D = diag(res$values)

Y = K %*% V
plot(Y[,1], Y[,2],  main = "Kernel PCA (Poly)"
     , xlab="First component", ylab="Second component")



#P <- 1:10/sum(1:10)
#Q <- 20:29/sum(20:29)
#x <- rbind(P,Q)
#x
#JSD(x)

