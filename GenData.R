library(mvtnorm)

#####交互项处理两个函数
index.interac <- function(index.main, p, m){
  if (!all(index.main %in% 1:p)) stop("elements of index.main is wrong.")
  index.main <- sort(index.main)
  low <- as.vector(lower.tri(matrix(0, p, p)))
  low[low] <- 1:(p * (p - 1) / 2)
  Matr <- matrix(low, p, p, byrow = T)
  index.int <- as.vector(t(Matr[index.main, index.main]))
  index.int <- index.int[index.int != 0]
  index.int <- ((index.int - 1) * m + 1) : (index.int * m) 
  return(index.int)
}

index.main = function(index.interac, p){
  if (!all(index.interac %in% 1:(p*(p-1)/2))) stop("elements of index.interac is wrong.")
  index.interac = sort(index.interac)
  low = as.vector(lower.tri(matrix(0,p,p)))
  low[low] = 1:(p*(p-1)/2)
  Matr = matrix(low,p,p,byrow=T)
  index.int = unique(c(which(apply(Matr,1,function(x) any(x %in% index.interac))),
                       which(apply(Matr,2,function(x) any(x %in% index.interac)))))
  return(index.int)
}

####X交互
inter<-function(X){
  p = ncol(X)
  XX = NULL
  for (k in 1:(p-1)) XX = cbind(XX,X[,(k+1):p,drop=F]*X[,k])
  return(XX)
}

S1情况
b1=rep(c(2,1.5,1,1,0.8,0),c(1,1,1,1,1,p-5))
r1=rep(0,p*(p-1)/2)
r1[c(index.interac(1:3,p),index.interac(4:5,p))] = rep(c(2,1.5,1,0.5))
b2=rep(c(0.8,1,1,1.5,2,0),c(1,1,1,1,1,p-5))
r2=rep(0,p*(p-1)/2)
r2[c(index.interac(1:3,p),index.interac(4:5,p))] = rep(c(2,1,0.5,1.5))
b3=rep(c(2,1,0.8,1,1.5,0),c(1,1,1,1,1,p-5))
r3=rep(0,p*(p-1)/2)
r3[c(index.interac(1:3,p),index.interac(4:5,p))] = rep(c(1.5,1,0.5,2))

x1=rmvnorm(n1,mean=rep(0,p),sigma=diag(p)) 
xx1=inter(x1)
y1=x1%*%b1+xx1%*%r1+0.1*rnorm(n1,0,1)

x2=rmvnorm(n2,mean=rep(0,p),sigma=diag(p)) 
xx2=inter(x2)
y2=x2%*%b2+xx2%*%r2+0.1*rnorm(n2,0,1)

x3=rmvnorm(n3,mean=rep(0,p),sigma=diag(p)) 
xx3=inter(x3)
y3=x3%*%b3+xx3%*%r3+0.1*rnorm(n3,0,1)

S2、S3情况改b1,r1等的参数

S4情况
####AR协方差函数
cresigma<-function(p){
  sigma=matrix(nrow=p,ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      sigma[i,j]=0.5^(abs(i-j))
    }
  }
  return(sigma)
}

x1=rmvnorm(n1,mean=rep(0,p),sigma=cresigma(p)) 
xx1=inter(x1)
y1=x1%*%b1+xx1%*%r1+0.1*rnorm(n1,0,1)

x2=rmvnorm(n2,mean=rep(0,p),sigma=cresigma(p)) 
xx2=inter(x2)
y2=x2%*%b2+xx2%*%r2+0.1*rnorm(n2,0,1)

x3=rmvnorm(n3,mean=rep(0,p),sigma=cresigma(p)) 
xx3=inter(x3)
y3=x3%*%b3+xx3%*%r3+0.1*rnorm(n3,0,1)

S5情况(主效应X1-X5，交互X1X2,X1X3,X2X3,X1X6，不满足强分层)
b1=rep(c(2,1.5,1,1,0.8,0),c(1,1,1,1,1,p-5))
r1=rep(0,p*(p-1)/2)
r1[c(1,2,100,5)] = rep(c(2,1.5,1,0.5))
b2=rep(c(0.8,1,1,1.5,2,0),c(1,1,1,1,1,p-5))
r2=rep(0,p*(p-1)/2)
r2[c(1,2,100,5)] = rep(c(2,1,0.5,1.5))
b3=rep(c(2,1,0.8,1,1.5,0),c(1,1,1,1,1,p-5))
r3=rep(0,p*(p-1)/2)
r3[c(1,2,100,5)] = rep(c(1.5,1,0.5,2))

x1=rmvnorm(n1,mean=rep(0,p),sigma=diag(p)) 
xx1=inter(x1)
y1=x1%*%b1+xx1%*%r1+0.1*rnorm(n1,0,1)

x2=rmvnorm(n2,mean=rep(0,p),sigma=diag(p)) 
xx2=inter(x2)
y2=x2%*%b2+xx2%*%r2+0.1*rnorm(n2,0,1)

x3=rmvnorm(n3,mean=rep(0,p),sigma=diag(p)) 
xx3=inter(x3)
y3=x3%*%b3+xx3%*%r3+0.1*rnorm(n3,0,1)
