beta <- read.csv('timb.csv')[, -1]
p <- ncol(beta)
head(beta)
for(i in 1:nrow(beta)){
  temp <- beta[i, ]
  barplot(as.matrix(temp))
}

alpha <- read.csv('tima.csv')[, -1]
pp <- ncol(alpha)
head(alpha)
for(i in 1:nrow(alpha)){
  temp <- alpha[i, ]
  barplot(as.matrix(temp))
}
