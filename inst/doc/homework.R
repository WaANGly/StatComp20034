## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(StatComp20034)
library(MASS)
library(Ball)
library(energy)
library(DAAG)
library(RANN)
library(boot)
library(bootstrap)
library(inline)

## ----echo=FALSE---------------------------------------------------------------
x<-rnorm(20)
y<-rpois(20,10)
plot(x,y)


## -----------------------------------------------------------------------------
n<-100
u<-runif(n)
x<-2*(1-u)^{-1/2}
hist(x, prob = TRUE, main = expression(f(x)==8*x^{-3}))
y<-seq(2,10000,.01)
lines(y,8*y^{-3})

## -----------------------------------------------------------------------------
n<-10000
u1<-runif(n,-1,1)
u2<-runif(n,-1,1)
u3<-runif(n,-1,1)
data<-data.frame(u1,u2,u3)
random<-function(data,c1,c2,c3){
  d<-ifelse(abs(data[c3])>=max(abs(data[c2]),abs(data[c3])>=abs(data[c1])),data[c2],data[c3])
return(d)}
u<-apply(data,1,random)
hist(u,prob=TRUE,main=NULL)


## -----------------------------------------------------------------------------
u<-runif(1000)
y<-2*(1-u)^{-1/4}-1
hist(y,prob=TRUE,main=expression(f(x)==64*(2+x)^{-5}))
y=seq(0,10000,0.1)
lines(y,64*(2+y)^{-5})

## -----------------------------------------------------------------------------
set.seed(333)
t<-runif(1e6,0,pi/3)
theta.hat<-mean(sin(t))*pi/3
f<-c(theta.hat,cos(0)-cos(pi/3))
f

## -----------------------------------------------------------------------------
set.seed(1234)
n<-1e5
u<-runif(n)
U<-exp(u)
x<- (exp(u)+exp(1-u))/2
theta1<-mean(U)# the simple Monte Carlo method
theta2<-mean(x)#the antithetic variate approach
theta<-exp(1)-1#the real value
cat("the simple MC:",theta1,"\nthe antithetic variate approach:",theta2,"\nthe real value:",theta,"\n")

v1<-var(U)
v2<-var(x)
reduction<-(v1-v2)/v1# the percent reduction in variance 
cat("empirical percent reduction in variate:",reduction,"\n")

tv1<--exp(2)+4*exp(1)-3#theoretical variance using simple MC
tv2<--3*exp(2)+10*exp(1)-5#theoretical variance using the antithetic variate approach
d<-(tv1-2*tv2)/tv1
cat("theoretical percent reduction in variance:",d,"")

## -----------------------------------------------------------------------------
set.seed(333)
t<-runif(1e6,0,pi/3)
theta.hat<-mean(sin(t))*pi/3
f<-c(theta.hat,cos(0)-cos(pi/3))
f

## -----------------------------------------------------------------------------
set.seed(1234)
n<-1e5
u<-runif(n)
U<-exp(u)
x<- (exp(u)+exp(1-u))/2
theta1<-mean(U)# the simple Monte Carlo method
theta2<-mean(x)#the antithetic variate approach
theta<-exp(1)-1#the real value
cat("the simple MC:",theta1,"\nthe antithetic variate approach:",theta2,"\nthe real value:",theta,"\n")

v1<-var(U)
v2<-var(x)
reduction<-(v1-v2)/v1# the percent reduction in variance 
cat("empirical percent reduction in variate:",reduction,"\n")

tv1<--exp(2)+4*exp(1)-3#theoretical variance using simple MC
tv2<--3*exp(2)+10*exp(1)-5#theoretical variance using the antithetic variate approach
d<-(tv1-2*tv2)/tv1
cat("theoretical percent reduction in variance:",d,"")

## -----------------------------------------------------------------------------
sk <- function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
 }#computes the sample skewness coeff.

n<-30
m<-1000
alpha<-seq(0,100,1)
p<- numeric(length(alpha))
cv<-qnorm(0.975,0,sqrt(6*(n-2)/((n+1)*(n+3))))

for (j in 1:length(alpha)) { #for each alpha
  a<-alpha[j]
  sktests <- numeric(m)
for (i in 1:m) { #for each replicate
  x <- rbeta(n, shape1=a, shape2=a)
  sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
p[j] <- mean(sktests)
}
plot(alpha, p, type = "l",
xlab = "alpha", ylim = c(0,0.15))
abline(h=0.05,lty=3)

## -----------------------------------------------------------------------------
sk <- function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
 }#computes the sample skewness coeff.

n<-30
m<-1000
v<-seq(1,100,1)
p<- numeric(length(v))
cv<-qnorm(0.975,0,sqrt(6*(n-2)/((n+1)*(n+3))))

for (j in 1:length(v)) { #for each v
  a<-v[j]
  sktests <- numeric(m)
for (i in 1:m) { #for each replicate
  x <- rt(n,a)
  sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
p[j] <- mean(sktests)
}
plot(v, p, type = "l",xlab = "v",ylim=c(0,1))
abline(h=0.05,lty=3)


## -----------------------------------------------------------------------------
sigma1 <- 1
sigma2 <- 1.5
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
power <- mean(replicate(m, expr={
x <- rnorm(20, 0, sigma1)
y <- rnorm(20, 0, sigma2)
count5test(x, y)
}))

print(power)

## -----------------------------------------------------------------------------
alpha<-0.055
n1<-20
n2<-20
d1<-qf(alpha/2,n1-1,n2-1)
d2<-qf(1-alpha/2,n1-1,n2-1)
m<-1000
b<-numeric(m)
for(i in 1:m){
  x <- rnorm(n1, 0, 1)
  y <- rnorm(n2, 0, 1.5)
  s<-var(x)/var(y)
  b[i]<-as.integer((s>d1)&&(s<d2))
}
mean(b)

## -----------------------------------------------------------------------------
alpha <- 0.055
n <- c(20,200,2000,20000)
m <- 1000 
cf <- f <-b1<-b2<-  vector()
# estimate power
for (j in 1:4) {
  a<-n[j]
for(i in 1:m){
  x <- rnorm(a, 0, 1)
  y <- rnorm(a, 0, 1.5)
  outx <- sum(x > max(y)) + sum(x < min(y))
  outy <- sum(y > max(x)) + sum(y < min(x))
  b1[i]<-as.integer(max(c(outx, outy)) > 5)#the Count Five test
  s<-var(x)/var(y)
  b2[i]<-as.integer((s>d1)&&(s<d2))#  the F test
}
  cf[j]=mean(b1)
  f[j]=mean(b2)
  }
d = data.frame(cf,f,row.names = c("n=20","n=200","n=2000","n=20000"))
knitr::kable(d)

## -----------------------------------------------------------------------------

sigma<-matrix(c(5,1,1,4),2,2)
mu<-c(1,4)
d<-2
c1<-qchisq(0.025,d*(d+1)*(d+2)/6)
c2<-qchisq(0.975,d*(d+1)*(d+2)/6)
m<-1000
n<- c(10,20,30,50,100,500)#样本量
p<-numeric(length(n))

for(i in 1:length(n)){
  
b<-replicate(m,expr={
               x<-mvrnorm(n[i],mu,sigma)
               sigmahat<-var(x)
               y<-matrix(c(mean(x[,1]),mean(x[,2])),n[i],2,byrow=TRUE)
               xc<-x-y
               a<-(xc%*%solve(sigmahat)%*%t(xc))^3
               mean(a)
             })
p[i]<-mean(as.integer((n[i]*b/6>c2)|(n[i]*b/6<c1)))
           } 
g<-data.frame(p,row.names=c("n=10","n=20","n=30","n=50","n=100","n=500"))
knitr::kable(g)

## -----------------------------------------------------------------------------
sigma1<-matrix(c(5,1,1,4),2,2)
sigma2<-matrix(c(8,2,2,1),2,2)
mu<-c(0,0)
d<-2
c1<-qchisq(0.025,d*(d+1)*(d+2)/6)
c2<-qchisq(0.975,d*(d+1)*(d+2)/6)
m<-1000
n<-30
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
for (j in 1:N) { 
e <- epsilon[j]
b<-replicate(m,expr={
  a <- sum(sample(c(0,1), replace = TRUE,size = n, prob = c(1-e, e)))
           if(a==0){x<-mvrnorm(n,mu,sigma1)}
  else if(a==n){x<-mvrnorm(n,mu,sigma2)}
else                {x1 <- mvrnorm(n-a,mu,sigma1)
                x2 <- mvrnorm(a,mu,sigma2)
                x <- rbind(x1,x2)
}
               sigmahat<-var(x)
               y<-matrix(c(mean(x[,1]),mean(x[,2])),n,2,byrow=TRUE)
               xc<-x-y
               bhat<-(xc%*%solve(sigmahat)%*%t(xc))^3
               mean(bhat)
             })
 pwr[j] <- mean(as.integer((n*b/6<c1)|(n*b/6>c2)))}

plot(epsilon, pwr, type="b",xlab=bquote(epsilon),ylim=c(0,0.5))
abline(h=0.05, lty=3)


## -----------------------------------------------------------------------------
n<-nrow(law)
thetahat<-cor(law$LSAT, law$GPA)
thetaj<-LSAT1<-GPA1 <- numeric(n)
for (i in 1:n){
  thetaj[i]<-cor(law$LSAT[-i],law$GPA[-i])
}
bias <- (n - 1) * (mean(thetaj) - thetahat)
cat("jackknife estimate of bias:",bias,"\n") #jackknife estimate of bias

se <- sqrt((n-1) *mean((thetaj - mean(thetaj))^2))
cat("the standard error of the correlation:",se,"")

## -----------------------------------------------------------------------------
theta.boot <- function(dat,ind) {
#function to compute the statistic
  y <- dat[ind,1]
  mean(y) 
}
data(aircondit,package="boot")
boot.obj <- boot(aircondit, statistic =theta.boot , R = 2000)
print(boot.ci(boot.obj,type = c("basic","norm","perc","bca")))
#calculations for bootstrap confidence intervals
alpha <- c(.025, .975)

#normal
print(boot.obj$t0 + qnorm(alpha )* sd(boot.obj$t))


#basic
print(2*boot.obj$t0 -quantile(boot.obj$t, rev(alpha), type=1))

#percentile
print(quantile(boot.obj$t, alpha, type=6))


## -----------------------------------------------------------------------------
 boot.BCa <-function(x, th0, th,  conf = .95) {
# bootstrap with BCa bootstrap confidence interval
# th0 is the observed statistic
# th is the vector of bootstrap replicates
 x <- as.matrix(x)
 n <- nrow(x) #observations in rows
 N <- 1:n
 alpha <- (1 + c(-conf, conf))/2
 zalpha <- qnorm(alpha)
# the bias correction factor
 z0 <- qnorm(sum(th < th0) / length(th))
# the acceleration factor (jackknife est.)
 th.jack <- numeric(n)
 for (i in 1:n) {
 
 th.jack[i] <- mean(x[-i, ])
}
 L <- mean(th.jack) - th.jack
 a <- sum(L^3)/(6 * sum(L^2)^1.5)
# BCa conf. limits
 adj.alpha <- pnorm(z0 + (z0+zalpha)/(1-a*(z0+zalpha)))
 limits <- quantile(th, adj.alpha, type=6)
 return(list("est"=th0, "BCa"=limits))
 }

data(aircondit,package="boot")
B <- 2000
n<-nrow(aircondit)
theta.b <- numeric(B)
theta.hat <- mean(aircondit[,1])

#bootstrap
for (b in 1:B) {
i <- sample(1:n, size = n, replace = TRUE)
y <- aircondit[i,]
theta.b[b] <- mean(y) }

#compute the BCa interval
boot.BCa(aircondit, th0 = theta.hat, th = theta.b)


## -----------------------------------------------------------------------------
data(scor,package="bootstrap")
sigma<-cov(scor)# MLE of covariance matrix
eig<-eigen(sigma)
e<-sum(eig$values)
thetahat<-(eig$values[1])/e

for (j in 1:nrow(scor)){
  c<-cov(scor[-j,])
  thetaj[j]<-(eigen(c)$values[1])/sum(eigen(c)$values)
}
bias <- (n - 1) * (mean(thetaj) - thetahat)
cat("jackknife estimate of bias:",bias,"\n") #jackknife estimate of bias

se <- sqrt((n-1) *mean((thetaj - mean(thetaj))^2))
cat("the standard error of the correlation:",se,"")

## -----------------------------------------------------------------------------
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1<-e2<-e3<-e4<-e11<-e22<-e33<-e44<-matrix(0,n,n-1)

for (k in 1:(n-1)) {
  x<-chemical[-k]
  y<-magnetic[-k]
  for(j in k:(n-1)){
  x1<-x[-j]
  y1<-y[-j]
  
  J1 <- lm(y1 ~ x1)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  yhat11<-J1$coef[1] + J1$coef[2] * x[j]
  e1[k,j] <- magnetic[k] - yhat1
  e11[k,j]<-x[j]-yhat11
  
  J2 <- lm(y1 ~ x1 + I(x1^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +J2$coef[3] * chemical[k]^2
  yhat22 <- J2$coef[1] + J2$coef[2] * x[j] +J2$coef[3] * x[j]^2
  e2[k,j] <- magnetic[k] - yhat2
  e22[k,j]<-x[j]-yhat22
  
  J3 <- lm(log(y1) ~ x1)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  logyhat33 <- J3$coef[1] + J3$coef[2] * x[j]
  yhat3 <- exp(logyhat3)
  yhat33 <- exp(logyhat33)
  e3[k,j] <-chemical[k]-yhat3
  e33[k,j]<-x[j]-yhat33
  
  J4 <- lm(log(y1) ~ log(x1))
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    logyhat44 <- J4$coef[1] + J4$coef[2] * log(x[j])
  yhat4 <- exp(logyhat4)
  yhat44 <-exp(logyhat44)
  e4[k,j]<- magnetic[k] - yhat4
  e44[k,j]<-x[j]-yhat44
  }
  }
cat("the prediction error of J1:",mean(e1^2+e11^2),"\nthe prediction error of J2:",mean(e2^2+e22^2),"\nthe prediction error of J3:",mean(e3^2+e33^2),"\nthe prediction error of J4:",mean(e4^2+e44^2),"\n")


## -----------------------------------------------------------------------------
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
set.seed(1)
R<-1000
n1<-20
n2<-30
x<-rnorm(n1,0,1)
y<-rnorm(n2,0,1)
c<-numeric(R)
for(i in 1:R){
  a<-sample(y,size=5,replace=FALSE)
  xx<-c(x,a)
  yy<-setdiff(y,a)
  c[i]<-count5test(xx,yy)
}
p<-mean(c)
p

## -----------------------------------------------------------------------------
#方差不等的情况下
set.seed(1)
R<-1000
n1<-20
n2<-30
x1<-rnorm(n1,0,1)
y1<-rnorm(n2,0,5)
c1<-numeric(R)
for(i in 1:R){
  a<-sample(y1,size=5,replace=FALSE)
  xx1<-c(x1,a)
  yy1<-setdiff(y1,a)
  c1[i]<-count5test(xx1,yy1)
}
p1<-mean(c1)
p1

## -----------------------------------------------------------------------------
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) 
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
  

m <- 1000; k<-3
n1 <- n2 <- 50; R<-99
N<-c(n1,n2)
p1<-p2<-p3<-p4<-p5<-p6<-matrix(NA,m,3)

#Unequal variances and equal expectations
for(i in 1:m){
  x1 <- matrix(rnorm(100,mean=0.5,sd=1), ncol=2)
  y1 <- matrix(rnorm(100,mean=0.5,sd=1.5),ncol=2)
  z1 <- rbind(x1, y1)
  p1[i,1]<-eqdist.nn(z1,N,k)$p.value
  p1[i,2]<-eqdist.etest(z1,sizes=N,R=R)$p.value
  p1[i,3]<-bd.test(x=x1,y=y1,R=99,seed=i*12345)$p.value
}
pow1 <- colMeans(p1<0.1)
pow1

## -----------------------------------------------------------------------------
#Unequal variances and unequal expectations
for(i in 1:m){
  x2 <- matrix(rnorm(100,mean=0.5,sd=1.2), ncol=2)
  y2 <- matrix(rnorm(100,mean=1,sd=1.5),ncol=2)
  z2 <- rbind(x2, y2)
  p2[i,1]<-eqdist.nn(z2,N,k)$p.value
  p2[i,2]<-eqdist.etest(z2,sizes=N,R=R)$p.value
  p2[i,3]<-bd.test(x=x2,y=y2,R=99,seed=i*12345)$p.value
}
pow2 <- colMeans(p2<0.1)
pow2

## -----------------------------------------------------------------------------
#Non-normal distributions: t distribution with 1 df，bimodel distribution
for(i in 1:m){
  x3 <- matrix(rt(100,df=1), ncol=2)
  y<-numeric()
  for(j in 1:100){
      r<-sample(c(0,1),1,replace=TRUE,prob=c(0.3,0.7))
      y[j]<-ifelse(r==0,rnorm(1,0.5,2),rnorm(1,1,3))
  }
  y3<-matrix(y,ncol=2)
  z3 <- rbind(x3, y3)
  p3[i,1]<-eqdist.nn(z3,N,k)$p.value
  p3[i,2]<-eqdist.etest(z3,sizes=N,R=R)$p.value
  p3[i,3]<-bd.test(x=x3,y=y3,R=99,seed=i*12345)$p.value
}
pow3 <- colMeans(p3<0.1)
pow3

## -----------------------------------------------------------------------------
#t distribution 
for(i in 1:m){
  x4 <- matrix(rt(100,df=1), ncol=2)
  y4<-matrix(rt(100,df=5),ncol=2)
  z4 <- rbind(x4, y4)
  p4[i,1]<-eqdist.nn(z4,N,k)$p.value
  p4[i,2]<-eqdist.etest(z4,sizes=N,R=R)$p.value
  p4[i,3]<-bd.test(x=x4,y=y4,R=99,seed=i*12345)$p.value
}
pow4 <- colMeans(p4<0.1)
pow4


## -----------------------------------------------------------------------------
#unbalanced samples
n3<-25;n4<-250;n5 <- n3+n4;N1<-c(n3,n4)
for(i in 1:m){
  x6 <- matrix(rnorm(50,mean=0.5,sd=1), ncol=2)
  y6 <- matrix(rnorm(500,mean=0.5,sd=1.5),ncol=2)
  z6 <- rbind(x6, y6)
  p6[i,1]<-eqdist.nn(z6,N1,k)$p.value
  p6[i,2]<-eqdist.etest(z6,sizes=N1,R=R)$p.value
  p6[i,3]<-bd.test(x=x6,y=y6,R=99,seed=i*12345)$p.value
}
pow6 <- colMeans(p6<0.1)
pow6


## -----------------------------------------------------------------------------
wl <- function(sigma, x0, N) {
 x <- numeric(N)
 x[1] <- x0
 u <- runif(N)
 k <- 0
 for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    r<-exp(abs(x[i-1])-abs(y))
    if (u[i] <= r)
     x[i] <- y else {
      x[i] <- x[i-1]
      k <- k + 1
  }
 }
  return(list(x=x, k=k))
}

N <- 1000
sigma <- c(.5, 2, 10, 16)
x0 <- 10
wl1 <- wl( sigma[1], x0, N)
wl2 <- wl( sigma[2], x0, N)
wl3 <- wl( sigma[3], x0, N)
wl4 <- wl( sigma[4], x0, N)
#number of candidate points accepted
print(c(wl1$k, wl2$k, wl3$k, wl4$k))
accept<-c((N-1-wl1$k)/N,(N-1-wl2$k)/N,(N-1-wl3$k)/N,(N-1-wl4$k)/N)
m<-matrix(c(sigma,accept),ncol=2,dimnames=list(c("wl1","wl2","wl3","wl4"),c("variance","acceptance rate")))
m

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

LA.chain <- function(sigma, N, X1) {
#generates a Metropolis chain for standard Laplace distribution
#with Normal(X[t], sigma) proposal distribution
#and starting value X1
  x <- rep(0, N)
  x[1] <- X1
  u <- runif(N)
  for (i in 2:N) {
     xt <- x[i-1]
     y <- rnorm(1, xt, sigma) #candidate point
     r1 <- exp(-abs(y)) * dnorm(xt, y, sigma)/2
     r2 <- exp(-abs(xt)) * dnorm(y, xt, sigma)/2
     r <- r1 / r2
     if (u[i] <= r) x[i] <- y else
     x[i] <- xt
}
return(x)
}

sigma <- 5 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 15000 #length of chains
b <- 1000 #burn-in length
#choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
X[i, ] <- LA.chain(sigma, n, x0[i])

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))






## -----------------------------------------------------------------------------
root=function(k){
s1=function(a){
  1-pt(sqrt((k-1)*a^2/(k-a^2)),k-1)
}
s2=function(a){
  1-pt(sqrt(k*a^2/(k+1-a^2)),k)
}
f<-function(a){
  s1(a)-s2(a)
}
return(uniroot(f,interval = c(1e-6,sqrt(k)-1e-6))$root)
}
r = sapply(c(4:25, 100, 500, 1000), function (k) {
  root(k)
  })
r



## ----echo=FALSE---------------------------------------------------------------
r<-matrix(c('Frequency','p^2','q^2','r^2','2pr','2qr','2pq','1','Count','nAA','nBB','nOO','nAO','nBO','nAB','n'),nrow=2,byrow=TRUE,dimnames = NULL)
knitr::kable(r,row.names = NA,col.names =c('Genotype','AA','BB','OO','AO','BO','AB','Sum') )


## -----------------------------------------------------------------------------

na<-444
nb<-132
noo<-361
nab<-63
n<-na+nb+nab+noo
lg<-numeric(11)
p<-0.6
q<-0.3#设置初始值
   for(i in 1:15){
   nao<-na*(1-p/(2-p-2*q))
   nbo<-nb*(1-q/(2-2*p-q))
   lg[i]<-na*log(p^2+2*p*(1-p-q))+nb*log(q^2+2*q*(1-p-q))+noo*log((1-p-q)^2)+nab*log(2*p*q)#the corresponding log-maximum likelihood values (for observed data)
   p<-(nab+2*na-nao)/(2*n)
   q<-(nab+2*nb-nbo)/(2*n)
   }
cat("p=",p,"\nq=",q,"\n")
print(lg)

## -----------------------------------------------------------------------------
 formulas <- list(
       mpg ~ disp,
       mpg ~ I(1 / disp),
       mpg ~ disp + wt,
       mpg ~ I(1 / disp) + wt
       )

# use for loops
L <- vector("list", length(formulas))
for (i in seq_along(formulas)) {
  L[[i]] <- lm(formulas[[i]],data=mtcars)
}
print(L)

#use lapply()
lapply(formulas,lm,data=mtcars)

## -----------------------------------------------------------------------------

 trials <- replicate(100,
      t.test(rpois(10, 10), rpois(7, 10)),
      simplify = FALSE
      )
#use sapply() and the anonymous function
sapply(trials,function(x){x$p.value})

# get rid of the anonymous function
sapply(trials,"[[","p.value")


## -----------------------------------------------------------------------------

L <- list(mtcars,trees)

lapply(L, function(x) vapply(x, mean, numeric(1)))

Mapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE){return(simplify2array(out))}
  out
}

Mapply(L, mean, numeric(1))


## -----------------------------------------------------------------------------
rwr <- function(sigma, x0, N){
 x = numeric(N)
 x[1] = x0
 u = runif(N)
 k = 0
 for (i in 2:N) {
  y = rnorm(1, x[i-1], sigma)
  if (u[i] <= exp(abs(x[i-1])-abs(y))) x[i] = y 
  else {
  x[i] = x[i-1]
  k = k+1
  }
 }
 return(list(x = x, k = k))
}


