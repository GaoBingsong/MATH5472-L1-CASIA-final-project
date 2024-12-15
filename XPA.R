library(Matrix)
library(MASS)

set.seed(1)

# 随机逼近trace
estimate_trace <- function(A) {
  B=50
  trace_estimates <- rep(0, B)
  for (b in 1:B) {
    d <- rnorm(nrow(A))  # 生成随机向量
    trace_estimates[b] <- t(d) %*% A %*% d  # 计算二次型
  }
  mean(trace_estimates)  # 返回迹的估计值
}

estimate_effect <- function(X1,X2,Z1,Z2,y1,y2){
  
  Z <- bdiag(Z1,Z2)
  X <- bdiag(X1,X2)
  y <- c(y1,y2)
  
  M <- diag(n1+n2) - Z%*%solve(t(Z) %*% Z)%*%t(Z)
  M1 <- M[1:n1,1:n1]
  M2 <- M[(n1+1):(n1+n2),(n1+1):(n1+n2)]
  
  K1 <- X1%*%t(X1)
  K2 <- X2%*%t(X2)
  K12 <- X1%*%t(X2)
  
  tilde_K1 <- M1%*%K1%*%t(M1)
  tilde_K2 <- M2%*%K2%*%t(M2)
  tilde_K12 <- M1%*%K12%*%t(M2)
  
  tilde_y1 <- M1%*%y1
  tilde_y2 <- M2%*%y2
  
  #随机逼近计算S
  S <- matrix(rep(0,25),nrow=5)
  
  S[1,1]=estimate_trace(tilde_K1%*%tilde_K1)
  S[1,2]=sum(diag(tilde_K1))
  S[2,1]=S[1,2]
  S[2,2]=sum(diag(M1))
  
  S[3,3]=estimate_trace(tilde_K2%*%tilde_K2)
  S[3,4]=sum(diag(tilde_K2))
  S[4,3]=S[3,4]
  S[4,4]=sum(diag(M2))
  
  S[5,5]=estimate_trace(tilde_K12%*%t(tilde_K12))
  
  #计算q
  q <- rep(0,5)
  
  q[1]=t(tilde_y1)%*%tilde_K1%*%tilde_y1
  q[2]=t(tilde_y1)%*%tilde_y1
  
  q[3]=t(tilde_y2)%*%tilde_K2%*%tilde_y2
  q[4]=t(tilde_y2)%*%tilde_y2
  
  q[5]=t(tilde_y1)%*%tilde_K12%*%tilde_y2
  
  #得到参数估计
  para <- solve(S)%*%q
  sigma1_square <- para[1]
  sigma_epsilon_square <- para[2]
  sigma2_square <- para[3]
  sigma_xi_square <- para[4]
  delta <- para[5]
  
  #估计Omega
  
  Sigma_beta <- matrix(rep(0,4),nrow=2)
  Sigma_beta[1,1]=sigma1_square
  Sigma_beta[1,2]=delta
  Sigma_beta[2,1]=delta
  Sigma_beta[2,2]=sigma2_square
  
  Sigma_e <- bdiag(diag(sigma_epsilon_square,n1),diag(sigma_xi_square,n2))
  
  Omega <- X%*%(kronecker(Sigma_beta,diag(p)))%*%t(X) + Sigma_e
  Omega_inverse <- solve(Omega) 
  
  #估计fixed effect
  w <- solve(t(Z)%*%Omega_inverse%*%Z)%*%t(Z)%*%Omega_inverse%*%y
  w1 <- w[1:c1]
  w2 <- w[(c1+1):(c1+c2)]
  
  #估计random effect
  mu1 <- t(rbind(sigma1_square*X1,delta*X2))%*%Omega_inverse%*%rbind((y1-Z1%*%w1),(y2-Z2%*%w2))
  
  return(c(w,mu1))
}





n1=100
n1test=50
nall=n1+n1test

p=10
c1=2
c2=4

rou=0.8  # 0, 0.4
sigma1=1
sigma2=2
sigmae=0.5
sigmax=1

omega1=rep(1,c1)
omega2=rep(2,c2)

Rsquare=c()

for(i in 1:7){
  Rsquare_sum=c()
  for(r in 1:10){
    
    n2=n1*i
    
    Zall <- matrix(rnorm(nall*c1), nrow = nall)
    Z1 <- Zall[1:n1,]
    Z1test <- Zall[(n1+1):nall,]
    Z2 <- matrix(rnorm(n2*c2), nrow = n2)
    
    Gall=matrix(sample(0:1,nall*p,replace = T),nrow=nall)
    G1 <- Gall[1:n1,]
    G1test <- Gall[(n1+1):nall,]
    G2=matrix(sample(0:1,n2*p,replace = T),nrow=n2)
    
    X1 <- scale(G1)/sqrt(p)
    X2 <- scale(G2)/sqrt(p)
    
    g_sd=c()
    g_mean=c()
    for (j in 1:p) {
      g_sd=c(g_sd,sqrt(var(G1[,j])))
      g_mean=c(g_mean,mean(G1[,j]))
    }
    
    
    Sigmab <- matrix(c(sigma1^2, rou*sigma1*sigma2,
                       rou*sigma1*sigma2 , sigma2^2), nrow = 2, byrow = TRUE)
    
    beta <- mvrnorm(p, c(0,0), Sigmab)
    
    beta1 <- beta[,1]
    beta2 <- beta[,2]
    
    y1 <- Z1%*%omega1+X1%*%beta1+rnorm(n1,sd=sigmae)
    y1test <- Z1test%*%omega1+X1test%*%beta1+rnorm(n1test,sd=sigmae)
    y2 <- Z2%*%omega2+X2%*%beta2+rnorm(n2,sd=sigmax)
    
    effect <- estimate_effect(X1,X2,Z1,Z2,y1,y2)
    estimate_omega1<-effect[[1]][1:c1]
    estimate_beta1<-effect[[2]]
    tilde_mu1<-estimate_beta1/g_sd/sqrt(p)
    predict_y1<-Z1test%*%estimate_omega1+G1test%*%(tilde_mu1)-sum(g_mean*tilde_mu1)
    Rsquare_sum <- c(Rsquare_sum, 1-sum((y1test-predict_y1)^2)/sum((y1test-mean(y1test))^2))
  }
  Rsquare<-c(Rsquare,mean(Rsquare_sum))
}





