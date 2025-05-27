{rm(list=ls())
  library(MASS)
  library(dlm)
  library(flexclust)
  library(Ckmeans.1d.dp)
  T<-500
  L<-100
  sigma<-10^(-3)*5
  p<-100
  q<-5
  n<-400
  K<-3
  tau<-1.345
  gamma<-3.7
  ALPHA<-matrix(0,T,n)
  BETA<-matrix(0,T,p)
  MU<-matrix(0,T,n)
  RIs<-matrix(0,T)
  mu0<-matrix(0,K,1)
  for (k in (1:K)){
    mu0[k]<-3*k
  }
  alpha0<-matrix(0,n,1)
  beta0<-matrix(c(2,2,2,3,3,matrix(0,1,p-q)),p,1)
  for (i in (1:n)){
    alpha0[i]<-mu0[ceiling(i*K/n)]
  }
  STEP<-matrix(0,1,T)
  TIMES<-matrix(0,T)
  RMSE<-matrix(0,T,3)
  MAE<-matrix(0,T,3)
  Q<-matrix(0,T)
  KKs<-matrix(0,T)
}
for (t in (1:T)){
  start_time <- Sys.time()
  set.seed(t)
  #数据生成
  epsilon<-matrix(rnorm(n,0,0.5),n,1)
  #epsilon<-0.5*matrix(rt(n,5),n,1)
  #epsilon<-0.95*matrix(rnorm(n,0,0.5),n,1)+0.05*matrix(rnorm(n,0,5),n,1)
  x<-mvrnorm(n,rep(0,p),diag(p))
  y<-alpha0+x%*%beta0+epsilon
  #初值
  D<-matrix(0,n*(n-1)/2,n)
  ll<-1
  for (i in (1:(n-1))){
    for (j in ((i+1):n)){
      D[ll,i]<-1
      D[ll,j]<--1
      ll<-ll+1
    }
  }
  lambda01<-0.02
  lambda02<-0.3
  A1<-diag(n)
  B1<-lambda01*t(D)%*%D
  A<-cbind(A1,x);
  B<-bdiag(B1,lambda01*diag(p))
  the<-solve(t(A)%*%A+B)%*%t(A)%*%y
  alp<-matrix(the[1:n],n,1)
  bet<-matrix(the[(n+1):(n+p)],p,1)
  for (l in (1:20)){
    alp<-solve(A1+lambda01*B1)%*%(y-x%*%bet)
    for (j in (1:p)){
      LL<-y-alp-x%*%bet+x[,j]*bet[j]
      a<-2*t(x[,j])%*%x[,j]
      b<-t(LL)%*%x[,j]/(t(x[,j])%*%x[,j])
      if(abs(b)<=lambda02*(1+1/a)){
        bet[j]<-sign(b)*(abs(b)-lambda02)*(abs(b)>=lambda02)
      }
      else if(abs(b)>gamma*lambda02){
        bet[j]<-b
      }
      else{
        bet[j]<-sign(b)*(abs(b)-gamma*lambda02/((gamma-1)*a))/(1-1/((gamma-1)*a))*(abs(b)>=gamma*lambda02/((gamma-1)*a))
      }
    }
  }
  B1<-lambda01*t(D)%*%D/n
  alp<-solve(diag(n)+B1)%*%(y-x%*%bet)
  #norm(alp-alpha0)/sqrt(n)
  
  KK<-c(2,3,4,5)
  lambda1<-seq(0.4,2,0.4)/4
  lambda2<-seq(0.15,0.75,0.15)
  BIC<-5000#初始非常大
  rho<-2
  alp0<-alp
  bet0<-bet
  for (kk1 in KK){
    for (kk2 in lambda1){
      for (kk3 in lambda2){
        l<-1
        alpp<-0
        bett<-0
        center<-Ckmeans.1d.dp(alp0,kk1)
        mu<-matrix(center$centers[center$cluster],n,1)
        s<-y-alp0-x%*%bet0
        v<-matrix(0,n,1)
        
        while(((l>=1)&&(l<=L))&&((norm(alp-alpp)/sqrt(n)+norm(alp-alpp)/sqrt(p))>sigma)){
          alpp<-alp
          bett<-bet
          #更新alpha
          alp<-(rho*(y-x%*%bet-s)+v+2*kk2*mu)/(rho+2*kk2)
          #更新mu
          center<-Ckmeans.1d.dp(alp,kk1)
          mu<-matrix(center$centers[center$cluster],n,1)
          #更新beta
          for (j in (1:p)){
            LL<-y-alp-x%*%bet+x[,j]*bet[j]-s+v/rho
            a<-rho*t(x[,j])%*%x[,j]
            b<-t(LL)%*%x[,j]/(t(x[,j])%*%x[,j])
            if(abs(b)<=kk3*(1+1/a)){
              bet[j]<-sign(b)*(abs(b)-kk3)*(abs(b)>=kk3)
            }
            else if(abs(b)>gamma*kk3){
              bet[j]<-b
            }
            else{
              bet[j]<-sign(b)*(abs(b)-gamma*kk3/((gamma-1)*a))/(1-1/((gamma-1)*a))*(abs(b)>=gamma*kk3/((gamma-1)*a))
            }
          }
          #更新s
          for (i in (1:n)){
            a<-y[i]-alp[i]-x[i,]%*%bet+v[i]/rho
            if(a>(tau+tau/(n*rho))){
              s[i]<-y[i]-alp[i]-x[i,]%*%bet+v[i]/rho-tau/(n*rho)
            }
            else if(a<(-(tau+tau/(n*rho)))){
              s[i]<-y[i]-alp[i]-x[i,]%*%bet+v[i]/rho+tau/(n*rho)
            }
            else{
              s[i]<-a/(1+1/(n*rho))
            }
          }
          #更新v
          v<-v+rho*(y-alp-x%*%bet-s)
          #更新l
          l<-l+1
        }
        
        bet[bet<0.005]<-0
        qq<-sum(abs(bet)>0)
        bic<-log(sum((y-alp-x%*%bet)^2/2*(abs(y-alp-x%*%bet)<=tau)+(tau*abs(y-alp-x%*%bet)-tau^2/2)*(abs(y-alp-x%*%bet)>tau))/n)+4*log(n*kk1+p)*(qq+kk1)*log(n)/n                              
        
        if(bic<=BIC){
          BIC<-bic
          STEP[t]<-l
          ALPHA[t,]<-alp
          BETA[t,]<-bet
          MU[t,]<-mu
          RIs[t]<-randIndex(table(alpha0,mu),correct = F)
          MAE[t,]<-c(mean(abs(alpha0-alp)),mean(abs(alpha0-mu)),mean(abs(beta0-bet)))
          RMSE[t,]<-c(norm(alpha0-alp)/sqrt(n),norm(alpha0-mu)/sqrt(n),norm(beta0-bet)/sqrt(p))
          Q[t,]<-qq
          KKs[t,]<-kk1
        }
      }
    }
  }
  print(t)
  end_time <- Sys.time()
  TIMES[t]<-as.numeric(end_time - start_time, units = "secs")
}
{
  MEAN_MAE_alpha<-round(mean(MAE[,1]), 4)
  MEAN_MAE_mu<-round(mean(MAE[,2]),4)
  MEAN_MAE_beta<-round(mean(MAE[,3]),4)
  VAR_MAE_alpha<-round(var(MAE[,1]),4)
  VAR_MAE_mu<-round(var(MAE[,2]),4)
  VAR_MAE_beta<-round(var(MAE[,3]),4)
  barQ<-mean(Q)
  titleQ<-median(Q)
  k2<-sum(KKs== 2)/T
  k3<-sum(KKs== 3)/T
  k4<-sum(KKs== 4)/T
  k5<-sum(KKs== 5)/T
  barK<-mean(KKs)
  titleK<-median(KKs)
  barRI<-mean(RIs)
  varRI<-var(RIs)
  Time<-mean(TIMES)
  result<-round(rbind(MEAN_MAE_alpha,VAR_MAE_alpha,MEAN_MAE_mu,VAR_MAE_mu,MEAN_MAE_beta,VAR_MAE_beta,barQ,titleQ,barK,titleK,barRI,varRI,Time),4)
  print(result)
}
