 rm(list=ls())

dyn.load("~/Documents/spatial_temporal/c_code/m2.dll")
require(MASS)
require(doMC)
registerDoMC(cores=6)
idx<-as.numeric(Sys.getenv("LSB_JOBINDEX"))
filname<-paste("m2_s",idx,".RData", sep="")

###likelihood###
LL=function(b0,b1,b2,b3,sig,W,y,s,x,x2,year,region_code){
    m=b0+b1*x+W[region_code,1]+b2*x2+b3*(year-1)+W[region_code,2]*(year-1)
    #lam=exp(-m/sig)
    #a=1/sig
    #ll=sum(s*log(lam*a*y^(a-1))-lam*y^a)
    
    ll=sum(dpois(y,m,log=T))
  }
  

  nsim=500
  p=2
  cancer_weight=read.table("~/Documents/spatial_temporal/paper_in_preparation/weight.txt",header=T)[,4]
  #adj_max=t(t(c(1,1,1)))%*%c(1,1,1)-diag(1,3,3)
  adj_max=as.matrix(read.csv("~/Documents/spatial_temporal/paper_in_preparation/adjacent_matrix_PA.csv"))
  region=nrow(adj_max)
  total_patient=5/min(cancer_weight)
  total_region_pat=as.integer(total_patient*cancer_weight)
  each_year_pat=as.integer(total_region_pat/5)
  region_code=rep(1:region,each_year_pat*5)
  year=unlist(sapply(each_year_pat,function(x) rep(1:5,each=x)))
  year1=(year==2)
  year2=(year==3)
  year3=(year==4)
  year4=(year==5)

  n=length(region_code)
  m=rowSums(adj_max)
  D_w=diag(m)
  rho=0.99
  cov_max_inv=D_w-rho*adj_max
  Sigma=matrix(c(1,0.3,0.3,1),nrow=2)
  Sigma_inv=solve(Sigma)
  mcar_sigma=kronecker(cov_max_inv,Sigma_inv)
  estimate=matrix(0,nrow=nsim,ncol=3)
  
  Sigma2=1
  Sigma_inv2=solve(Sigma2)
  mcar_sigma2=kronecker(cov_max_inv,Sigma_inv2)
  
  t1=proc.time()
estimate=foreach(n.sim=1:nsim,.errorhandling=c('remove')) %dopar% {
    ####MCMC###
    set.seed(n.sim)
    W_t=mvrnorm(1,rep(0,region*2),solve(mcar_sigma))
    W_init=scale(matrix(W_t,ncol=2,byrow=T),center=T,scale=F)####random spatial effect initial value###
    
    W_t2=mvrnorm(1,rep(0,region),solve(mcar_sigma2))
    W_init2=scale(matrix(W_t2,ncol=1,byrow=T),center=T,scale=F)####random spatial effect initial value###
    
    spatial_temporal_W=W_init
    x=rnorm(region)
    x=x[region_code]
    x2=runif(n)<0.5
    beta0=0.5
    beta1=1
    beta2=0.5
    beta3=0.5
    sigma=1
    if(idx==1) mu=beta0+beta1*x+W_init2[region_code,1]+beta2*x2+beta3*(year-1)
    if(idx==2) mu=beta0+beta1*x+W_init2[region_code,1]+beta2*x2+beta3*(year-1)^0.5
    if(idx==3) {
      year_effect=rep(NA,5)
      year_effect[1]=0.1
      for(year.j in 2:5){
        year_effect[year.j]=2*year_effect[year.j-1]+rnorm(1,0,1)
      }
      mu=beta0+beta1*x+spatial_temporal_W[region_code,1]+beta2*x2+year_effect[year]
    }
    if(idx==4) mu= mu=beta0+beta1*x+W_init2[region_code,1]+beta2*x2+0.5*year1-0.5*year2+0.6*year3-0.8*year4
    
  lambda_t=exp(-mu/sigma)
  a_t=1/sigma
  b_t=(1/lambda_t)^(1/a_t)
  time=sapply(b_t,function(x) rweibull(1,shape=a_t,scale=x))
  censor=runif(n,0,80)
  y=ifelse(time<=censor,time,censor)
  status=(time<=censor)
  
  eta_init=rWishart(1,2,diag(100,2,2))[,,1]###variance covariance matrix inverse####
  mc_sample=2000
  n_burn=1000
  b0_sample=rep(0,mc_sample)
  b1_sample=rep(0,mc_sample)
  b2_sample=rep(0,mc_sample)
  b3_sample=rep(0,mc_sample)
  sigma_sample=rep(1,mc_sample)
  phi_sample=rep(0,mc_sample*p*region) 
  result=.C("weibul",b0=as.double(b0_sample), b1=as.double(b1_sample), b2=as.double(b2_sample), b3=as.double(b3_sample), sig=as.double(sigma_sample), lamda=as.double(eta_init), w1=as.double(W_init[,1]), w2=as.double(W_init[,2]), D_w=as.double(D_w), C_w=as.double(adj_max), y=as.double(y), s=as.integer(status), N=as.integer(length(y)), x1=as.double(x), x2=as.double(x2), year=as.double(year+1), Cores=as.integer(region_code-1), nend=as.integer(mc_sample), no_regions=as.integer(region),dims=as.integer(c(p,p)),lamda_prec=as.double(diag(0.01,2,2)),deviance=as.double(0),W_sample=as.double(phi_sample),n_burn=as.integer(n_burn))
  W_sample=matrix(result$W_sample,nrow=mc_sample)[-(1:n_burn),]
  W_means=matrix(colMeans(W_sample),nrow=region)
  param_mean=lapply(result,function(x) mean(x[-(1:n_burn)],na.rm=T) )
  D_theta_bar=-2*LL(param_mean$b0,param_mean$b1,param_mean$b2,param_mean$b3,param_mean$sig,W_means,y,status,x,x2,year+1,region_code)
  D_bar=result$deviance/(mc_sample-n_burn)
  DIC=2*D_bar-D_theta_bar
  print(n.sim)
  c(param_mean$b1,param_mean$b2,DIC)
}
  save.image(filname)
  t2=t1-proc.time()
