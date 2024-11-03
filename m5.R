rm(list=ls())
#require(servr)
#servr::make(dir="C:/Users/mliang23/Box/sptime-master/src/")
setwd("/Users/mliang23/Box/sptime-master/")
dyn.load("src/weibull_fivemcar.dll")
#dyn.load("src/liblapack.dll")
library(pracma)
library(xtable)
library(RgoogleMaps)
library(CARBayes)
library(foreign)
library(epitools)
library(data.table)
library(km.ci)
library(survival)
library(DCluster)
library(spdep)
library(ctv)
library(data.table)
library(sp)
library(gstat)
library(lattice)
library(maps)
library(plm)
library(splancs)
library(haven)
library(sas7bdat)
library(ngspatial)
library(sf)
library(terra)
#library(maptools)
require(MASS)
#require(doMC)
#registerDoMC(cores=2)
#idx<-as.numeric(Sys.getenv("LSB_JOBINDEX"))
idx<-3
filname<-paste("m5s",idx,".RData", sep="")
#us_county <- read_sf("/Users/mengluliang/Desktop/sptime-master/updated shape files/4 countydataV4_06_12_2020.shp")
#project2<-"+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +ellps=GRS80    
#   +datum=NAD83 +units=m +no_defs" #USA Contiguous Equidistant Conic Projection
#data.shape<-readShapePoly("/Users/mengluliang/Library/CloudStorage/Dropbox/biostatisitcs_PSU/RA/2020summer/Menglu's files/updated shape files/4 countydataV4_06_12_2020.shp",IDvar="FIPS",proj4string=CRS(project2))
#county_nbr <- poly2nb(us_county, queen=TRUE)
#county_w <- nb2listw(county_nbr,zero.policy=TRUE)
#county_wb <- nb2listw(county_nbr, style="B",zero.policy=TRUE)
#### Create the dissimilarity metric based on county_w;
#us_county$id <- us_county$GEOID ###there are two missing data in fips 01001 01003, so we use GEOID;
#W.mat <- array(0, c(dim(us_county)[1],dim(us_county)[1]))
#for(i in 1:dim(us_county)[1]){
#  temp <- county_w$neighbours[[which(unique(us_county$id)==us_county$id[i])]]
#  ind <- us_county$id %in% unique(us_county$id)[temp]
#  W.mat[i,which(ind==TRUE)] <- rep(1,length(which(ind==TRUE)))
#} 

library(tidyverse)
library(spdep)
#install.packages("devtools")
#devtools::install_github("UrbanInstitute/urbnmapr")
library(urbnmapr)
counties_sf <- get_urbn_map("counties", sf = TRUE)

counties <- counties_sf[counties_sf$state_name %in% c("Florida"),]

counties_polylist <- poly2nb(counties)

Adj_mat <- nb2mat(counties_polylist, style = "B", zero.policy = T)           # Adjacency matrix
m1 <- apply(Adj_mat, 2, sum)
##likelihood###
LL=function(offset,b0,b1,b2,b3,b4,b5,b6,sig,W,y,s,x,x2,year1,year2,year3,year4,region_code,stratums){
    m=offset+b0+b1*x+W[region_code,1]+b2*x2+
      b3*year1[stratums]+b4*year2[stratums]+b5*year3[stratums]+
      b6*year4[stratums]+W[region_code,2]*year1[stratums]+
      W[region_code,3]*year2[stratums]+
    W[region_code,4]*year3[stratums]+W[region_code,5]*year4[stratums]
    #lam=exp(-m/sig)
    #a=1/sig
    #ll=sum(s*log(lam*a*y^(a-1))-lam*y^a)
    ll=sum(na.omit(dpois(y,m,log=F)))
    #lam=exp(m-offset)
    #ll=sum(y*(m-offset)-log(lam))
    #ll=sum(-lam+y*log(lam)-lgamma(y + 1))
    ll
  }


  nsim=100
  p=5
  #cancer_weight=read.table("~/Documents/spatial_temporal/paper_in_preparation/weight.txt",header=T)[,4]
  #adj_max=t(t(c(1,1,1)))%*%c(1,1,1)-diag(1,3,3)
  #adj_max=as.matrix(read.csv("~/Documents/spatial_temporal/paper_in_preparation/adjacent_matrix_PA.csv"))
  #cancer_weight=c(0.5,0.2,0.1,0.01,0.03,0.01,0.15)
  adj_max=as.matrix(Adj_mat)
  region=nrow(adj_max)
  #total_patient=5/min(cancer_weight)
  total_patient=1000
  #total_region_pat=as.integer(total_patient*cancer_weight)
  total_region_pat=as.integer(total_patient*region)
  each_year_pat=as.integer(total_region_pat/5)
  region_code=rep(1:region,max(each_year_pat*5))
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
  library(doParallel)
  estimate<-list()
  #estimate=foreach(n.sim=1:nsim,.errorhandling=c('remove')) %dopar% {
  for(n.sim in 1:nsim) {
  ####MCMC###
  set.seed(12345*n.sim)
  W_t=mvrnorm(1,rep(0,region*2),solve(mcar_sigma))
  
  W_init=scale(matrix(W_t,ncol=2,byrow=T),center=T,scale=F)####random spatial effect initial value###
  
  W_t2=mvrnorm(1,rep(0,region),solve(mcar_sigma2))
  W_init2=scale(matrix(W_t2,ncol=1,byrow=T),center=T,scale=F)####random spatial effect initial value###
  
  spatial_temporal_W=W_init
  x=rnorm(n,0,2)
  offset<-100
  x2=runif(n)<0.5
  beta0=0.5
  beta1=1
  beta2=0.5
  beta3=0.5
  sigma=1
  if(idx==1) mu= offset+beta0+beta1*x+W_init[region_code,1]+beta2*x2+beta3*(year[region_code]-1)
  if(idx==2) mu= offset+beta0+beta1*x+W_init[region_code,1]+beta2*x2+beta3*(year[region_code]-1)^0.5
  
  if(idx==3) {
    year_effect=rep(NA,5)
    year_effect[1]=0.1
    for(year.j in 2:5){
      year_effect[year.j]=2*year_effect[year.j-1]+rnorm(1,0,1)
    }
    mu=offset+beta0+beta1*x+spatial_temporal_W[region_code,1]+beta2*x2+year_effect[year]
  }
  if(idx==4)mu=offset+beta0+beta1*x+W_init[region_code,1]+beta2*x2+0.5*year1-0.5*year2+0.6*year3-0.8*year4
  
  #lambda_t=exp(-mu/sigma)
  #a_t=1/sigma
  #b_t=(1/lambda_t)^(1/a_t)
  b_t=exp(mu)
  #time=sapply(b_t,function(x) rweibull(1,shape=a_t,scale=x))
  #time=sapply(b_t,function(x) rpois(1,x))
  #time=ifelse(is.na(time),0,time)
  censor=runif(n,0,10)
  #y=round(ifelse(time<=censor,time,censor))
  #y<-rbinom(length(mu),offset,exp(mu-offset)/(offset+exp(mu-offset)))
  y=sapply(mu,function(x) rpois(1,x))
  y=ifelse(is.na(y),0,y)
  time<-y
  #status=(time<=censor)
  status=rep(1,length(y))
  #y=sapply(mu,function(x) rpois(1,x))
  #y<-round(na.omit(y))
  eta_init=rWishart(1,2,diag(1,2,2))[,,1]###variance covariance matrix inverse####
  mc_sample=110
  n_burn=10
  b0_sample=rep(0,mc_sample)
  b1_sample=rep(0,mc_sample)
  b2_sample=rep(0,mc_sample)
  b3_sample=rep(0,mc_sample)
  b4_sample=rep(0,mc_sample)
  b5_sample=rep(0,mc_sample)
  b6_sample=rep(0,mc_sample)
  sigma_sample=rep(1,mc_sample)
  phi_sample=rep(1,mc_sample*p*region) 
  estimate[[n.sim]]=.C("weibul",offset =as.double(offset),b0=as.double(b0_sample), b1=as.double(b1_sample), b2=as.double(b2_sample), 
                       b3=as.double(b3_sample), b4=as.double(b4_sample),b5=as.double(b5_sample), b6=as.double(b6_sample),
                       sig=as.double(sigma_sample), lamda=as.double(eta_init), 
                       w1=as.double(rep(0,region)), w2=as.double(rep(0,region)), 
                       w3=as.double(rep(0,region)),w4=as.double(rep(0,region)),
                       w5=as.double(rep(0,region)), 
                       D_w=as.double(D_w), C_w=as.double(adj_max), 
                       y=as.double(y), s=as.integer(status), 
                       N=as.integer(length(y)), x1=as.double(x), 
                       x2=as.double(x2), year1=as.double(year1),year2=as.double(year2),
                       year3=as.double(year3),
                       year4=as.double(year4),  
                       Cores=as.integer(region_code-1),Stratums=as.integer(region_code-1),
                       nend=as.integer(mc_sample), 
                       no_regions=as.integer(region),
                       dims=as.integer(c(p,p)),
                       lamda_prec=as.double(diag(1,2,2)),
                       deviance=as.double(0),
                       W_sample=as.double(phi_sample),n_burn=as.integer(n_burn))
   W_sample=matrix(estimate[[n.sim]]$W_sample,nrow=mc_sample)[-(1:n_burn),]
  W_means=matrix(colMeans(W_sample),nrow=region)
  param_mean=lapply(estimate[[n.sim]],function(x) mean(x[-(1:n_burn)],na.rm=T) )
  param_mean$offset<-ifelse(is.na(param_mean$offset),offset,param_mean$offset)
  D_theta_bar=-2*LL(param_mean$offset,param_mean$b0,param_mean$b1,param_mean$b2,param_mean$b3,param_mean$b4,param_mean$b5,param_mean$b6, param_mean$sig,W_means,y,status,x,x2,year1,year2,year3,year4,region_code,region_code)
  estimate[[n.sim]]$deviance<-ifelse(is.na(estimate[[n.sim]]$deviance)|is.infinite(estimate[[n.sim]]$deviance),0,estimate[[n.sim]]$deviance)
  D_bar=estimate[[n.sim]]$deviance/(mc_sample-n_burn)
  DIC=2*D_bar-D_theta_bar
  print(c(param_mean$b1,param_mean$b2,DIC))
  }
  setwd("C:/Users/mliang23/Box/sptime-master/simulation")
save.image(filname)
  t2=t1-proc.time()
  #load("mliang23/Box/sptime-master/simlation/loglogis_m21.RData")
  