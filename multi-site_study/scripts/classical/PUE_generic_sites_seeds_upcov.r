
if(length(commandArgs(TRUE))!=0){
  arguments <- commandArgs(TRUE)
  cat("arguments:",arguments,"\n")
  
  id<-as.character(arguments[[1]])
  output_dir<-as.character(arguments[[2]])
}

if(Sys.info()["sysname"]=="Linux"){
  wd<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2"
} else {
  wd<-"C:/Users/kim079/Documents/model_optimisation_framework"
}


setwd(wd)
# output_dir<-"output/gibbs_sampler_param_uncertainty_on_state_errors_real_data_all_sites_v2"
dir.create(output_dir,showWarnings = F)


preprocess_dir<-"output/gr4j.calib.param.state.all.sites.preprocess"
# id<-415214
theta_files<-list.files(output_dir,pattern=paste0("theta_.*_",id,".csv.gz"),full.names = T)
if(length(theta_files)>0){
  restart_number<-5e5
  prev_seed<-as.numeric(gsub(paste0("theta_|_",id,".csv.gz"),"",basename(theta_files)))
  index_most_recent<-which(prev_seed==max(prev_seed))
  theta_file<-theta_files[index_most_recent]
  SD2_file<-paste0(output_dir,"/SD2_",prev_seed[index_most_recent],"_",id,".csv.gz")
  init_cov_file<-paste0(output_dir,"/CovPar_",prev_seed[index_most_recent],"_",id,".csv.gz")
}



# source("scripts/generate.rating.curve.error.r")
library(mvtnorm)
library(MASS)
library(lattice)
#library(mvnfast)
#library(parallel)
library(condMVNorm)
# library(matrixcalc)
# library(data.table)


source("packages/gr4j_with_routing_error/R/gr4j_sma_routing.r")
if(!is.loaded("routing_gr4j_sk")){
  if(Sys.info()["sysname"]=="Linux"){
    dyn.load("packages/gr4j_with_routing_error/src/gr4j_with_routing_error.so")
  } else {
    dyn.load("packages/gr4j_with_routing_error/src/gr4j_with_routing_error.dll")
  }
}

sd_zero_mean<-function(x) sqrt(mean(x^2))


logprior_fun<-function(params,min,max,data){
  
  # log_dens1<-dunif(params[1],min[1],max[1],log=T)
  # log_dens2<-dunif(params[2],min[2],max[2],log=T)
  # log_dens3<-dunif(params[3],min[3],max[3],log=T)
  # log_dens4<-dunif(params[4],min[4],max[4],log=T)
  # log_dens5<-dunif(params[5],min[5],max[5],log=T)
  # log_dens6<-dunif(params[6],min[6],max[6],log=T)
  # log_dens7<-dunif(params[7],min[7],max[7],log=T)
  # log_dens8<-dunif(params[8],min[8],max[8],log=T)

  # return(log_dens1+log_dens2+log_dens3+log_dens4+log_dens5+log_dens6+log_dens7+log_dens8)
  
  log_dens<-0
  for(ld in 1:length(params)){
    log_dens<-log_dens+dunif(params[ld],min[ld],max[ld],log=T)
  }
  
  return(log_dens)
}

loglikelihood_fun<-function(data,params){
  Q1<-gr4j.run(param=params[1:4], initial_state_S=params[5]*params[1], initial_state_R=params[6]*params[3], 
               state_error=rep(0,nrow(data$cal.input)-1), input=data$cal.input,
               state_error_R=rep(0,nrow(data$cal.input)-1),
               run_compiled = T,return_state=F)/1000*data$area_m2/86400 #convert to cumecs
  if(is.na(Q1[1])){
    return(-Inf)
  } else {
    # resid<-Q1-data$obs
    # plot(resid,type="l")
    
    # AR1<-resid*params[7]
    # AR1_removed<-resid[2:length(resid)]-AR1[1:(length(AR1)-1)]
    
    # AR1_removed<-c(0,AR1_removed)
    # test<-resid[1]
    # for(i in 2:length(resid)){
    #   test<-c(test,test[i-1]*params[7]+AR1_removed[i])
    # }
    # sum(abs(test-resid))
    
    # acf(resid)
    # acf(AR1_removed)
    # plot(AR1_removed,type="l")
    # hist(AR1_removed)
    # log_dens<-sum(dnorm(AR1_removed,0,params[8],log=T)) # sd=1
    # log_dens<-sum(dnorm(resid,0,1,log=T))
    thin_number<-30 #6 # TODO: adjust this by examining acf
    thin<-seq(1,length(Q1),by=thin_number)
    predict_flow_thinned<-Q1[thin]
    obs_thinned<-data$obs[thin]
    # acf(predict_flow_thinned)
    # acf(obs_thinned)
    
    # sigma_error<-params[7]+params[8]*predict_flow_thinned
    sigma_error<-params[7]*(predict_flow_thinned^params[8])
    # sigma_error<-exp(log(params[7])+log(predict_flow_thinned)*params[8])
    # acf(sigma_error)
    if(is.na(sigma_error[1]) | all(sigma_error<=0)){
      return(-Inf)
    } else {
      sigma_error[sigma_error<=0]<-min(sigma_error[sigma_error>0])
      
      threshold<--100 #0.09 # TODO: adjust this based on rating curve
      
      is_pos_sigma_error<-sigma_error>0
      indices_over_threshold<-which(obs_thinned>threshold & is_pos_sigma_error)
      
      resid_thinned<-predict_flow_thinned-obs_thinned
      # plot(resid_thinned,type="b")
      # plot(density(resid_thinned))
      # plot(y=resid_thinned,x=predict_flow_thinned,log="x")
      # quants<-quantile(predict_flow_thinned,seq(0,1,by=0.01))
      # quant_sds<-c()
      # for(q in 2:length(quants)){
      #   index_quant<-which(predict_flow_thinned<=quants[q])
      #   quant_sds<-c(quant_sds,sd_zero_mean(resid_thinned[index_quant]))
      # }
      # plot(x=quants[-1],y=quant_sds,log="x")
      # points(y=sigma_error,x=predict_flow_thinned,col=2)
      # 
      # plot(x=quants[-1],y=quant_sds,log="xy")
      # points(y=sigma_error,x=predict_flow_thinned,col=2)
      # 
      # plot(x=quants[-1],y=quant_sds)
      # points(y=sigma_error,x=predict_flow_thinned,col=2)
      # 
      # lm_loglog<-lm(log(quant_sds)~log(quants[-1]))
      # modelled_sds<-exp(lm_loglog$coefficients[1])*(quants[-1]^lm_loglog$coefficients[2])
      # lines(y=modelled_sds,x=quants[-1],col=3)
      # 
      # # recalculate sd
      # sigma_error_trial<-exp(lm_loglog$coefficients[1])*(predict_flow_thinned^lm_loglog$coefficients[2])
      # sum(dnorm(resid_thinned[indices_over_threshold],0,sigma_error_trial[indices_over_threshold],log=T))
      # 
      # plot(obs_thinned,type="l")
      # lines(predict_flow_thinned,col=2,lty=2)
      # browser()
      # layout(1)
      # acf(resid_thinned)
      # length(resid_thinned)
      # plot(y=abs(resid_thinned),x=predict_flow_thinned,log="xy")
      # plot(y=abs(resid_thinned),x=predict_flow_thinned)
      # plot(resid_thinned/sigma_error,type="l")
      # plot(density(resid_thinned/sigma_error))
      # acf(resid_thinned/sigma_error)
      # acf(resid_thinned)
      # 
      # acf(dnorm(resid_thinned[indices_over_threshold],0,sigma_error[indices_over_threshold],log=T))
      
      log_dens<-sum(dnorm(resid_thinned[indices_over_threshold],0,sigma_error[indices_over_threshold],log=T))
      # plot(sigma_error[indices_over_threshold],type="b",log="y")
      # points(abs(resid_thinned[indices_over_threshold]),col=2)
      # plot(dnorm(resid_thinned[indices_over_threshold],0,sigma_error[indices_over_threshold],log=T),type="l")
      if(is.na(log_dens)){
        browser()
      } else {
        # if(log_dens>6) browser()
        return(log_dens)
      }

    }

  }

  
}
plot_function<-function(){
  
  if(export_diagnostic){
    # write parameters
    output_name<-paste0(output_dir,"/diagnostic_",prefix,".png")
    png(output_name,height=1000)
  }
  cat("loglikelihood:",L[i],"\n")
  # layout(matrix(1:10,nrow=5,byrow=T))
  layout(matrix(1:10,ncol=2,byrow=T))
  
  logposterior<-L+Pr
  if(!all(is.infinite(logposterior[!is.na(logposterior)]))){
    plot(logposterior[1:i],type="l")
    
  }    
  
  acceptance_rate<-length(which(Jump[!is.na(Jump)]==1))/length(Jump[!is.na(Jump)])
  cat("total acceptance_rate =",acceptance_rate,"\n")
  
  if(i>10000){
    for(pp in 1:length(param.ranges)){
      plot(density(theta[10000:i,pp]))
    }
  } else {
    for(pp in 1:length(param.ranges)){
      plot(density(theta[1:i,pp]))
    }

  }
  
  
  
  if(export_diagnostic){
    dev.off()
  }
}


# from PyMC3 (metropolis.py)
tune<-function(scale, acc_rate){
  if(acc_rate < 0.001){
    # reduce by 90%
    scale<-scale*0.1
  } else if(acc_rate < 0.05){
    # reduce by 50%
    scale<-scale*0.5
  } else if(acc_rate < 0.2){
    # reduce by 10%
    scale<-scale*0.9
  } else if(acc_rate > 0.95){
    # increase by factor of 10
    scale<-scale*10
  } else if(acc_rate > 0.75){
    # increase by double
    scale<-scale*2
  } else if(acc_rate > 0.5){
    # increase by 10%
    scale<-scale*1.1
  }
  return(scale)
}


# Start Data setup ###################################
set.seed(12321)


export_diagnostic<-T

input_ts_file<-paste0(preprocess_dir,"/state_error_simulation_data_",id,".csv")
param_file<-paste0(preprocess_dir,"/gr4j_params_",id,".csv")

input_ts<-read.csv(input_ts_file,as.is=T)
orig_params<-read.csv(param_file,as.is=T)


# input_ts$Q_with_state_errors

sd_zero_mean<-function(x) sqrt(mean(x^2))

# End Data setup #####################################


data<-list(cal.input=data.frame(P=input_ts$P,E=input_ts$E),
           area_m2=orig_params$area_m2,obs=input_ts$obs) #obs=input_ts$Q_no_state_errors #obs=input_ts$Q_with_state_errors

# param.ranges<-list(c(1,10000),c(-30,30), c(1, 1000), c(0.5, 5), c(0,10000), c(0, 1000), c(0,1))
# param.ranges<-list(c(1,10000),c(-30,30), c(1, 1000), c(0.5, 5), c(0,1), c(0, 1), c(0,1))
# param.ranges<-list(c(1,10000),c(-30,30), c(1, 1000), c(0.5, 5), c(0,1), c(0, 1), c(0,1), c(1e-3,1000))
param.ranges<-list(c(1,10000),c(-30,30), c(1, 1000), c(0.5, 5), c(0,1), c(0,1), c(-1000,1000), c(-1000,1000))
params_min<-sapply(param.ranges,function(x) x[1])
params_max<-sapply(param.ranges,function(x) x[2])
# initial_params<-sapply(param.ranges,mean)
# initial_params<-c(as.numeric(orig_params[1,-ncol(orig_params)]),0.79)
numeric_orig_params<-as.numeric(orig_params[1,-ncol(orig_params)])
# initial_params<-c(numeric_orig_params[1:4],numeric_orig_params[5]/numeric_orig_params[1],numeric_orig_params[6]/numeric_orig_params[3],0.79)
# initial_params<-c(numeric_orig_params[1:4],numeric_orig_params[5]/numeric_orig_params[1],numeric_orig_params[6]/numeric_orig_params[3],0.79,10)
# initial_params<-c(numeric_orig_params[1:4],numeric_orig_params[5]/numeric_orig_params[1],numeric_orig_params[6]/numeric_orig_params[3],0,1)
initial_params<-c(numeric_orig_params[1:4],numeric_orig_params[5]/numeric_orig_params[1],numeric_orig_params[6]/numeric_orig_params[3],exp(-2.078),0.829)
# initial_params<-c(numeric_orig_params[1:4],numeric_orig_params[5]/numeric_orig_params[1],numeric_orig_params[6]/numeric_orig_params[3],10,10)


# initial_params<-c(588.8467, 17.43998, 430.1803, 4.553398, 1.255141, 0.1678335,  94.21935,  7.417855)


logprior_init<--Inf
loglike_init<--Inf

seed_number<-unlist(strsplit(output_dir,"_"))
seed_number<-as.numeric(seed_number[length(seed_number)])
set.seed(seed_number)

initial_params_counter<-0
while(is.infinite(logprior_init) | is.infinite(loglike_init)){
  initial_params_counter<-initial_params_counter+1
  if(initial_params_counter>1e6) stop("initialisation shouldn't take this long")
  for(i in 1:6){
    initial_params[i]<-runif(1,param.ranges[[i]][1],param.ranges[[i]][2])
  }
  
  logprior_init<-logprior_fun(initial_params,params_min,params_max,data)
  loglike_init<-loglikelihood_fun(data,initial_params)[[1]]
}



if(is.infinite(logprior_init+loglike_init)) stop("Initialisation failed!")


num_timesteps_cov<-Inf # 1000000
cor_plot_interval<-20000
tune_scale<-T # tunes the scaling factor (SD2) according to the acceptance rate over the last tune_interval. From pymc3 (metropolis.py)
tune_interval<-2000 # this should be divisible by update_covariance_interval
update_covariance_interval<-20 #100 # this should be divide into tune_interval
start_collecting_CovPar<-Inf #1000001
start_sampling_CovPar<-Inf #1000001
resample_thetas<-F
start_gibbs_iter<-1



# restart params ################
additional_iter<-1000000
seed<-as.numeric(Sys.time())

if(exists("SD2_file")) SD2<-as.numeric(readLines(SD2_file))

if(!exists("all_CovPar_dir")){
  all_CovPar_dir<-paste0(output_dir,"/all_CovPar")
}
dir.create(all_CovPar_dir,showWarnings = F)
all_CovPar_files<-list.files(all_CovPar_dir,pattern="CovPar_",full.names=T)
all_CovPar<-list()
all_CovPar_means<-NULL
if(length(all_CovPar_files)>0 & is.finite(start_collecting_CovPar)){
  for(i in 1:length(all_CovPar_files)){
    all_CovPar[[i]]<-as.matrix(read.csv(all_CovPar_files[i],as.is=T,header=F))
  }
  theta_means_file<-list.files(all_CovPar_dir,pattern="theta_means_",full.names = T)
  all_CovPar_means<-read.csv(theta_means_file,as.is=T)
  start_collecting_CovPar<-1
  start_sampling_CovPar<-1
}



if(!exists("ITER")) ITER = 1000000  #number of iterations to perform
i0   = 0.001 #0.10  #percentage of initial iterations before adaptation
# rm("SD2")
if(exists("SD2")){
  SD1  = SD2  #initial covariance matrix scaling factor (i0)
} else {
  SD1  = 0.50  #initial covariance matrix scaling factor (i0)
  SD2  = (2.4^2)/length(initial_params) #from Haario #0.30 #0.009 #0.15  #adaptive covariance matrix scaling factor (1-i0) (lower scaling is higher acceptance) )
}
if(!exists("stop_update_covariance")) stop_update_covariance<-0.5

if(export_diagnostic){
  diagnostic_plot_interval<-ITER
} else {
  diagnostic_plot_interval<-1000
}

#ncores<-detectCores(logical=F)

### 2 - Simulation model parameters
#PAR1 = 0.5	#place holder
#PAR2 = 15.00	#place holder
### 3 - Likelihood function parameters
#VARP = 0.185	#variance parameter
### 4 - Define parameter matrix
INIT_PAR<-initial_params
if(exists("theta_file")){
  theta<-as.matrix(read.csv(theta_file,as.is=T))
  prev_theta<-theta[!is.na(theta[,1]),]
  if(is.null(nrow(prev_theta))){
    INIT_PAR<-prev_theta
  } else if(nrow(prev_theta)==1) {
    INIT_PAR<-prev_theta[nrow(prev_theta),]
  } else {
    INIT_PAR<-prev_theta[nrow(prev_theta)-1,]
  }
}


# theta<-as.data.table(matrix(NA,nrow=ITER,ncol=length(INIT_PAR)))
theta<-matrix(NA,nrow=ITER,ncol=length(INIT_PAR))

if(exists("restart_number")){
  if(is.na(restart_number)){
    ITER<-nrow(prev_theta)+additional_iter
    theta<-rbind(prev_theta,theta)
    start_iter<-nrow(prev_theta)+1
  } else {
    indices_to_replace<-(nrow(prev_theta)-restart_number):(nrow(prev_theta)-1)
    theta[1:length(indices_to_replace),]<-prev_theta[indices_to_replace,]
    start_iter<-length(indices_to_replace)+1
  }
  
} else {
  start_iter<-2
  theta[1,] = INIT_PAR #first row using initial values
}

nPAR = ncol(theta) #determines number of parameters
### 5 - Define parameter names for output
#ParName = c('Slope','Intercept','Variance') #redefine as appropriate

## Define parameter priors -----------------------------------------------#
# Define (log) prior to be used for each calibrated parameter.
# Pr1 = dbeta(PAR1,0.5,1.0,log=T)
# Pr2 = dgamma(PAR2,0.5,scale=1.0,log=T)
# Pr3 = dexp(VARP,0.5,log=T)
# Pr<-rep(NA,ITER)
# Pr[1] = Pr1+Pr2+Pr3
Pr<-rep(NA,ITER)
Pr[start_iter-1] = logprior_fun(theta[start_iter-1,],params_min,params_max,data)
if(is.infinite(Pr[start_iter-1])) stop("Initial prior is infinite!!")

## Initialize covariance matrix ------------------------------------------#
### 1 - Define best guess variances for each calibrated parameter
VarPar = rep(1E-6,nPAR)
# A <- matrix(VarPar,ncol=nPAR) # allocates some memory for rmvn function

### 2 - Populate initial covariance matrix
if(!exists("init_cov_file")){
  CovPar<-matrix(0,nrow=length(VarPar),ncol=length(VarPar))
  for(i in 1:nPAR){
    CovPar[i,i] = SD1*VarPar[i]
  }
} else {
  CovPar<-as.matrix(read.csv(init_cov_file,as.is=T,header=F))
  if(!all(dim(CovPar)==c(length(VarPar),length(VarPar)))){
    stop("dimensions of initial covariance matrix is wrong")
    CovPar<-matrix(0,nrow=length(VarPar),ncol=length(VarPar))
  }
}

### 3 - Define other covariance terms needed later
epsilon = 1e-8 #1e-20 #small number preventing CovPar from becoming singular during adaptation
Id = diag(nPAR)  #identity matrix used in CovPar adaptation

recalculate_CovPar<-T
if(recalculate_CovPar & exists("prev_theta")){
  CovPar_unscaled<-cov(theta[max(1,start_iter-num_timesteps_cov):(start_iter-1),])
  CovPar = SD2*CovPar_unscaled+SD2*epsilon*Id
}


## Initial run of simulation model ---------------------------------------#
# Replace the following line with your function as appropriate, where
# y is the output, SimModel is the function name, PAR1 & PAR2 are
# simulation model parameters
# SimModel<-function(PAR1,PAR2){
#   1:10*PAR1+PAR2
# }
# y = SimModel(theta[1,],data)

# Replace the following line with the appropriate log likelihood
# function, where obs is the observed data, pred is the predicted data,
# and varp is the variance parameter.

# normal_loglikelihood<-function(obs,pred,varp){
#   sum(dnorm(obs-pred,mean=0,sd=sqrt(varp),log=T))
# }
# obs<-15:24
# varp<-VARP

logL = loglikelihood_fun(data,theta[start_iter-1,])
L<-rep(NA,ITER)

L[start_iter-1] = logL[[1]] #logL$combined_log_likehood #initial value for the log likelihood

## Run AM algorithm ------------------------------------------------------#
# theta_pro<-matrix(NA,nrow=ITER,ncol=nPAR)
L_pro<-rep(NA,ITER)
Pr_pro<-rep(NA,ITER)
psi<-rep(NA,ITER)
Jump<-rep(NA,ITER)
Jumps = 0 #jump counter
beginning_time<-Sys.time()
all_cor_theta<-matrix(NA,nrow=ITER/cor_plot_interval,ncol=nPAR*nPAR)
all_SD2<-rep(NA,ITER/tune_interval)
all_interval_acceptance<-rep(NA,ITER/tune_interval)
all_discharge_error_sd<-c()
first_covpar_calculated<-F
prefix<-paste(seed,id,sep="_")
iii<-0
param_indices<-1:nPAR

for(i in start_iter:ITER){
  if(i%%1000==0) cat(i,"/",ITER,"|","\n")
  # if(i%%1000==0) browser()
  iii<-iii+1
  if(iii>nPAR) iii<-1
  
  ### 1 - Adapt covariance matrix
  if (i > i0*ITER & i%%update_covariance_interval==0 & i<stop_update_covariance*ITER){  #adaptive covariance matrix routine
    thin_for_covar<-F
    if(thin_for_covar){
      theta_tmp<-theta[max(1,i-num_timesteps_cov):(i-1),]
      theta_tmp<-theta_tmp[seq(1,nrow(theta_tmp),by=100),]
      CovPar_unscaled<-cov(theta_tmp)
    } else {
      use_updateCovariance<-T
      if(use_updateCovariance & exists("CovPar_unscaled") & first_covpar_calculated & is.infinite(num_timesteps_cov)){
        library(onlinePCA)
        CovPar_unscaled_new<-updateCovariance(CovPar_unscaled,
                                              theta[max(1,i-update_covariance_interval):(i-1),],
                                              i-update_covariance_interval-1,
                                              colMeans(theta[1:max(1,i-update_covariance_interval-1),]))
        
        if(i%%(update_covariance_interval*10000)==0){
          actual_CovPar_unscaled<-cov(theta[1:(i-1),])
          if(!all.equal(CovPar_unscaled_new,actual_CovPar_unscaled)){
            cat("covariance update problem \n")
            browser()
            stop("covariance update problem")
            CovPar_unscaled_new<-actual_CovPar_unscaled
          }
        }
        CovPar_unscaled<-CovPar_unscaled_new
        
      } else {
        CovPar_unscaled<-cov(theta[max(1,i-num_timesteps_cov):(i-1),])
        first_covpar_calculated<-T
      }
      
      
    }
    
    CovPar = SD2*CovPar_unscaled+SD2*epsilon*Id #updates the covariance
    if(i>=start_collecting_CovPar){
      all_CovPar[[length(all_CovPar)+1]]<-CovPar
      all_CovPar_means<-rbind(all_CovPar_means,theta[i-1,])
    }
    if(i>=start_sampling_CovPar & i%%(update_covariance_interval*100)==0){
      sample_CovPar<-sample(length(all_CovPar),1)
      CovPar<-all_CovPar[[sample_CovPar]]
      if(resample_thetas){
        theta[i-1,]<-all_CovPar_means[sample_CovPar,]
      }
      
    }
  }
  # Delayed Rejection ###########
  number_delayed_rejections<-4
  attempt<-0
  DR_scaling_factor<-1
  while(attempt<=number_delayed_rejections){
    attempt<-attempt+1
    ### 2 - Generate parameter proposals
    if(i>=start_gibbs_iter){
      
      theta[i,]<-theta[i-1,]
      
      # check positive definiteness
      # if(!is.positive.definite(CovPar)){
      #   stop("CovPar not positive definite")
      # }
      theta[i,iii]<-rcmvnorm(n=1,
                             mean=theta[i-1,],
                             sigma=CovPar*DR_scaling_factor,
                             dependent.ind = iii,
                             given.ind = param_indices[-iii],
                             X.given = theta[i-1,-iii],
                             method="chol",check.sigma = F)
    } else {
      theta[i,]<-mvrnorm(1,theta[i-1,],CovPar*DR_scaling_factor) #proposed values - more robust faster
    }
    #theta[i,] = rmvnorm(1,theta[i-1,],CovPar) #proposed values - less robust slower
    
    DR_scaling_factor<-DR_scaling_factor*0.5 # for next delayed rejection stage
    
    # # adjust covpar
    # input_error_indices<-(length_ts+1):(length_ts*2)
    # for(ii in input_error_indices){
    #   # CovPar[ii,ii] = SD2*(input_error_sd^2)+SD2*epsilon
    #   CovPar[ii,ii] = (input_error_sd^2)
    # }
    # prop_mean<-theta[i-1,]
    # prop_mean[input_error_indices]<-0
    # theta[i,]<-mvrnorm(1,prop_mean,CovPar) #proposed values - more robust faster
    # # plot(mvrnorm(1,theta[i-1,],CovPar))
    
    
    #system.time({for(jj in 1:1000000) theta[i,]<-mvrnorm(1,theta[i-1,],CovPar)})
    #system.time({for(jj in 1:1000000) theta[i,]<-rmvn(1,theta[i-1,],CovPar,ncores=ncores)})
    #system.time({for(jj in 1:1000000) rmvn(1,theta[i-1,],CovPar,ncores=ncores,A=A); theta[i,]<-A})
    #rmvn(1,theta[i-1,],CovPar,ncores=ncores,isChol=T,A=A); theta[i,]<-A #proposed values - faster
    #rmvn(1,theta[i-1,],CovPar,ncores=ncores,isChol=F,A=A); theta[i,]<-A #proposed values - faster
    
    # # replace the input error proposal with sample from actual distribution
    # input_error_indices<-(length_ts+1):(length_ts*2)
    # input_error_sample<-rnorm(length(input_error_indices),mean=0,sd=input_error_sd)
    # # sd_zero_mean(input_error_sample)
    # # TODO: scale this
    # theta[i,input_error_indices]<-input_error_sample
    
    # Ensure within range
    #theta[i,theta[i,]<params_min]<-params_min[theta[i,]<params_min]
    #theta[i,theta[i,]>params_max]<-params_max[theta[i,]>params_max]
    #theta_restrict<-theta[i,]
    #theta_restrict[theta[i,]<params_min]<-params_min[theta[i,]<params_min]
    #theta_restrict[theta[i,]>params_max]<-params_max[theta[i,]>params_max]
    
    # theta_pro[i,] <- theta[i,] #matrix of proposed theta's (rejected & accepted)
    
    ### 3 - Run simulation model
    #y = SimModel(theta[i,],data)
    
    ### 5 - Compute log prior
    #   Pr1 = dbeta(theta[i,1],0.5,1.0,log=T)
    #   Pr2 = dgamma(theta[i,2],0.5,scale=1.0,log=T)
    #   Pr3 = dexp(theta[i,3],0.5,log=T)
    #   Pr[i] = Pr1+Pr2+Pr3
    Pr[i] = logprior_fun(theta[i,],params_min,params_max,data)
    #Pr[i] = logprior_fun(theta_restrict,params_min,params_max,data)
    Pr_pro[i] = Pr[i] #proposed prior
    
    ### 4 - Compute log likelihood
    if(is.infinite(Pr[i])){
      logL = -Inf
    } else {
      logL = loglikelihood_fun(data,theta[i,])
    }
    # logL = loglikelihood_fun(data,theta[i,])
    #logL = loglikelihood_fun(data,theta_restrict)
    
    L[i] = logL[[1]] #logL$combined_log_likehood
    L_pro[i] = L[i] #proposed likelihood
    
    ### 6 - Compute Metropolis ratio
    psi[i] = exp((L[i]-L[i-1])+(Pr[i]-Pr[i-1])) #Metropolis ratio
    if(is.nan(psi[i])) psi[i]<-0
    ### 7 - Determine accept/reject of proposal
    z = runif(1,0,1)
    if(z <= psi[i]){
      #theta[i,] = theta[i,] #jump to next theta
      Jumps = Jumps+1
      Jump[i] = 1
      break
    } else{
      theta[i,] = theta[i-1,] #remain at current theta
      L[i] = L[i-1]
      Pr[i] = Pr[i-1]
      Jump[i] = 0
    }
  }
  
  ## 8 - Iterate
  if(i%%diagnostic_plot_interval==0){
    plot_function()
  }
  
  if(tune_scale & i%%tune_interval==0){
    Jump_interval<-Jump[(i-tune_interval+1):i]
    interval_acc_rate<-length(which(Jump_interval==1))/length(Jump_interval)
    
    SD2<-tune(SD2,interval_acc_rate)
    cat("interval acceptance_rate =",interval_acc_rate,"SD2 =",SD2,"\n")
    all_SD2[ITER/tune_interval]<-SD2
    all_interval_acceptance[ITER/tune_interval]<-interval_acc_rate
    CovPar = SD2*CovPar_unscaled+SD2*epsilon*Id
  }
  
  
  
}
elapsed_time<-difftime(Sys.time(),beginning_time,units="hours")
cat("Elapsed time: ",elapsed_time," ",units(elapsed_time),"\n")

# write the output

if(Sys.info()[1]=="Windows"){
  gzip<-"software/gzip/gzip.exe -f"
} else {
  gzip<-"gzip -f"
}

# write parameters
output_name<-paste0(output_dir,"/theta_",prefix,".csv")
write.csv(theta,output_name,row.names=F,quote=F)
system(paste(gzip,output_name))

# write SD2
output_name<-paste0(output_dir,"/SD2_",prefix,".csv")
write.table(SD2,output_name,row.names=F,quote=F,col.names=F,sep=",")
system(paste(gzip,output_name))


# output_name<-paste0(output_dir,"/scaling_R_",prefix,".png")
# png(output_name)
# layout(1)
# plot(theta[1:i,length_ts+number_of_sampled_inputs+1],type="l")
# dev.off()

# write Covariance
output_name<-paste0(output_dir,"/CovPar_",prefix,".csv")
write.table(CovPar,output_name,row.names=F,quote=F,col.names=F,sep=",")
system(paste(gzip,output_name))

# write all covariances
# empty the folder
unlink(list.files(all_CovPar_dir,full.names = T))
if(length(all_CovPar)>0){
  for(i in 1:length(all_CovPar)){
    output_name<-paste0(all_CovPar_dir,"/CovPar_",i,"_",prefix,".csv")
    write.table(all_CovPar[[i]],output_name,row.names=F,quote=F,col.names=F,sep=",")
    system(paste(gzip,output_name))
  }
  output_name<-paste0(all_CovPar_dir,"/theta_means_",prefix,".csv")
  write.csv(all_CovPar_means,output_name,row.names=F,quote=F)
  system(paste(gzip,output_name))
}


# if(Sys.info()[1]=="Linux"){
#   subject<-prefix
#   message<-prefix
#   email_command<-sprintf("echo \"%s\" | mail -s \"%s\" shaunsanghokim@gmail.com",message,subject)
#   system(email_command)
# }

# stop()
# # write log likelihood
# output_name<-paste0(output_dir,"/loglikelihood_",prefix,".csv")
# write.table(L,output_name,row.names=F,quote=F,col.names=F,sep=",")
# system(paste(gzip,output_name))
# 
# # write log prior
# output_name<-paste0(output_dir,"/logprior_",prefix,".csv")
# write.table(Pr,output_name,row.names=F,quote=F,col.names=F,sep=",")
# system(paste(gzip,output_name))
# 
# #write log posterior
# output_name<-paste0(output_dir,"/logposterior_",prefix,".csv")
# write.table(L+Pr,output_name,row.names=F,quote=F,col.names=F,sep=",")
# system(paste(gzip,output_name))
# 
# # write jump
# output_name<-paste0(output_dir,"/Jump_",prefix,".csv")
# write.table(Jump,output_name,row.names=F,quote=F,col.names=F,sep=",")
# system(paste(gzip,output_name))


# # write correlation timeseries
# output_name<-paste0(output_dir,"/cor_timeseries_",prefix,".csv")
# write.table(all_cor_theta,output_name,row.names=F,quote=F,col.names=F,sep=",")
# system(paste(gzip,output_name))


# fraction of accepted jumps
acceptance_rate<-length(which(Jump[!is.na(Jump)]==1))/length(Jump[!is.na(Jump)])
cat("acceptance_rate =",acceptance_rate,"\n")
cat("Finished!! \n")

# 
# # stop()
# # write.table(CovPar,"scripts/CovPar.csv",row.names=F,quote=F,col.names=F,sep=",")
# # 
# # write.csv(theta,"output/state_uncertainty/AM/theta.csv",row.names=F,quote=F)
# # 
# # # fraction of accepted jumps
# # acceptance_rate<-length(which(Jump[!is.na(Jump)]==1))/length(Jump[!is.na(Jump)])
# # cat("acceptance_rate =",acceptance_rate,"\n")
# # 
# # 
# # correlations
# library(lattice)
# cor_theta<-cor(theta)
# layout(1)
# lp<-levelplot(cor_theta,ylim=c(nrow(cor_theta)+0.5,0.5),at=seq(-1,1,length.out=51))
# print(lp)
# # 
# # # library(car)
# # # scatterplotMatrix(theta[1:100,1:3])
# # # pairs(theta[,1:3])
# # 
# # library(PerformanceAnalytics)
# # chart.Correlation(theta,method="pearson",histogram=T,pch=16)
# # 
# # # get rid of nas
# # theta_nona<-theta[!is.na(theta[,1]),]
# # 
# # indices of different data types
# groups_indices<-list(1,2:(length(actual_state_error)+1),
#                      (length(actual_state_error)+2):(length(actual_state_error)+length(input_trial)+1))
# # 
# # # plot intitial state, state error 1, input error 1
# # layout(1:3)
# # plot(theta_nona[,groups_indices[[1]]],type="l")
# # plot(theta_nona[,groups_indices[[2]][1]],type="l",main="2nd time step state error")
# # plot(theta_nona[,groups_indices[[3]][1]],type="l")
# # 
# # 
# # # Likelihood removed nas
# # layout(1:3)
# L_no_na<-L[!is.na(L)]
# # plot(L_no_na,type="l",main="log likelihood")
# # 
# # plot(Pr,type="l",main="log prior")
# # 
# # Calculate posterior
# Pr_no_na<-Pr[!is.na(Pr)]
# logposterior<-L_no_na+Pr_no_na
# # plot(logposterior,type="l",main="log posterior")
# # write.csv(logposterior,"output/state_uncertainty/AM/posterior.csv",row.names=F,quote=F)
# # 
# # unnormalise theta
# #theta_unnorm<-t(apply(theta[!is.na(theta[,1]),],1,inv.normalise,factor=normalisers))
# state_normaliser<-10^theta[1:i,length_ts+length_ts+1]
# all_state_normaliser<-matrix(rep(state_normaliser,length_ts-1),ncol=length_ts-1)
# all_normaliser<-cbind(rep(data$normalisers[1],length(state_normaliser)),
#                       all_state_normaliser,
#                       matrix(1,nrow=length(state_normaliser),ncol=length_ts))
# 
# theta_unnorm<-cbind(inv.normalise(theta[1:i,-ncol(theta)],all_normaliser),state_normaliser)
# 
# # layout(1:3)
# # calculate state error sd
# state_sd_calc<-apply(theta_unnorm[,groups_indices[[2]]],1,sd_zero_mean)
# plot(state_sd_calc,type="l",main="state error sd",log="y")
# 
# burn.in<-100000
# thin<-100
# theta_p<-theta_unnorm[seq(from=min(nrow(theta_unnorm),burn.in),to=nrow(theta_unnorm),by=thin),]
# #acf(theta_p[,15])
# state_sd_calc_p<-apply(theta_p[,groups_indices[[2]]],1,sd_zero_mean)
# layout(1:2)
# boxplot(state_sd_calc_p,main="state error sd burn-in thinned")
# boxplot(state_sd_calc_p,main="state error sd burn-in thinned (log y)",log="y")
# 
# 
# 
# # 
# # # calculate input error sd
# # input_sd_calc<-apply(theta_unnorm[,groups_indices[[3]]],1,sd_zero_mean)
# # plot(input_sd_calc,type="l",main="input error sd")
# # 
# # # calculate discharge error sd
# # discharge_sd_calc<-c()
# # for(ii in 1:nrow(theta_unnorm)){
# #   if(ii%%1000==0) cat(ii,"/",nrow(theta_unnorm),"\n")
# #   input_run<-data.frame(P=input_trial-as.numeric(theta_unnorm[ii,groups_indices[[3]]]),E=E_input)
# #   model_run<-gr4j.sma(data$model_param,as.numeric(theta_unnorm[ii,1]),as.numeric(theta_unnorm[ii,groups_indices[[2]]]),input_run)
# #   error_discharge<-model_run-data$obs_discharge
# #   discharge_sd_calc<-c(discharge_sd_calc,sd_zero_mean(error_discharge))
# # }
# # plot(discharge_sd_calc,type="l",main="discharge error sd")
# # 
# # layout(1)
# # boxplot(state_sd_calc,main="state sd calc",log="y")
# # abline(h=sd_zero_mean(actual_state_error),col=2,lty=2)
# # 
# # cat("actual_state_error_sd=",sd_zero_mean(actual_state_error),"\n")
# # cat("input_error_sd=",sqrt(mean(input_error^2)),"\n")
# # cat("actual_discharge_error_sd=",sqrt(mean(actual_discharge_error^2)),"\n")
# 
# # colormap
# library(RColorBrewer)
# #library(colorRamps)
# #col.ramp<-colorRampPalette(c("white","blue","red"))(100)
# col.ramp<-rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))
# layout(1)
# smoothScatter(x=state_sd_calc,y=logposterior[1:i],nrpoints=0,colramp=colorRampPalette(col.ramp))
# 
# 
# #theta_unnorm<-t(apply(prev_theta[!is.na(prev_theta[,1]),],1,inv.normalise,factor=normalisers))
# #ppp<-prev_theta
