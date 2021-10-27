

if(Sys.info()["sysname"]=="Linux"){
  setwd("/data/kim079/phd")
} else {
  setwd("C:/Users/kim079/Documents/PhD")
}

source("gr4j/R/gr4j_sma.r")
if(!is.loaded("sma_gr4j")){
  if(Sys.info()["sysname"]=="Linux"){
    dyn.load("gr4j/src/gr4j.so")
  } else {
    dyn.load("gr4j/src/gr4j.dll")
  }
}

set.seed(12321)

length_ts<-7
true_input_trial<-runif(length_ts,0,20) #runif(length_ts,100,112)
E_input<-runif(length_ts,0,20) #runif(length_ts,100,105)
input_error_sd<-1e-1 #1e-1
input_error<-rnorm(length(true_input_trial),mean=0,sd=input_error_sd) #sd=0.5
input_trial<-true_input_trial+input_error
input_trial<-pmax(input_trial,0)
actual_state_error<-rnorm(length(input_trial)-1,mean=0,sd=1) #1e-30
discharge_error_sd<-0.01 #0.01
actual_discharge_error<-rnorm(length(input_trial),mean=0,sd=discharge_error_sd) #sd=1e-6
true_input<-data.frame(P=true_input_trial,E=E_input)
actual_model_param<-1000
actual_initial_state<-500

#actual_initial_state<--2000
tt1<-gr4j.sma(actual_model_param,actual_initial_state,actual_state_error,true_input,run_compiled = T,return_state = F)
tt2<-gr4j.sma(actual_model_param,actual_initial_state,actual_state_error,true_input,run_compiled = F,return_state = F)

#matplot(cbind(tt1,tt2),type="l")

tt1-tt2



