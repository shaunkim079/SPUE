

if(Sys.info()["sysname"]=="Linux"){
  setwd("/data/kim079/model_optimisation_framework_v2/packages")
} else {
  setwd("C:/Users/kim079/Documents/model_optimisation_framework/packages")
}

source("gr4j_with_routing_error/R/gr4j_sma_routing.r")
if(!is.loaded("routing_gr4j")){
  if(Sys.info()["sysname"]=="Linux"){
    dyn.load("gr4j_with_routing_error/src/gr4j_with_routing_error.so")
  } else {
    dyn.load("gr4j_with_routing_error/src/gr4j_with_routing_error.dll")
  }
}

set.seed(12321)

length_ts<-40
true_input_trial<-runif(length_ts,0,20) #runif(length_ts,100,112)
E_input<-runif(length_ts,0,20) #runif(length_ts,100,105)
input_error_sd<-1e-1 #1e-1
input_error<-rnorm(length(true_input_trial),mean=0,sd=input_error_sd) #sd=0.5
input_trial<-true_input_trial+input_error
input_trial<-pmax(input_trial,0)
actual_state_error<-rnorm(length(input_trial)-1,mean=0,sd=1) #1e-30
actual_state_error_R<-rnorm(length(input_trial)-1,mean=0,sd=1) #1e-30
discharge_error_sd<-0.01 #0.01
actual_discharge_error<-rnorm(length(input_trial),mean=0,sd=discharge_error_sd) #sd=1e-6
true_input<-data.frame(P=true_input_trial,E=E_input)
actual_model_param<-1000
actual_initial_state<-500

param<-c(actual_model_param,1,200,2)
#actual_initial_state<--2000
actual_state_error<-rep(0,length(actual_state_error))
actual_state_error_R<-rep(0,length(rep(0,length(actual_state_error))))

initial_condition_UH1<-rep(0,ceiling(param[4])-1)
initial_condition_UH2<-rep(0,floor(param[4])+length(ceiling(param[4]):ceiling(param[4]*2))-1)


tt1<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
              state_error=actual_state_error, input=true_input,
              initial_condition_UH1=initial_condition_UH1,
              initial_condition_UH2=initial_condition_UH2,
              return_state=F,
              state_error_R=actual_state_error_R,run_compiled=T)

tt2<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
              state_error=actual_state_error, input=true_input,
              initial_condition_UH1=initial_condition_UH1,
              initial_condition_UH2=initial_condition_UH2,
              return_state=T,
              state_error_R=actual_state_error_R,run_compiled=T)

tt1
tt2$output_flow

set.seed(1)
state_error_R_<-rnorm(nrow(true_input)-1,0,1)
tt3<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
              state_error=actual_state_error, input=true_input,
              initial_condition_UH1=initial_condition_UH1,
              initial_condition_UH2=initial_condition_UH2,
              return_state=F,
              state_error_R=state_error_R_,run_compiled=T)
tt4<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
              state_error=actual_state_error, input=true_input,
              initial_condition_UH1=initial_condition_UH1,
              initial_condition_UH2=initial_condition_UH2,
              return_state=F,
              state_error_R=state_error_R_,run_compiled=F)

tt5<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
              state_error=actual_state_error, input=true_input,
              initial_condition_UH1=initial_condition_UH1,
              initial_condition_UH2=initial_condition_UH2,
              return_state=T,
              state_error_R=state_error_R_,run_compiled=T)

tt6<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
              state_error=actual_state_error, input=true_input,
              initial_condition_UH1=initial_condition_UH1,
              initial_condition_UH2=initial_condition_UH2,
              return_state=T,
              state_error_R=state_error_R_,run_compiled=F)


tt5$output_flow
tt6$output_flow

new_state_R<-tt5$output_states_nonUH[-nrow(tt5$output_states_nonUH),6]+state_error_R_

tt7<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
              state_error=actual_state_error, input=true_input,
              initial_condition_UH1=initial_condition_UH1,
              initial_condition_UH2=initial_condition_UH2,
              return_state=F,
              state_error_R=rep(0,length(state_error_R_)),run_compiled=F,
              state_R_ts = new_state_R)
tt8<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
              state_error=actual_state_error, input=true_input,
              initial_condition_UH1=initial_condition_UH1,
              initial_condition_UH2=initial_condition_UH2,
              return_state=F,
              state_error_R=rep(0,length(state_error_R_)),run_compiled=T,
              state_R_ts = new_state_R)
tt9<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
              state_error=actual_state_error, input=true_input,
              initial_condition_UH1=initial_condition_UH1,
              initial_condition_UH2=initial_condition_UH2,
              return_state=T,
              state_error_R=rep(0,length(state_error_R_)),run_compiled=T,
              state_R_ts = new_state_R)
tt10<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
              state_error=actual_state_error, input=true_input,
              initial_condition_UH1=initial_condition_UH1,
              initial_condition_UH2=initial_condition_UH2,
              return_state=T,
              state_error_R=rep(0,length(state_error_R_)),run_compiled=F,
              state_R_ts = new_state_R)

tt3
tt4
tt5$output_flow
tt6$output_flow
tt7
tt8
tt9$output_flow
tt10$output_flow

stop()

timer<-system.time({
  for(i in 1:1000){
    tt1<-gr4j.run(param=param, initial_state_S=actual_initial_state, initial_state_R=100, 
                  state_error=actual_state_error, input=true_input, return_state=F, state_error_R=actual_state_error_R)
  }
})

cat("elapsed",timer[3])


#tt2<-gr4j.sma(actual_model_param,actual_initial_state,actual_state_error,true_input,run_compiled = F,return_state = F)

#matplot(cbind(tt1,tt2),type="l")

#tt1-tt2

test_hydromad<-function(){
  set.seed(12321)
  
  length_ts<-40
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
  
  param<-c(actual_model_param,1,200,2)
  
  
  library(hydromad)
  
  tt1<-gr4j.sim(true_input,actual_model_param,S_0=actual_initial_state/actual_model_param)
  
  tt2<-gr4jrouting.sim(tt1,param[2],param[3],param[4],R_0=100/param[3])
  
}


