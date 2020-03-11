for(seed_num in 1:4){
  seed_number<-seed_num
  for(ss in 1:11){
    source('~/model_optimisation_framework/scripts/parameter_uncertainty_analysis_on_state_errors_v2_seed_upcov.r')
  }
  
  rm(list = ls())
}
