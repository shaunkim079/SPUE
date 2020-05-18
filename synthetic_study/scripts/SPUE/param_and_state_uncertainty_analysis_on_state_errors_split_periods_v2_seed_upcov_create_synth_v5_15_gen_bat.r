batch.file.initial.lines<-c(
  "#!/bin/bash",
  "#SBATCH --job-name=SPcrsy15",
  "#SBATCH --time=72:00:00",
  "#SBATCH --nodes=1",
  #"#SBATCH --ntasks-per-node=20",
  "#SBATCH --cpus-per-task=1",
  "#SBATCH --mem=5000",
  "",
  "module load R/3.4.0"
)
# #SBATCH --time= 96:00:00


if(Sys.info()["sysname"]=="Linux"){
  wd<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2"
  rscript_runner<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2/scripts/param_and_state_uncertainty_analysis_on_state_errors_split_periods_v2_seed_upcov_create_synth_v5_15.r"
} else {
  # wd<-"C:/Users/kim079/Documents/model_optimisation_framework"
  wd<-"Y:/"
  # rscript_runner<-"C:/Users/kim079/Documents/model_optimisation_framework/scripts/Gibbs_state_uncertainty_bates_routing_error_adjusted_likelihood_runner_v10.r"
  rscript_runner<-"Y:/scripts/param_and_state_uncertainty_analysis_on_state_errors_split_periods_v2_seed_upcov_create_synth_v5_15.r"

}

remove_for_linux<-"//pearceydata.csiro.au|//pearceyflush1.csiro.au|//pearceyflush2.csiro.au"

# replace_this<-"//gpfs2-cbr.san.csiro.au/lw_tvd_mdba_work" #"//lw-osm-03-cdc.it.csiro.au/OSM_CBR_LW_TVD_MDBA_work"
# replace_with<-"/datasets/work/LW_TVD_MDBA_WORK" #"/OSM/CBR/LW_TVD_MDBA/work"
replace_this<-"Y:"
replace_with<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2"

batch.write.dir<-"scripts/syn_SPUE_create_synth_15"




setwd(wd)

# dir.create(output_dir, showWarnings = F)
dir.create(batch.write.dir, showWarnings = F)

all_seeds<-1:4

all_batch_names<-c()
for(iseeds in all_seeds){
  command<-rep(paste("Rscript",rscript_runner,iseeds),11)
  batch.file.lines<-c(batch.file.initial.lines,command)
  batch_fn<-paste0(batch.write.dir,"/",iseeds)
  all_batch_names<-c(all_batch_names,basename(batch_fn))
  batch.file.lines<-gsub(remove_for_linux,"",batch.file.lines)
  batch.file.lines<-gsub(replace_this,replace_with,batch.file.lines)
  writeLines(batch.file.lines,batch_fn)
  
}





batch_runner_fn<-paste0(batch.write.dir,"/batch_runner.sh")
all_cluster_lines<-c(paste("dos2unix",all_batch_names),paste("sbatch",all_batch_names))
writeLines(all_cluster_lines,batch_runner_fn)

