batch.file.initial.lines<-c(
  "#!/bin/bash",
  "#SBATCH --job-name=pueuc",
  "#SBATCH --time=72:00:00",
  "#SBATCH --nodes=1",
  #"#SBATCH --ntasks-per-node=20",
  "#SBATCH --cpus-per-task=1",
  "#SBATCH --mem=2500",
  "",
  "module load R/3.4.0"
)


if(Sys.info()["sysname"]=="Linux"){
  wd<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2"
} else {
  # wd<-"C:/Users/kim079/Documents/model_optimisation_framework"
  wd<-"//gpfs2-cbr.san.csiro.au/lw_tvd_mdba_work/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2"
}
remove_for_linux<-"//pearceydata.csiro.au|//pearceyflush1.csiro.au|//pearceyflush2.csiro.au"

replace_this<-"//gpfs2-cbr.san.csiro.au/lw_tvd_mdba_work" #"//lw-osm-03-cdc.it.csiro.au/OSM_CBR_LW_TVD_MDBA_work"
replace_with<-"/datasets/work/LW_TVD_MDBA_WORK" #"/OSM/CBR/LW_TVD_MDBA/work"

preprocess_dir<-"output/gr4j.calib.param.state.all.sites.preprocess"
# output_dir<-"output/gibbs_sampler_param_uncertainty_on_state_errors_real_data_all_sites_v2"
batch.write.dir<-"scripts/PUE_generic_sites_seeds_upcov"

setwd(wd)
simulation_data_files<-list.files(preprocess_dir,pattern = "state_error_simulation_data_",full.names = T)

all_ids<-gsub("state_error_simulation_data_|.csv","",basename(simulation_data_files))
remove_sites<-"110014|415202|415214"
all_ids<-all_ids[grep(remove_sites,all_ids,invert = T)]

# dir.create(output_dir, showWarnings = F)
dir.create(batch.write.dir, showWarnings = F)

all_batch_names<-c()
for(ii in 1:4){
  for(i in 1:length(all_ids)){
    output_dir<-paste0("output/PUE_generic_sites_upcov_seeds_",ii)
    command<-rep(paste0("Rscript /datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2/scripts/PUE_generic_sites_seeds_upcov.r ",
                   "\"",all_ids[i],"\" ",output_dir),11)
    batch.file.lines<-c(batch.file.initial.lines,command)
    batch_fn<-paste0(batch.write.dir,"/",all_ids[i],"_seed_",ii)
    all_batch_names<-c(all_batch_names,basename(batch_fn))
    batch.file.lines<-gsub(remove_for_linux,"",batch.file.lines)
    batch.file.lines<-gsub(replace_this,replace_with,batch.file.lines)
    writeLines(batch.file.lines,batch_fn)
  }
}

batch_runner_fn<-paste0(batch.write.dir,"/batch_runner.sh")
all_cluster_lines<-c(paste("dos2unix",all_batch_names),paste("sbatch",all_batch_names))
writeLines(all_cluster_lines,batch_runner_fn)

# batch_cancel<-paste("scancel -u kim079 -n",all_job_names)
# writeLines(batch_cancel,paste0(batch.file.write.dir,"/kill_all_jobs"))
