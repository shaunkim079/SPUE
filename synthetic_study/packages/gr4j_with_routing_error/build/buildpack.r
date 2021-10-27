library(devtools)
library(testthat)
library(roxygen2)
library(RCurl)


if(Sys.info()["sysname"]=="Linux"){
  # setwd("/data/kim079/model_optimisation_framework_v2/packages")
  setwd("/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2/packages")
} else {
  setwd("C:/Users/kim079/Documents/model_optimisation_framework/packages")
}

# uninstall function
source("gr4j_with_routing_error/build/uninst.r")

if(is.loaded("sma_gr4j_sk")){
  if(Sys.info()["sysname"]=="Linux"){
    dyn.unload("gr4j_with_routing_error/src/gr4j_with_routing_error.so")
  } else {
    dyn.unload("gr4j_with_routing_error/src/gr4j_with_routing_error.dll")
  }
}


# Function used to compile the C and Fortran code
compile <-function(pcknm){
  path <- paste(pcknm,"/src",sep="")
  lf <- list.files(path,full.names=TRUE)
  unlink(lf[grep("\\.o$|\\.dll$|\\.so$",lf)],force=TRUE)
  cffiles <- list.files(path,pattern=".*\\.c$|.*\\.f$|.*\\.f90$")
  
  ext <- ".so "                                  
  if(Sys.info()["sysname"]!="Linux") ext <- ".dll "
  
#   cmd <- paste("R CMD SHLIB -o ",path,"/",pcknm,ext,
#                paste(path,"/",cffiles,sep="",collapse=" "),
#                " 2>../scripts/compile.log",collapse="",sep="")
  cmd <- paste("R CMD SHLIB -o ",path,"/",pcknm,ext,
               paste(path,"/",cffiles,sep="",collapse=" "),
               collapse="",sep="")
  
  system(cmd)
  return(cmd)
}


pcknm <- "gr4j_with_routing_error"

# Detach and uninstall package if already installe
uninst(pcknm)

compile(pcknm)
# document(pcknm)

lfo <- list.files(pcknm,patter="\\.o$",recursive=TRUE,full.names=TRUE)
unlink(lfo)

# Build  (binaries in windows and sources in :inux)
# build(pcknm,path="../binaries",binary=TRUE) # final build
