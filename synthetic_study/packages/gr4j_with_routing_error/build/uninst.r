#' Uninstall package function
#'
#' Return warnings if unsuccesful
#' @param pcknm Package name
#' @return Issue warnings if unsuccessful
uninst<-function(pcknm){

	cat(sprintf("\n... Removing package %s ..\n",pcknm))

	# Detach package if loaded
	if(length(which(search()==paste("package:",pcknm,sep="")))>0) 
		eval(parse(text=paste("detach(package:",pcknm,")",sep="")))

	# Look for installed packages
	ipck <- which(installed.packages()[,1]==pcknm)
	if(length(ipck)>0) for(i in ipck){
		err<- try(remove.packages(pcknm,lib=installed.packages()[i,2]),silen=TRUE)
		warnings(as.character(err))
	}
}

