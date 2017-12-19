
#' A unified loader for github and CRAN packages
#'
#' This function installs, updates and loads packages provided by CRAN or github
#'
#'@param required.packages a vector of characters
#'@param mode character defining wether the packages are provided by CRAN or github, use "CRAN" or "GIT"
#'@param update logical wether to force updating previously installed packages
#'
#' @return none
#'
#' @keywords packages, load
#'
#' @examples
#' load.and.install(c("ggplot2", "raster"), mode = "CRAN", update = T)
#'
#'@export


load.and.install <- function(required.packages, mode = "CRAN", update = F){
  for(lib.file in .libPaths()) assign(lib.file, list.files(lib.file))
  rm(lib.file)
  all.packages <- do.call(c, mget(ls()))
  
  if(mode == "CRAN"){
    for(package in required.packages){
      if(!package%in%all.packages | update) install.packages(package)
      library(package, character.only = T)
    }
  }
  
  if(mode == "GIT"){
    
    if(!'devtools'%in%all.packages) install.packages('devtools')
    library(devtools)
    
    package.names <- unlist(lapply(strsplit(required.packages, "/"), function(x) tail(x, n = 1)))
    
    for(package in 1:length(required.packages)){
      if(!package.names[package]%in%all.packages | update) install.packages(required.packages[package])
      library(package.names[package], character.only = T)
    }
  }
}
