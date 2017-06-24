pkg_list = pkg_list =c("XML","tidyverse","magrittr","corrplot","tidyr","DT","tibble","stringr","ggplot2","cowplot")

#----Dependency checking function----#
check_pkg <- function(pkg){
  if(!require(pkg)){
    install.packages(pkg)
  }
}

sapply(pkg_list,check_pkg)