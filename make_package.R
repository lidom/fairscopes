## ################################
##
## Make Package Codes
##
###################################

## Remove pkg 
remove.packages("fairscopes")

## Create/update documentation and (re-)write NAMESPACE
devtools::document()

## CRAN-check pkg
#devtools::check(remote = TRUE)  
#devtools::check_built(path = "../fairscopes", remote = TRUE)  

## Install
devtools::install_local(force = TRUE)
##
library("fairscopes")
citation("fairscopes")

help("fairscopes")
