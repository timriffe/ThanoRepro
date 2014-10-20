
# for Tim, this will choke
if (system("hostname",intern=TRUE)=="triffe-N80Vm"){
  # if I'm on the laptop
  setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
  # in that case I'm on Berkeley system, and other people in the dept can run this too
  setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
}


# use devtools package for easy documentation and building.

library(devtools)

devtools::document("R/RiffeFunctions")
#devtools::load_all("R/RiffeFunctions")

library("TimUtils")
IncrementVersion("R/RiffeFunctions","1",package.origin="2014-01-01")
