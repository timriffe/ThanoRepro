
# Author: tim
###############################################################################

# for Tim, this will choke
if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/ThanoRepro/ThanoRepro")
} else {
	# in that case I'm on Berkeley system, and other people in the dept can run this too
	setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/ThanoRepro/ThanoRepro"))
}

library(HMDHFDplus)
source("R/Functions.R")

# 1) read in example survival curves. How about Sweden 1970, just because?

flt <- readHMDweb("SWE","fltper_1x1",password=pw,username=us)
mlt <- readHMDweb("SWE","mltper_1x1",password=pw,username=us)





