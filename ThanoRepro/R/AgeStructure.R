
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
devtools::load_all("R/RiffeFunctions")

# 1) read in example survival curves. How about Sweden 1970, just because?

Data <- local(get(load("Data/DataAll.Rdata")))

yr   <- 2010
Name <- "USA"

# friendly tools
# mx2dxHMD()
# da2fya()
head(Data)
mxm <- with(Data, mx[Sex == "m" & Code == Name & Year == yr])
mxf <- with(Data, mx[Sex == "f" & Code == Name & Year == yr])

# closeout:
malejump <- diff(mxm[c(110,111)])
mxm[111] + 1:10 * malejump
mxme <- c(mxm, mxm[111] + 1:10 * malejump)
femalejump <- diff(mxf[c(110,111)])
mxf[111] + 1:10 * femalejump
mxfe <- c(mxf, mxf[111] + 1:10 * femalejump)
#e0jf <- with(Data, ex[Sex == "f" & Code == "JPN" & Year == 2012 & Age == 0])
#plot(0:110,with(Data, mx[Sex == "f" & Code == "JPN" & Year == 2012 ]),log="y")
#lines(0:110,mxf * coefs[4])



# Note, mxm and mxf were smoothed by HMD into Kannisto. We need a higher closeout
# age for these plots, like 120. Just assume a straight line to asymptote of 1 from last two points.



mx2e0 <- function(mx,sex="m"){
	mx                  <- Mna0(as.numeric(mx))
	# mean proportion of interval passed at death
	ax                  <- mx * 0 + .5                      # ax = .5, pg 38 MPv5
	
	ax[1]               <- AKm02a0(mx[1], sex)              # v6 a0 protocol
	
	ages.mids <- 1:length(mx) - 1 + ax
	
	
	dx <- mx2dxHMD(mx, sex)
	dx <- dx / sum(dx)

	sum(dx * ages.mids)
}

mx2lx <- function(mx, sex){
	dx <- mx2dxHMD(mx, sex)
	rev(cumsum(rev(dx)))
}


# this is the median proportion female of each year/pop of births in the HMD
PF <- 0.486540388060525

getcoef <- function(coef, mx, sex){
	e01 <- round(mx2e0(mx * coef, sex))
	coefout <- c(rep(coef, e01),seq(coef,1,length = length(mx) - e01))
	coefout
}

e0avg <- function(coef, mxm, mxf, PF){
	
	e0m <- mx2e0(mxm * getcoef(coef, mxm, "m"), "m")
	e0f <- mx2e0(mxf * getcoef(coef, mxf, "f"), "f")
	PF * e0f + (1 - PF) * e0m
}

e0min <- function(coef = 1, e0goal = 75, mxm, mxf, PF){
	abs(e0avg(coef,mxm,mxf,PF) - e0goal)
}


# which e0 do we put in columns?
e0vec <- c(55,65,75,85,95)
e0avg(-86, mxm, mxf, PF)

coefs <- c()
for (i in 1:length(e0vec)){
	coefs[i] <-	optimise(
			e0min, 
			interval=c(0,10),
			e0goal=e0vec[i],
			mxm=mxme,
			mxf=mxfe,
			PF=PF, 
			maximum = FALSE, 
			tol = 1e-7)$minimum
		}
e0avg(-10.7133251,mxm,mxf,PF)
plot(0:120, mx2lx(mxme * getcoef(coefs[1], mxme, "m"),"m"), type = 'l', ylim = c(0,1))
lines(0:120, mx2lx(mxme * getcoef(coefs[2], mxme, "m"),"m"))
lines(0:120, mx2lx(mxme * getcoef(coefs[3], mxme, "m"),"m"))
lines(0:120, mx2lx(mxme * getcoef(coefs[4], mxme, "m"),"m"))
lines(0:120, mx2lx(mxme * getcoef(coefs[5], mxme, "m"),"m"))
# so, coefs are what we need to use to derive the right avg e0, just multiply in mxm, mxf
# to get back corresponding lx, rescale using PF and make pyramids.
age <- 0:120
rs <- c(-.02,-.01,0,.01,.02)
par(mai = c(0,0,0,0))

# ready!!!
plot(NULL, type = "n", axes = FALSE, xlab = "", ylab = "",xlim=c(-2,2),ylim=c(0,115))
PyramidOutline(males = mx2lx(mxme*getcoef(coefs[1], mxme, "f"),"m") * (1 - PF) * exp(-age*.02),
		females = mx2lx(mxfe*getcoef(coefs[1], mxfe, "f"),"f") * PF * exp(-age*.02),
		scale = 100,
		col = "gray",
		border =  NA)




