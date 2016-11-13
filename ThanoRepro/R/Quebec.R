
# Author: tim
###############################################################################
library(foreign)
library(lubridate)
setwd("/home/tim/Data/Quebec/")

# this is the MPIDR version from 2012, currently checking to see if it has been updated.
path      <- "FrenchCanadian.individuals.2012-01-27/RPQA.MarcKlemp.individus.2012-01-27.sav"
# TR: line to read file
Q         <- read.spss(path)
# TR: line to read file on MG's system
# Hm <-   read.spss("U:/quebec/Quebec/Quebec/FrenchCanadian.individuals.2012-01-27/RPQA.MarcKlemp.individus.2012-01-27.sav")
Q         <- as.data.frame(Q)

# ---------------------------------------#
# begin data prep                        #
# ---------------------------------------#

# compose Birthday
Q$BD   	  <- as.Date(with(Q, paste(dateNaissAnnee, 
						              dateNaissMois, 
						              dateNaissJour, 
						              sep = "-")))
#compose Deathday, need to hangle NAs
DD 		  <- with(Q, paste(dateDecesAnnee, 
				         dateDecesMois, 
				         dateDecesJour, 
				         sep = "-"))
		 
# missing days/months appear to have been imputed, good!
#any(is.na(Q$dateDecesMois) & !is.na(Q$dateDecesAnnee))
#any(is.na(Q$dateDecesJour) & !is.na(Q$dateDecesMois))
DD[is.na(Q$dateDecesAnnee)] 	<- NA
Q$DD 							<- as.Date(DD)

# calculate decimal lifespan
Q$L       <- decimal_date(Q$DD) - decimal_date(Q$BD)

#hist(Hm$L[Hm$L > 1])

# OK we have idPere and idMere for parent linkage.
# let's find everyone whose mother is linked.
# idPere would have been cool for my dissertation
# 4 years ago...

# everyone needs a birthdate for this analysis
# (but not everyone needs a death date)
Q         <- Q[!is.na(Q$BD),]

# specifically, mothers need both birth and death dates,
# but children only need birth dates, since we're measuring
# the fertility of potential mothers. So, for cleanliness, let's split
# Q into two, M, for mothers, and Ch, for children.
M         <- Ch      <- Q

# filter potential mothers down to females with non-NA L:
M         <- M[!is.na(M$L) & M$sexe == "f", ]

# children need to be filtered down to only those
# whose mother is still in M:
Ch        <- Ch[Ch$idMere %in% M$idIndividu, ]

# now Ch needs two new columns, MBD, and MDD, 
# mothers birth and death dates, respectively.

# make named vectors 
# (not efficient, but intuitive, and avoids loop, which is less efficient)
BD        <- M$BD
names(BD) <- M$idIndividu
DD        <- M$DD
names(DD) <- M$idIndividu
# we're using name-selection here. 
Ch$MBD    <- BD[as.character(Ch$idMere)]
Ch$MDD    <- DD[as.character(Ch$idMere)]
#any(is.na(Ch$MBD)) # works
#any(is.na(Ch$MDD))

# now we have anough information to calculate
# mothers age and time-to-death when id is born:
Ch$MT <- decimal_date(Ch$MDD) - decimal_date(Ch$BD) # time to death
Ch$MA <- decimal_date(Ch$BD) - decimal_date(Ch$MBD) # chronological age
#range(Ch$MT)
## how many negative mothers' remaining lives are there?
#sum(Ch$MT < 0)
#mean(Ch$MT < 0)
#sum(Ch$MT < -1)
## this appears to be an issue of date precision and not outright error.
## it looks like some mothers died in pregnancy, and somehow births
## end up getting coded earlier than maternal deaths. Check that Ch$L is
## 0ish for these cases:
#Ch$L[Ch$MT < 0] # indeed, mostly stillbirths. 
# This won't change things
# either way, but it seems better to not throw these cases out,
# and instead assume, Ch$MT == 0 in these cases.
Ch$MT[Ch$MT < 0] <- 0
# that's an executive decision of sorts. Note that it's only a problem because
# we're on a boundary. Other ages and TTD estimates are probably jittered to
# a similar degree, but this appears to be trivial for such aggregate estimates.

#hist(Ch$MT -Ch$MA) # very cool: most have more years left than lived :-)

# ---------------------------------------#
# end data prep                          #
# ---------------------------------------#


# ---------------------------------------#
# a couple functions to calc rates       #
# ---------------------------------------#
replaceNANaN0 <- function(x){
	x[is.na(x) | is.nan(x)] <- 0
	x
}

simple.smooth <- function(ASFR, age){
	TFR                  <- sum(ASFR)
	
	# detect if this is chronological or not:
	mult                 <- rep(1,length(ASFR))
	if (sum(ASFR[1:10]) == 0){
		mult[1:10]       <- 0
	}
		
	ASFRhat              <- predict(loess(ASFR ~ age, span = .25),
			                        newdata = data.frame(age = age))
	# don't allow negatives
	ASFRhat[ASFRhat < 0] <- 0
	ASFRhat              <- ASFRhat * mult
	# rescale to force sum
	ASFRhat * (TFR / sum(ASFRhat))
}

get.birth.rates    <- function(
		Ch,          # data.frame of potential children
		M,           # data.frame of potential mothers
		startyear,   # lower bound of cohort
		cohN = 5,    # cohort width (5 default)
		ageN = 1,    # age interval (1 default)
		omega = 110, # maximum age to fill out to, default 110
		sex = NULL) { # optional, sex of birth.
	# choose the cohort of your choice
	endyear             <- startyear + cohN
	cohi                <- M$dateNaissAnnee >= startyear & M$dateNaissAnnee <= endyear
	cohort              <- M[cohi, ]
	
	# radix is number born
	l0                  <- nrow(cohort)
	
	# calc age:
	age <- seq(0, omega, by = ageN)
	
	# calculate Lx  (denominator)
	denom               <- rep(0, length(age))
	for (i in 1:length(age))  {
		# ceiling() will work for single, but need fancy otherwise
		ind             <- (cohort$L - (cohort$L %% ageN) + ageN) == age[i]
		# one more step with %% to get right age group.
		dec             <- sum(ageN - (age[i] - cohort$L[ind]))
		denom[i]        <- sum(cohort$L >= age[i]) + dec
	} 
	
	# optionally, just select female births:
	if (!is.null(sex)){
		Ch <- Ch[Ch$sexe == sex, ]
	}
	# get the number of children, which were born by moms of this cohort
		
	chi                 <- Ch$idMere %in% cohort$idIndividu
	children            <- Ch[chi, ]
	numerator           <- table(floor(children$MT), floor(children$MA))
	# image(as.integer(colnames(numerator)),as.integer(rownames(numerator)),t(numerator),
	# asp=1, xlab = "years lived", ylab = "years left")
	
	# get birth counts crosstab
	BA                  <- colSums(numerator)
	BT                  <- rowSums(numerator)
	
	out                 <- data.frame(Cohort = startyear,
			                          Age = age, 
			                          Exposure = denom, 
									  BA = NA, 
									  BT = NA, 
									  ASFR = 0, 
									  TSFR = 0,
									  Lx = denom / l0)
	# assign birth counts to right ages
	out$BT              <- BT[as.character(out$Age)]
	out$BA              <- BA[as.character(out$Age)]
	out$BT              <- replaceNANaN0(out$BT)
	out$BA              <- replaceNANaN0(out$BA)
	
	# get rates
	out$ASFR            <- out$BA / out$Exposure
	out$TSFR            <- out$BT / out$Exposure
	out$ASFR            <- replaceNANaN0(out$ASFR)
	out$TSFR            <- replaceNANaN0(out$TSFR)
	
	# do some simple smoothing
	out$ASFRs           <- simple.smooth(out$ASFR,age)
	out$TSFRs           <- simple.smooth(out$TSFR,age)
	
	return(out) 
}

# ---------------------------------------#
# end function defs                      #
# ---------------------------------------#

# ---------------------------------------#
# begin calc fertility rates             #
# ---------------------------------------#


#set lower bounds for cohorts
coh5 <- seq(from = 1665, to = 1740, by = 5)  

#like usual... basic represents the first examined cohort
#basic2 is getting added over and over again to create the final data.frame which is still called "basic"
#basic <- get.birth.rates(Ch=Ch,M=M, 1660, 5,110,"f")

Rates1x5 <- do.call(rbind,
		lapply(coh5, function(cohi, Ch, M){
					get.birth.rates(Ch = Ch,
							M = M, 
							startyear = cohi, 
							cohN = 5,
							ageN = 1,
							omega = 110,
							sex = "f")
				},Ch = Ch , M = M))


save(Rates1x5,file = "/home/tim/git/ThanoRepro/ThanoRepro/Data/QuebecRates1x5.Rdata")

Rates2x5 <- do.call(rbind,
		lapply(coh5, function(cohi, Ch, M){
					get.birth.rates(Ch = Ch,
							M = M, 
							startyear = cohi, 
							cohN = 5,
							ageN = 2,
							omega = 110,
							sex = "f")
				},Ch = Ch , M = M))


save(Rates2x5,file = "/home/tim/git/ThanoRepro/ThanoRepro/Data/QuebecRates2x5.Rdata")

Rates3x5 <- do.call(rbind,
		lapply(coh5, function(cohi, Ch, M){
					get.birth.rates(Ch = Ch,
							M = M, 
							startyear = cohi, 
							cohN = 5,
							ageN = 3,
							omega = 110,
							sex = "f")
				},Ch = Ch , M = M))


save(Rates3x5,file = "/home/tim/git/ThanoRepro/ThanoRepro/Data/QuebecRates3x5.Rdata")

Rates5x5 <- do.call(rbind,
		lapply(coh5, function(cohi, Ch, M){
					get.birth.rates(Ch = Ch,
							M = M, 
							startyear = cohi, 
							cohN = 5,
							ageN = 5,
							omega = 110,
							sex = "f")
				},Ch = Ch , M = M))


save(Rates5x5,file = "/home/tim/git/ThanoRepro/ThanoRepro/Data/QuebecRates5x5.Rdata")
coh10 <- seq(from = 1660, to = 1740, by = 10) 
Rates1x10 <- do.call(rbind,
		lapply(coh, function(cohi, Ch, M){
					get.birth.rates(Ch = Ch,
							M = M, 
							startyear = cohi, 
							cohN = 10,
							omega = 110,
							sex = "f")
				},Ch = Ch , M = M))


save(Rates1x10,file = "/home/tim/git/ThanoRepro/ThanoRepro/Data/QuebecRates1x10.Rdata")
Rates2x10 <- do.call(rbind,
		lapply(coh, function(cohi, Ch, M){
					get.birth.rates(Ch = Ch,
							M = M, 
							startyear = cohi, 
							cohN = 10,
							ageN = 2,
							omega = 110,
							sex = "f")
				},Ch = Ch , M = M))


save(Rates2x10,file = "/home/tim/git/ThanoRepro/ThanoRepro/Data/QuebecRates2x10.Rdata")
# ----------------------------------------------#
# experiment smoothing                          #
# ----------------------------------------------#
#library(reshape2)
#library(RColorBrewer)
#Rates <- local(get(load("/home/tim/git/ThanoRepro/ThanoRepro/Data/QuebecRates1x5.Rdata")))
#head(Rates)
#
#ASFR <- acast(Rates, Age~Cohort, value.var = "ASFRs")
#TSFR <- acast(Rates, Age~Cohort, value.var = "TSFR")
#
#cols <- colorRampPalette(brewer.pal(9,"Blues"),space="Lab")(ncol(TSFR))
#png("/home/tim/git/ThanoRepro/ThanoRepro/QASFR.png",width=800,height=400)
#par(mfrow=c(1,2))
#matplot(0:90,ASFR[1:91, ], type = 'l', lty = 1, col = cols, lwd = seq(3,1,length=16),
#		xlab = "years left", ylab = "fertility rate", main = "ASFR")
#matplot(0:90,TSFR[1:91, ], type = 'l', lty = 1, col = cols, lwd = seq(3,1,length=16),
#		xlab = "years left", ylab = "fertility rate", main = "TSFR")
#dev.off()
#
#LX <- acast(Rates, Age~Cohort, value.var = "Lx")
#
#png("/home/tim/git/ThanoRepro/ThanoRepro/Lx.png")
#matplot(0:90,LX[1:91, ], type = 'l', lty = 1, col = cols, lwd = seq(3,1,length=16),
#		xlab = "years left", ylab = "survivors", main = "Survivorships", ylim=c(0,1))
#dev.off()
#
#LX[16:91, ] %*% diag(1/LX[16,])
#png("/home/tim/git/ThanoRepro/ThanoRepro/Lx.png")
#matplot(15:90,LX[16:91, ] %*% diag(1/LX[16,]), type = 'l', lty = 1, col = cols, lwd = seq(3,1,length=16),
#		xlab = "years left", ylab = "survivors", main = "Survivorships", ylim=c(0,1))
#dev.off()

#TSFR <- acast(Rates2x5, Age~Cohort, value.var = "TSFR")
#matplot(as.integer(rownames(TSFR)),TSFR[,-c(1,2)],lty=1,type='l')
#TSFR <- acast(Rates2x10, Age~Cohort, value.var = "TSFR")
#matplot(as.integer(rownames(TSFR)),TSFR[,-c(1,2)],lty=1,type='l')
#
#TSFR <- acast(Rates1x5, Age~Cohort, value.var = "TSFR")
#TSFRs <- TSFR*0
#for (i in 1:ncol(TSFR)){
#	TSFRs[,i] <- simple.smooth(TSFR[,i],age)
#    plot(age, TSFR[,i],main=colnames(TSFR)[i], xlab = "T", ylab = "TSFR")
#    lines(age, TSFRs[,i])
#    locator(1)
#}

#simple.smooth <- function(ASFR, age){
#	TFR     <- sum(ASFR)
#	ASFRhat <- smooth.spline(ASFR~age)$y
#	# don't allow negatives
#	ASFRhat[ASFRhat < 0] <- 0
#	# rescale to force sum
#	ASFRhat * (TFR / sum(ASFRhat))
#}



#image(as.integer(colnames(TSFR)),age,t(TSFR),asp=1)
#dev.new()
#image(as.integer(colnames(TSFRs)),age,t(TSFRs),asp=1)
#
#matplot(age, TSFR, type='l', col = "#00000040", lty =1)
#matplot(age, TSFRs, type='l', col = "#0000FFAA", lty =1,add=TRUE)
#


#
#for (i in 1:ncol(BT)){
##	TSFR2[1:91,i] <- Mort1Dsmooth(x = ages[1:91],
##			y = BT[1:91,i], offset = EX[1:91,i], w = W[1:91, i])
#}

#head(Rates)
#library(reshape2)
#matplot(0:110,acast(Rates,Age~Cohort,value.var="ASFR"),type='l',lty=1,col="#00000050")
#matplot(0:110,acast(Rates,Age~Cohort,value.var="TSFR"),type='l',lty=1,col="#00000050")
#matplot(0:110,acast(Rates,Age~Cohort,value.var="Lx"),type='l',lty=1,col="#00000050")
#
## suspect immigration survival bias.
## would need to filter out immigrants.
## check how many.
#LX <- acast(Rates,Age~Cohort,value.var="Lx")
#plot(colSums(LX))
#
## but it could also be entirely due to fertility,
## higher fert means more infant deaths which means
## lower e0, by a lot. Test e10.
#plot(rowSums(t(LX[11:111, ]) / LX[11,]))
## appears stationary after 8th cohort
#coh[8]
# should probably pick things up in 1665 then...

# ---------------------------------------#
# end calc fertility rates               #
# ---------------------------------------#

#-----------------------------------------
# the rest of this is chunks from earlier versions
# and attempts at doing what the above code does
# files names may have changed, etc. There are 
# also some diagnostics that led to choosing the 
# bounds for which cohorts to include here and not.
# from here until the bottom of the script is
# no-man's land, enter at your own peril
#-----------------------------------------


# TR: make this more efficient
#for (i in 1:length(lower)) {
#	basic2 <-get.birth.rates(number.children, Q.with.L, lower[i], upper[i])
#	basic  <-rbind(basic, basic2)
#}


#save(basic, file = "U:/quebuec2/quebuec_rates.Rdata")
# test<-basic[basic$Cohort=="1755-1759",]
# sum(test$TSFR)
# plot (test$ASFR)
# 
# sum(test$BT)
# sum(test$BA)
# head(number.children)


#----------------------
# code replaced by TR
## numerators can only come from people
## with mothers that are also in the data
#HasMother <- Q[!is.na(Q$idMere), ]
#moms      <- HasMother$idMere
#
## most moms are in the data!
##length(unique(moms)) - 
##		sum(unique(moms) %in% Hm$idIndividu)
#
## unique mothers
#moms             <- unique(moms)
#
## mothers that are also ego in the data.
## (necessary because we need birth and death dates of mothers)
#momsin           <- moms[moms %in% Q$idIndividu]
#
## get a cut-down version of the data, save this for use later
#all.useful.moms  <- Q[Q$idIndividu %in% momsin,]
#head(all.useful.moms)
## head(all.useful.moms)
## dim(all.useful.moms)
#
## moms in data that have a death date.
#momsDD           <- momsin[!is.na(Q[Q$idIndividu %in% momsin,"DD"])]
#
#
################################################################################################
############## calculating number of births ###########################################################
#
#
## select all useful moms which have a death date and are linked
## check if they are linked
#all.useful.moms  <- Q[Q$idIndividu %in% momsin,]
##check if the death date is available
#all.useful.moms  <- all.useful.moms[!is.na(all.useful.moms$DD), ] 
##these are the moms with death date and linking
#
##we need people with available birthdate
##so lets delete people with no valid birth date
#
##     table(is.na(Hm$BD)) # 55012 people with NA's for birthdate
#Q.with.BD       <- Q[!is.na(Q$BD),]
#
##filter all.useful.moms again
##     table(is.na(all.useful.moms$BD)) #5307 moms with no birthdate - 35498 should remain
##filter by available death dates for moms and individuals 
#all.useful.moms <- all.useful.moms[all.useful.moms$idIndividu %in% Q.with.BD$idIndividu,]
#all.useful.moms <- all.useful.moms[all.useful.moms$idIndividu %in% Q.with.BD$idMere,]
#
##list of ID's from each useful mother
#list.useful.moms<-unique(all.useful.moms$idIndividu)
#
#
#
#######################################################################################################################################
#### now calculate it with all children of each mother#################################################################################
#
##in this step... I'm doing a basic dataframe (as usual :) to use rbind in the next step to build one big file
##so the beginning is called "one.mom" and "one.moms" is getting added over and over again to create "number.children"
##which is used for the birthcounts later in this script
#
#one.mom                       <- Q.with.BD[Hm.with.BD$idMere %in% list.useful.moms[1],]
#one.mom                       <- one.mom[order(one.mom$BD),] # sort from first birth to last birth    
##choose the correct mom 
#right.mom                     <- Q.with.BD[Q.with.BD$idIndividu %in% one.mom$idMere,]
## and get the mom's age at childbearing
#one.mom$childbearing          <- decimal_date(one.mom$BD) - decimal_date(right.mom$BD)
## and the remaining years of life
#one.mom$remaining             <- decimal_date(right.mom$DD) - decimal_date(one.mom$BD)
#
#for (i in 2:length(list.useful.moms)) {
#	one.moms                       <- Q.with.BD[Q.with.BD$idMere %in% list.useful.moms[i],]
#	one.moms                       <- one.moms[order(one.moms$BD),] # sort from first birth to last birth    
#	#choose the correct mom 
#	right.mom                     <- Q.with.BD[Q.with.BD$idIndividu %in% one.moms$idMere,]
#	# and get the mom's age at childbearing
#	one.moms$childbearing          <- decimal_date(one.moms$BD) - decimal_date(right.mom$BD)
#	# and the remaining years of life
#	one.moms$remaining             <- decimal_date(right.mom$DD) - decimal_date(one.moms$BD)
#	# combine the datasets to produce on large file
#	one.mom                       <- rbind(one.mom,one.moms)
#}
## str(number.children)
## head(one.mom)
## dim(one.mom)
#
## round the age at childbearing
#one.mom$childbearing              <-floor(one.mom$childbearing)
#
#
#number.children$Tfloor <- floor(number.children$remaining)
#head(number.children)
#
#
##################################################################################################
##################### calculating exposures #####################################################
#
## we need individual data with complete lifespan values
#Q.with.L                   <- Q[!is.na(Q$L),]
#
## we need all potential mom's -> so just pick the females
#Q.with.L                   <- Q.with.L[Q.with.L$sexe=="f",]
# ---------------------------------------------------
# TR: commented out because not used. Some diagnostics in here

# table(one.mom$childbearing)
# minimum age at childbearing: 12
# maximum age at childbearing: 59

# hist(one.mom$childbearing)

#save the file because the loop takes a while...
#save(one.mom, file = "C:/Users/Markus/Documents/MPI/quebuec/number_children.Rdata")

#load the data
#number.children<-get(load("U:/quebuec2/quebuec/number_children.Rdata"))
 
# dim(Hm.with.L) #we have 105519 potential moms
# table(Hm.with.L$dateNaissAnnee)

# #let's try to calculate fertility rates for differeny cohorts
# #exposures
# cohort.1660.1700       <- Hm.with.L[Hm.with.L$dateNaissAnnee >= 1660 & Hm.with.L$dateNaissAnnee <= 1700,]
# 
# #numerator
# #need to find all children with known moms from these cohorts to calculate correct birthcounts
# 
# children.1660.1700         <-number.children[number.children$idMere %in% cohort.1660.1700$idIndividu,]
# table(children.1660.1700$childbearing)
# 
# numeratorT                                 <- table(children.1660.1700$remaining,children.1660.1700$childbearing)
# blabla<-as.numeric(rownames(numeratorT))
# remaining.years<-rep(NA,ncol(numeratorT))
# for (i in 1:ncol(numeratorT)) {
#   remaining.years[i]<-sum(blabla*numeratorT[,i])
# }
# numeratorT<-data.frame(T=as.numeric(colnames(numeratorT)),remaining=remaining.years)


# rowSums(numeratorT)
# colSums(numeratorT)

#exposures        == Hm.with.L
#startyear        == begin of the cohort (given as number)
#endyear          == end of the cohort (given as number)
#birthcounts      == number.children  

#length(momsDD)
# lots of mothers w known death dates!
#length(momsDD)

# let's see how people in different cohorts are observed from birth until death
#CohDD  <- table(Q$dateNaissAnnee[Q$idIndividu %in% momsDD])
#Cohmom <- table(Q$dateNaissAnnee[Q$idIndividu %in% moms])

# take a look at proportion of birth cohorts (single year) that have
# a death date assigned:
#plot(as.integer(names(CohDD)),CohDD/Cohmom[names(CohDD)], col = "#00000080", pch = 16)
#grid()
#axis(2)
#abline(v=1735)         # a not-bad cutoff
#abline(h=.9,col="red") # 90% of birth cohort has death date

# sooo, 20-year cohorts?
# 1620-39, 1640-59, 1660-79, 1680-99, 1700-19, 1720-39 ??
# nah, we can calc for 5-year cohorts and aggregate as needed

# OK, cut down to just mothers from these cohorts w known death dates.
#MomCD        <- Q[Q$idIndividu%in% momsDD & Q$dateNaissAnnee >= 1620 & Q$dateNaissAnnee <= 1739,]
#MomCD$CohBig <- MomCD$dateNaissAnnee - MomCD$dateNaissAnnee %% 20

# might have to start in 1660? will see.
#table(MomCD$CohBig)
#head(MomCD)
#LS <- table(MomCD$CohBig,round(MomCD$L))
#LS <- LS[rowSums(LS)>1000,]
#ages <- as.integer(colnames(LS))
#LS/rowSums(LS)
#matplot(ages,t(LS/rowSums(LS)), type = 'l')
#
#hist(Q$dateNaissAnnee[Q$idIndividu %in% momsDD])
#hist(Q$dateNaissAnnee[Q$idIndividu %in% moms], add = TRUE)
#png("/home/tim/git/ThanoRepro/ThanoRepro/Figures/QuebecL.png",width=12000,height=4500)
#par(xaxs="i", yaxs="i",mai=c(.8,.8,.2,.2))
#plot(NULL, xlim=c(1572,1870),ylim=c(0,110),asp=1,xlab="Year",ylab="Age")
#segments(decimal_date(Hm$BD),0,decimal_date(Hm$DD),Hm$L,col = "#00000020",lwd=.2)
#dev.off()

# momsCBD<-Hm[momsDD %in% Hm$idIndividu,]

## get mother's remaing years of life by taking difference between death date of mom 
## and birth date of child

# only 69.934 different mom's in this dataset (because some of them have more than one child)

##########################################################################################################################
# MG: it's commented out because it just calculates everything for the first birth
#so never mind
#the real calculation of births is coming next


#     #create a dataframe for remaing years of life (of the mom since the first birth) 
#     final.data<-all.useful.moms
#     final.data$remaining<-NA
#     head(final.data)
# 
#     #Loop for remaining years after the first birth
#     for (i in 1:length(list.useful.moms)) {
#       #select every child of one specific mother
#       one.mom             <- Hm.with.BD[Hm.with.BD$idMere %in% list.useful.moms[i],]
#       one.mom             <- one.mom[order(one.mom$BD),] # sort from first birth to last birth
#       #take the first birth
#       first.birth         <- one.mom[1,]
#       #choose the correct mom and get the death date to calculate the remaining years
#       right.mom           <- Hm.with.BD[Hm.with.BD$idIndividu %in% first.birth$idMere,]
#       remaining.years     <- decimal_date(right.mom$DD) - decimal_date(first.birth$BD)
#       right.mom$remaining <- remaining.years
#       final.data[i,]       <- right.mom
#     }
# 
#     head(final.data)
# 
# 
# 
#   #now get the mother's age at childbearing
#   final.data$childbearing<-NA
#   for (i in 1:length(list.useful.moms)) {
#     one.mom                       <- Hm.with.BD[Hm.with.BD$idMere %in% list.useful.moms[i],]
#     one.mom                       <- one.mom[order(one.mom$BD),] # sort from first birth to last birth    
#     #take the first birth
#     first.birth                   <- one.mom[1,]
#     #choose the correct mom and get the mom's age at childbearing
#     right.mom                     <- Hm.with.BD[Hm.with.BD$idIndividu %in% first.birth$idMere,]
#     final.data[i,]$childbearing   <- decimal_date(first.birth$BD) - decimal_date(right.mom$BD) 
#   }
# 
# #round the mother's age at childbearing to calculate the number of birth's at a certain age
# final.data$childbearing           <-floor(final.data$childbearing)
#    
# #save the file because the loop takes a while...
# save(final.data, file = "U:/quebec/remaining_years_moms.Rdata")
# 
# 
# 
# #load final data 
# 
# first.child<-get(load("U:/quebec/remaining_years_moms.Rdata"))
# head(first.child)

# table(first.child$childbearing)
#minimum age at childbearing: 12 
#maximum age at childbearing: 48


#plot motherŽs age at first child

# hist(first.child$childbearing)
##########################################################################################################################
