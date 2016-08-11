
# Author: tim
###############################################################################
library(foreign)
library(lubridate)
Hm <- 
		read.spss("/home/tim/Data/Quebec/FrenchCanadian.individuals.2012-01-27/RPQA.MarcKlemp.individus.2012-01-27.sav")
Hm <- as.data.frame(Hm)
range(Hm$dateNaissAnnee, na.rm=TRUE)



Hm$BD   						<- as.Date(with(Hm, paste(dateNaissAnnee, 
													dateNaissMois, 
													dateNaissJour, 
													sep = "-")))
DD 								<- with(Hm, paste(	dateDecesAnnee, 
													dateDecesMois, 
													dateDecesJour, 
													sep = "-"))
DD[is.na(Hm$dateDecesAnnee)] 	<- NA
Hm$DD 							<- as.Date(DD)

Hm$L <- decimal_date(Hm$DD) - decimal_date(Hm$BD)
hist(Hm$L[Hm$L > 1])

# OK we have idPere and idMere for parent linkage.
# let's find everyone whose mother is linked.

HasMother <- Hm[!is.na(Hm$idMere), ]
moms      <- HasMother$idMere
# most moms are in the data!
length(unique(moms)) - 
sum(unique(moms) %in% Hm$idIndividu)
moms <- unique(moms)
momsin <- moms[moms %in% Hm$idIndividu]
# moms in data that have a death date.
momsDD <- momsin[!is.na(Hm[Hm$idIndividu %in% momsin,"DD"])]

# lots of mothers w known death dates!
length(momsDD)

# 
head(Hm)
CohDD  <- table(Hm$dateNaissAnnee[Hm$idIndividu %in% momsDD])
Cohmom <- table(Hm$dateNaissAnnee[Hm$idIndividu %in% moms])

plot(as.integer(names(CohDD)),CohDD/Cohmom[names(CohDD)], col = "#00000080", pch = 16)
grid()
axis(2)
abline(v=1735)
abline(h=.9,col="red")

# sooo, 20-year cohorts?
# 1620-39, 1640-59, 1660-79, 1680-99, 1700-19, 1720-39 ??

# OK, cut down to just mothers from these cohorts w known death dates.
MomCD <- Hm[Hm$idIndividu%in% momsDD & Hm$dateNaissAnnee >= 1620 & Hm$dateNaissAnnee <= 1739,]
MomCD$CohBig <- MomCD$dateNaissAnnee - MomCD$dateNaissAnnee %% 20

# might have to start in 1660? will see.
#table(MomCD$CohBig)
head(MomCD)
LS <- table(MomCD$CohBig,round(MomCD$L))
LS <- LS[rowSums(LS)>1000,]
ages <- as.integer(colnames(LS))
LS/rowSums(LS)
matplot(ages,t(LS/rowSums(LS)), type = 'l')




hist(Hm$dateNaissAnnee[Hm$idIndividu %in% momsDD])
hist(Hm$dateNaissAnnee[Hm$idIndividu %in% moms], add = TRUE)




#png("/home/tim/git/ThanoRepro/ThanoRepro/Figures/QuebecL.png",width=12000,height=4500)
#par(xaxs="i", yaxs="i",mai=c(.8,.8,.2,.2))
#plot(NULL, xlim=c(1572,1870),ylim=c(0,110),asp=1,xlab="Year",ylab="Age")
#segments(decimal_date(Hm$BD),0,decimal_date(Hm$DD),Hm$L,col = "#00000020",lwd=.2)
#dev.off()
