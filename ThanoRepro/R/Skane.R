
# Author: tim
###############################################################################

data.folder <- "/home/tim/Data/SEDD_4_0_Public"
metadata <- read.table(file.path(data.folder,"VariableSetup.txt"),header=TRUE,sep="\t")
head(metadata)
Dat <- read.table(file.path(data.folder,"SEDD_4_0_Public.txt"),header=TRUE,sep="\t")
nrow(Dat)
table(Dat$Type)
length(unique(Dat$ID_I))
Dat$DayFrac <- NULL
Dat$Date    <- as.Date(paste(Dat$Year, Dat$Month, Dat$Day, sep = "-"))
Dat$Month   <- NULL
Dat$Day     <- NULL

# keep:
keep <- with(Dat, Type %in% 
				c("Sex", "DeathDate", "BirthDate", "FamilyID", 
						"FamilyNumber", "HouseholdID", "HouseholdNo", "BirthChildIndicator"))
Dat <- Dat[keep, ]
Dat <- Dat[order(Dat$ID_I, Dat$Date), ]

Dat[Dat$ID_I == 5174, ]
head(Dat,10)
# not looking too promising. just not enough people in cohorts
