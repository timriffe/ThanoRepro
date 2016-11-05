 DEPRECATED. JUST SOME FUN TESTING HERE


# data from http://www.demog.berkeley.edu/croatia/datadir/Linkfiles/linkfiles.htm
setwd(paste0("/home/",system("whoami"),"/git/ThanoEmpirical/ThanoEmpirical"))
#
CRO <- read.table("Data/croatdata4.txt.gz",
        stringsAsFactors = FALSE, sep = "\t", col.names = c('bid','bdate','sex','motherid','fatherid','mid1','mid2','mid3','mid4','mid5',
                'sid1','sid2','sid3','sid4','sid5','mdate1','mdate2','mdate3','mdate4','mdate5','did',
                'ddate','idk1','idk2','idk3','idk4','idk5','idk6','idk7','idk8','idk9','idk10','idk11',
                'idk12','idk13','idk14','dobk1','dobk2','dobk3','dobk4','dobk5','dobk6','dobk7','dobk8','dobk9',
                'dobk10','dobk11','dobk12','dobk13','dobk14','sidk1','sidk2','sidk3','sidk4','sidk5','sidk6','sidk7',
                'sidk8','sidk9','sidk10','sidk11','sidk12','sidk13','sidk14','remark1','remark2','remark3','remark4','remark5',
                'remark6','remark7','remark8','remark9','remark10','remark11','remark12','remark13',
                'remark14','remarok1','remarok2','remarok3','remarok4','remarok5','remarok6','remarok7',
                'remarok8','remarok9','remarok10','remarok11','remarok12','remarok13','remarok14','park1',
                'park2','park3','park4','park5','park6','park7','park8','park9','park10','park11','park12','park13','park14','lastdate')
)
# add origin to decimal dates:
DateIDS       <- which(grepl("date", colnames(CRO)) | grepl("dob", colnames(CRO)))

CRO[DateIDS]  <- lapply(CRO[,DateIDS], function(x){
            as.numeric(x) + 1700
        })

# only known dates of death
CRO <- CRO[!is.na(CRO$ddate), ]

dobkcols <- paste0("dobk",1:14)
# tricky. Make ever-had a child column, anyk:
CRO$anyk <- apply(CRO,1,function(x){
            any(!is.na(unlist(x[dobkcols])))
        })
# lifespan
CRO$LS <- CRO$ddate - CRO$bdate
# take just necessary columns for thano fertility
# put in long format
CROL <- do.call(rbind,lapply(1:14, function(cc, CRO){
                    
            DAT              <- CRO[,c("bdate","ddate", "sex", "anyk","LS", paste0("dobk", cc))]
            colnames(DAT)[5] <- "dobk"
            DAT$bo           <- cc
            DAT
        }, CRO = CRO))


# remove NAs of higher-than 1 order births if
# dim(CROL)
CROL     <- CROL[!with(CROL, bo > 1 & is.na(dobk)),]

# hist(CROL$bo[CROL$sex == "m" & !is.na(CROL$dobk)])

CROL$tak <- CROL$ddate - CROL$dobk
CROL$cak <- CROL$dobk - CROL$bdate
# uh oh. maybe linkage in wrong direction somehow?
#hist(CROL$tak[CROL$tak < 0 & !is.na(CROL$tak)])

# should remove cases more than -2:
CROL <- CROL[(!CROL$tak < -2) | is.na(CROL$dobk), ]

# now any remaining negatives get set to zero:
CROL$tak[CROL$tak < 0] <- 0

# round tak
CROL$takfloor <- floor(CROL$tak)
CROL$cakfloor <- floor(CROL$cak)

#matplot(t(table(CROL$sex,CROL$takfloor)),type='l',lty=1)

# make exposure table, based on LS
# when calculating cohort exposures by aggregating element-wise over
# lifelines, chronological and thanatological exposures are 
# IDENTICAL
# wow.
# note: exposures don't come from the long file, otherwise they get
# included redundantly for each birth...
ME <- rowSums(sapply(CRO$LS[CRO$sex=="m"], function(x){
                    c(rep(1,floor(x)),x-floor(x),rep(0,101-floor(x)-1))    
                }))
FE <- rowSums(sapply(CRO$LS[CRO$sex=="f"], function(x){
                    c(rep(1,floor(x)),x-floor(x),rep(0,101-floor(x)-1))    
                }))
E <- cbind(Female = FE, male = ME)
rownames(E) <- 0:100
#library(Pyramid)
#Pyramid(ME,FE)
TBirths <- table(CROL$takfloor,CROL$sex)
CBirths <- table(CROL$cakfloor,CROL$sex)
CB <- TB      <- E * NA
TB[rownames(TBirths),] <- TBirths
CB[rownames(CBirths),] <- CBirths

yticks <- pretty(CB/E)
png("Figures/CROfert.png")
par(xpd=TRUE,xaxs="i",yaxs="i",mai=c(1,.4,.4,.4))
plot(NULL,type="n",xlim=c(0,80),ylim=range(yticks),axes=FALSE,
        xlab="",ylab="",panel.first=list(
                rect(0,0,80,max(yticks),col=gray(.95),border=NA),
                segments(0,yticks,80,yticks,col="white"),
                segments(seq(10,70,by=10),0,seq(10,70,by=10),max(yticks),col="white"),
                text(0,yticks,yticks,pos=2,cex=.8),
                text(seq(0,80,by=20),0,seq(0,80,by=20),pos=1,cex=.8)))
matplot(0:100,CB/E,type = 'l',add=TRUE,lty=1,col=c("red","blue"))
matplot(0:100,TB/E,type = 'l',add=TRUE,lty=2,col=c("red","blue"))
arrows(40,-.025,60,-.025,length=.1)
arrows(40,-.018,20,-.018,length=.1,lty=2)
text(40,-.018,"years left",pos=4)
text(40,-.025,"years lived",pos=2)
legend("topright",bty="n",col=c("red","red","blue","blue"),lty=c(1,2,1,2),legend=c("females chronological","females thanatological", "males chronological", "males thanatological"))
text(2,max(yticks)*1.02,"Fertility Rate",pos=3)
dev.off()



