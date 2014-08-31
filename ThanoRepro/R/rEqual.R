source("/home/tim/git/ThanoRepro/ThanoRepro/R/Functions.R")
#  thanatological age 

rThanoCoale <- compiler::cmpfun(function(Fy, da, .a = .5:110.5, maxit = 1e2, tol = 1e-15){  
            
            # Based on Coale (1957), modified.
            N    <- length(Fy)
            dxM  <- matrix(0, ncol = N, nrow = N)
            dxi  <- da
            for (i in 1:N){
                dxM[i, 1:length(dxi)  ] <- dxi 
                dxi <- dxi[2:length(dxi) ]
            }     
            R0      <- sum(dxM * Fy)
            T.guess <- sum(.a * dxM * Fy) / R0 # assuming r = 0
            r.i     <- log(R0) / T.guess
            # be careful to discount Fy by SRB appropriately for males / females
            # prior to specification
            for (i in 1:maxit){ # 15 is more than enough!
                #cat(r2,i,"\n")
                deltai <- 1 - sum(rowSums(t(t(dxM) * exp(-r.i * .a))) * Fy)
                # the mean generation time self-corrects 
                # according to the error produced by the Lotka equation
                r.i <- r.i - (deltai / (T.guess - (deltai / r.i))) 
                if (abs(deltai) <= tol){
                    break
                }
            }
            return(r.i)  
        })
# chronological age
Rmomentn <- compiler::cmpfun(function(fx,Lx,x,n=0){
            sum((x^n)*fx*Lx)
        })
rLotkaCoale <- compiler::cmpfun(function(fx,Lx,x=.5:110.5, maxit = 1e2, tol = 1e-15){
            # from Coale, Ansley J. (1957) A New Method for Calculating Lotka's r- the Intrinsic Rate of Growth in a Stable Population.
            # Population Studies, Vol. 11 no. 1, pp 92-94
            R0 <- Rmomentn(fx,Lx,x,0)
            # first assuming a mean generation time of 29
            ri <- log(R0)/29
            for (i in 1:maxit){ # 10 is more than enough!
                deltai <- sum(exp(-ri*x)*fx*Lx)-1
                # the mean generation time self-corrects 
                # according to the error produced by the Lotka equation
                ri <- ri+(deltai/(29-(deltai/ri)))
                if (abs(deltai) <= tol){
                    break
                }
            }
            return(ri)  
        })
# took me a while to get to this...
#install.packages("data.table")
library(data.table)
# read in data and select females
Data                        <- local(get(load("/home/tim/git/ThanoRepro/ThanoRepro/Data/DataAll.Rdata")))
Data                        <- Data[Data$Sex == "f", ]
Data$Fy[is.na(Data$Fy)]     <- 0
Data$Fyf[is.nan(Data$Fyf)]  <- 0
DATA                        <- data.table(Data)

# head(Data)
# f(a) is the fertility function, not a density function!
fa <- with(Data, Fxf[Year == 2000 & Code == "SWE"])
# lifetable exposure, not same as lx, but approximately. radix =1 
La <- with(Data, Lx[Year == 2000 & Code == "SWE"])
# d(a) is the lifetable deaths distribution, here it sums to 1.
da <- with(Data, dx[Year == 2000 & Code == "SWE"])
# this is mu(a), don't confuse it with the maternity function... we use 'f' for fertility here
ma <- with(Data, mx[Year == 2000 & Code == "SWE"])

# Cohen's identity:
# sum(da*1/ma)
# with(Data, ex[Year == 2000 & Code == "SWE"])[1]
#####

# stable age given survival (lifetable exposure),r, and age
getca <- function(La,r,a=.5:110.5){
    La * exp(a * -r) / sum(La * exp(a * -r))
}

# the standard approximations:
a   <- .5:110.5
(rL  <- rLotkaCoale(fa,La))
ca  <- getca(La,rL) # stable age distribution
ba  <- fa*ca        # stable age-specific births (hypothetical)

sum(fa*La*exp(-rL*a)) # should equal 1, the lotka equation.

# implied stable thanatological age structure (equal when r = 0)
.cy <- rowSums(Thano(ca,da))
# implied thanatological-age-specific birth distribution (hypothetical)
.by <- rowSums(Thano(ba,da))
# stable thanatological fertility rates.
.fy <- .by/.cy

# this is one kind of equality (tautology):
sum(.fy*.cy)
sum(fa*ca)
# and the same would hold of deaths...
# which means that b-d = r would be the same for both...

# another way would be to prove that the sum is equal if r is the 
# same and .fy is calculated from the stable pop
rL
rThanoCoale(.fy,da,.a=a) # 4th decimal place... is it a rounding error or not?
# hard to tell when r is close to zero anyway. We suspect that they are equal,
# and that this is due to rounding. But in order to show this, we need to ensure that
# inputs are perfectly conformable and do away with some approximations that entail
# rounding error. 

# 1) use l(a) instead of L(a)
# further make age strictly integer...
la2 <- with(Data, lx[Year == 2000 & Code == "SWE"])
# 2) force d(a) to entail l(a) and vice versa.
da2 <- -diff(c(la2,0)) # extra step to force perfect consistency
                      # between d(a) and l(a)

# 3) do away with mid-interval age. Use integer.
a2      <- 0:110
# 4) now recalculate the resulting r
(rL2    <- rLotkaCoale(fa,la2,x=a2))
rL # compare, close but not same.

# 5) recalc stable age and births using these 'clean' inputs
ca2     <- getca(la2,rL2,a2)
ba2     <- fa*ca2 

# 6) and convert these to thanatological equivalents *without*
# staggering the d(a) vector.
.cy2    <- rowSums(Thano(ca2,da2,stagger=FALSE))
.by2    <- rowSums(Thano(ba2,da2,stagger=FALSE))
# these is the 'clean' stable thanatological fertility schedule
.fy2    <- .by2 / .cy2

# now compare Chronos with Thanatos:
rL2
rThanoCoale(.fy2,da2,.a=a2)
# and there we have it! perfect identity! this means that the previous difference
# was due to rounding, in order to get a better approximation.

# easier to prove on paper if we know this is true...

# one thing to verify:
#f^\star(y) is the c(a)d(y|a)-weighted average of f(a)...

sum(.fy2)
sum(fa)

dya <- da2fya(da2,FALSE)
dxM  <- matrix(0, ncol = 111, nrow = 111)
dxi  <- da2
for (i in 1:111){
    dxM[i, 1:length(dxi)  ] <- dxi 
    dxi <- dxi[2:length(dxi) ]
} 
# these are identical, save for some machine error:
all(abs(dya * la2 - dxM) < 1e-12)
# demonstration:
all(colSums(fa * ca2 * dya) / colSums(ca2 * dya) == .fy2)

# plug into thano renewal eqn:
sum(
t(exp(-a2*rL2) * t(dya * la2)) * 
        colSums(fa * ca2 * dya) / colSums(ca2 * dya)
)
#wya <- ca2 * dya
## or define it as weights:
#sum(
#        t(exp(-a2*rL2) * t(dya * la2)) * 
#                colSums(fa * wya) / colSums(wya)
#)

# compare again w Lotka:
sum(exp(-a2*rL2)*la2*fa)

# note:
all(la2 == colSums(dya * la2))    # 1) d(y|a) * l(a) is just d(a) ...
all(la2 == rev(cumsum(rev(da2)))) # 2) d(a) and l(a) are commensurable here

# this is because d(y|a) is rescaled by dividing out l(a), 
# so we can just multiply it right back in


# this gives c(a):
all(la2 * exp(-a2 * rL2) / sum(la2 * exp(-a2 * rL2)) == ca2)

# plug it in:
sum(
        colSums(exp(-a2*rL2) * t(dya * la2)) * 
                colSums(fa * dya * la2 * exp(-a2 * rL2) / sum(la2 * exp(-a2 * rL2))) / 
                colSums(dya * la2 * exp(-a2 * rL2) / sum(la2 * exp(-a2 * rL2)))
)

# back to dxM, which is straight dx, unscaled:
sum(
        colSums(exp(-a2*rL2) * dxM) * 
                colSums(fa * dxM * exp(-a2 * rL2) / sum(la2 * exp(-a2 * rL2))) / 
                colSums(dxM * exp(-a2 * rL2) / sum(la2 * exp(-a2 * rL2)))
)
# let's define the exp() bit as growth:
growth <- exp(-a2*rL2) 
# this is still a thanatological setup:
sum(
        t(growth * dxM) * 
                colSums(fa * dxM * growth / sum(la2 * growth)) / 
                colSums(dxM * growth / sum(la2 * growth))
)

# back to Lotka, but increasing dimensions of l(a) to d(a) perspective
sum(growth * dxM * fa)

all(dxM == t(dxM))

# see about making comparable:
sum(
        t(growth * dxM) * 
                colSums(fa * dxM * growth / sum(dxM * growth)) / 
                colSums(dxM * growth / sum(dxM * growth))
)

# f^\star(y)
colSums(fa * dxM * growth / sum(dxM * growth)) / 
        colSums(dxM * growth / sum(dxM * growth))

sum(growth * dxM * fa)
sum(t(growth * dxM) * .fy2)

# and again:
sum(colSums(growth * dxM) * .fy2) # Thano
sum(rowSums(growth * dxM) * fa)   # Lotka

# and again:
sum(rowSums(growth * dxM) * fa)   # Lotka
                                  # Thano

sum(
        colSums(growth * dxM) * 
                colSums(fa * dxM * growth / sum(dxM * growth)) / 
                colSums(dxM * growth / sum(dxM * growth))
)
# cancel out divisors
sum(
        colSums(growth * dxM) * 
                colSums(fa * dxM * growth) / 
                colSums(dxM * growth)
)
# cancel out (multiply and divide
sum(colSums(fa * dxM * growth))
# bam!
sum(fa * dxM * growth)
# bam!
sum(fa * la2 * growth)
# this should serve to figure it out on paper
