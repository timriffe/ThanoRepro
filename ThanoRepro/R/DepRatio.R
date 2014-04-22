# note, Robert's neat dependency ratio parlor trick:

# P65+ / P18-64 = alpha
# how does alpha change with r?

# change in alpha with r = (xbar18-64 - xbar65+) / 100
# (or something like that)


# example data, HMD, I forget which country, year, sex.... legit tho
Lx <- c(0.99273, 0.99194, 0.99147, 0.99114, 0.99089, 0.99068, 0.99048, 
        0.99031, 0.99015, 0.98999, 0.98982, 0.98963, 0.98942, 0.98918, 
        0.98888, 0.98846, 0.98786, 0.98705, 0.98601, 0.98478, 0.98347, 
        0.9821, 0.98074, 0.97941, 0.97809, 0.9768, 0.97556, 0.97431, 
        0.97304, 0.97176, 0.97047, 0.96916, 0.96779, 0.96634, 0.96476, 
        0.9631, 0.96135, 0.95946, 0.9574, 0.9552, 0.95283, 0.95027, 0.94753, 
        0.94453, 0.94127, 0.93777, 0.93401, 0.92996, 0.92553, 0.9208, 
        0.91584, 0.91063, 0.9051, 0.899, 0.89265, 0.88585, 0.87828, 0.87005, 
        0.86113, 0.85155, 0.8412, 0.83022, 0.8183, 0.80528, 0.79125, 
        0.77627, 0.76047, 0.74349, 0.72533, 0.70625, 0.68594, 0.66435, 
        0.64136, 0.61696, 0.59165, 0.56534, 0.53797, 0.50965, 0.4803, 
        0.44982, 0.41786, 0.38513, 0.35184, 0.31817, 0.2848, 0.2517, 
        0.21943, 0.18854, 0.1596, 0.13276, 0.10839, 0.08691, 0.06812, 
        0.05206, 0.03879, 0.02823, 0.02006, 0.01387, 0.00932, 0.00608, 
        0.00385, 0.00237, 0.00141, 0.00081, 0.00046, 0.00025, 0.00013, 
        7e-05, 3e-05, 2e-05, 1e-05)

Plotka <- function(Lx, r, a){
    sum(Lx * exp(-r*a))
}
xbarLotka <- function(Lx, r, a){
    sum(a * Lx * exp(-r*a)) /  sum(Lx * exp(-r*a))
}

age <- 0:110
i.65plus <- age >= 65
i.18to64 <- age >= 18 & age < 65
a.65plus <- age[i.65plus] + .5
a.18to64 <- age[i.18to64] + .5

# difference between mean ages
Diffa <- function(Lx, r){
    xbarLotka(Lx[i.18to64], r, a = a.18to64) -
            xbarLotka(Lx[i.65plus], r, a = a.65plus)
}
# the old-age dependency ratio
alphaLotka <- function(Lx,r){
    Plotka(Lx[i.65plus], r, a = a.65plus) /
            Plotka(Lx[i.18to64], r, a = a.18to64)
}

r <- 0
alphaLotka(Lx,r)
# change in old-age dependency from 1% increase in r
# the difference divided by the original value
(alphaLotka(Lx, .01) - alphaLotka(Lx, 0)) / alphaLotka(Lx, 0)
# change in old-age dependency from 1% decrease in r
(alphaLotka(Lx, -.01) - alphaLotka(Lx, 0)) / alphaLotka(Lx, 0) 
# or should I put some midpoint in the denominator?
# 1% increase: closer to mean age differences below, but not exact
(alphaLotka(Lx, .01) - alphaLotka(Lx, 0)) / alphaLotka(Lx, 0.005)

# could also use the log of the ratio:
log(alphaLotka(Lx, .01) / alphaLotka(Lx, 0))
Diffa(Lx,0.005) / 100 # differs in 6th digit! closer, but still not quite

# also rather close, but to cigar:
(Diffa(Lx,0) / 100 +
Diffa(Lx,0.01) / 100)/2
# exact same as r=0.01 / r=0:
  
# still not convinced that this is the proper setup. 
# none of these match, although we get an approximate equivalence,
# depending on how we space things out.
# difference in mean age of these age-classes as percent? r = 0
Diffa(Lx,0) / 100 
Diffa(Lx,0.005) / 100 
# r = 0.01
Diffa(Lx,.01) / 100
# r = -0.01
Diffa(Lx,-.01) / 100

# my question: I was expecting an exact relationship, but see that we've got some 
# discrete interference going on. I still expect there to be a way to set this up 
# so that things come out exact. Any tips?