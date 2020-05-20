###########################################################
####################### The Analysis ######################
###########################################################

alpha = 0.025
p = 0.05
M = 10000 # M simulations to calculate the critical values of the test statistics of the ES backtests by Acerbi and Szekely (2014,2017)
N = 10000 # N simulations to calculate the significances/power in the examples
set.seed(20180117)

###########################################################
# Significance | Will the ES Backtests accept a correct ES prediction? 
#

#
# A first impression
#

T = 250 # T = 500

# Choose a H0 assumption: 
H0 = c("normal",0,0.1,0,0)
H0 = c("stdt",0,0.1,5,0)
H0 = c("skewn",0,0.1,5,0)
H0 = c("skewn",0,0.1,-5,0)
H0 = c("skewt",0,0.1,5,0.5)
H0 = c("skewt",0,0.1,5,1.5)
H0 = c()
# H0 = c(...)
var_es_analytic(H0,alpha)

set.seed(2018)

# Calculate the significance 
crit_mc = crit_H0(T,H0,M,alpha,p)[[1]] 
significance(H0,T,"random",M,alpha,p,crit_mc) 
#


#
# Significance for different H0 assumptions with a fixed number of VaR exceedances
# 

T = 250 # T = 500
varreject = qbinom(1-p,T,0.01) - 1 # For more outliers, the VaR forecast is rejected by a simple binomial test
maxN = varreject + 10

significance_m = matrix(0,6,maxN)
rownames(significance_m) = c("A/S Test1", "A/S Test2", "A/S Test3", "A/S Test4", "C/P Test 1", "C/P Test 2")

# Choose a H0 assumption: 
H0 = c("normal",0,1,0,0)
H0 = c("stdt",0,1,5,0)
H0 = c("skewn",0,1,-2,0)
# H0 = c(...)

set.seed(2018)

# Calculation of the critical values for A/S tests 1-4 under H0:
# Note that if one has chosen a skewed distribution, sometimes there is an error which disappears if one loads the function crit_H0 again
crit_mc = crit_H0(T,H0,M,alpha,p)[[1]]

# Calculation of the significance with respect to the number of VaR exceedances:
for (i in 1:maxN){
  significance_m[,i] = significance(H0,T,i,M,alpha,p,crit_mc)
}

# Result: 
significance_m

# Plots
par(mfrow=c(2,3))
plot(seq(1,maxN,1),significance_m[1,], type="l", ylim = range(0:1), xlab="n", ylab="significance A/S 1")
plot(seq(1,maxN,1),significance_m[2,], type="l", ylim = range(0:1), xlab="n", ylab="significance A/S 2")
plot(seq(1,maxN,1),significance_m[3,], type="l", ylim = range(0:1), xlab="n", ylab="significance A/S 3")
plot(seq(1,maxN,1),significance_m[4,], type="l", ylim = range(0:1), xlab="n", ylab="significance A/S 4")
plot(seq(1,maxN,1),significance_m[5,], type="l", ylim = range(0:1), xlab="n", ylab="significance C/P 1")
plot(seq(1,maxN,1),significance_m[6,], type="l", ylim = range(0:1), xlab="n", ylab="significance C/P 2")
par(mfrow=c(1,1))



###########################################################
# Power | Will underestimated ES predictions be rejected? 
#

#
# A first impression
#

T=250 # 250 | 500

H0 = c("stdt",0,0.01,5,0)
var_es_analytic(H0,alpha)
H1 = c("stdt",0,0.01,3,0)
var_es_analytic(H1,alpha)
help = crit_H0(T,H0,M,alpha,p) 
crit_mc=help[[1]]
var_es_H0 = help[[2]]
power_calc(T,crit_mc,var_es_H0,H0,H1,N,alpha,p)

#
# Underestimated SD 
#

#
# H0: N(0,0.01) | T(0,0.01,df) | ...
# H1: The same distribution but SD(-scale) between (0.008 and 0.02)
#
T = 250 # 250 | 500

# Selection of null hypothesis
H0 = c("normal",0,0.01,0,0)
H0 = c("stdt",0,0.01,5,0)
H0 = c("skewt",0,0.01,5,1.5)
H0 = c("skewn",0,0.01,-5)
# H0 = c(...)

sd = seq(0.008,0.02,by=0.0005)
power = matrix(0,length(sd),6) 
mc_H0 = crit_H0(T,H0,M,alpha,p) 
crit_mc = mc_H0[[1]] # Critical values of A/S under H0
var_es_H0 = mc_H0[[2]]
for (i in 1:length(sd)){
  power[i,]=power_calc(T,crit_mc,var_es_H0,H0,c(H0[1],as.numeric(H0[2]),sd[i],as.numeric(H0[4]),as.numeric(H0[5])),N=1000,alpha,p)
}
# Plot: 
par(mfrow=c(1,1))
plot_power(sd,power,"sd","Title")


#
# H0: N(0,0.01)
# H1: T(0,0.01,df)
#
H0 = c("normal",0,0.01,0,0)
df = rev(seq(1,10,by=1))
power = matrix(0,length(df),6)
mc_H0 = crit_H0(T,H0,M,alpha,p)
crit_mc = mc_H0[[1]]
var_es_H0 = mc_H0[[2]]
for (i in 1:length(df)){
  power[i,]=power_calc(T,crit_mc,var_es_H0,H0,c("stdt",0,0.01,df[i],0),N=1000,alpha,p)
}
# Plot: 
plot_power(df,power,"df","T=250 | N(0,0.01) vs T(0,0.01,df)")


# 
# A fixed number of VaR exceedances
# 

H0 = c("normal",0,0.01,0,0)
H1 = c("normal",0,0.012,0,0)
var_es_analytic(H0,alpha)
mc_H0 = crit_H0(T,H0,M,alpha,p) 
crit_mc = mc_H0[[1]] 
var_es_H0 = mc_H0[[2]]
# Propability for an (H0)-VaR-exceedance under H1:
p_exc=pnorm(-var_es_H0[1,1],0,0.012)
# Expected VaR exceedances under H1:
floor(p_exc*T)
# Power in a setting with random exceedances:
power_calc(T,crit_mc,var_es_H0,H0,H1,N=1000,alpha,p)
power = matrix(0,25,6)
# Power in the case of i VaR exceedances:
for (i in 1:25){
  power[i,]=power_calc_fixn_H0H1(T,i,crit_mc,var_es_H0,H0,H1,N=10000,alpha,p)
}
power
pbinom(9,250,p_outlier)

# Plots
# 
par(mfrow=c(2,3))
plot(seq(1,25,1),power[,1], type="l", ylim = range(0:1), xlab="n", ylab="power A/S 1")
plot(seq(1,25,1),power[,2], type="l", ylim = range(0:1), xlab="n", ylab="power A/S 2")
plot(seq(1,25,1),power[,3], type="l", ylim = range(0:1), xlab="n", ylab="power A/S 3")
plot(seq(1,25,1),power[,4], type="l", ylim = range(0:1), xlab="n", ylab="power A/S 4")
plot(seq(1,25,1),power[,5], type="l", ylim = range(0:1), xlab="n", ylab="power C/P 1")
plot(seq(1,25,1),power[,6], type="l", ylim = range(0:1), xlab="n", ylab="power C/P 2")
par(mfrow=c(1,1))


###########################################################
# The influence of the VaR on the ES-backtests
#

# Consider the following 2 examples:
# 1) A correct ES prediction but a wrong VaR prediction 
# 2) An underestimated ES with a correct VaR prediction

############# Example 1: Correct ES, but wrong VaR prediction ##############
# F: The returns are iid distributed with c("stdt",0,0.01,3,0)
# H0: As described in the thesis


H1 = c("stdt",0,0.01,3,0);
T = 250 # 250 | 500
var_es_H1 = matrix(0,T,2)
es_H1 = var_es_analytic(H1,alpha)[[2]]
var_es_H1[,1] = rep(var_es_analytic(H1,alpha)[[1]],T)
var_es_H1[,2] = rep(var_es_analytic(H1,alpha)[[2]],T)

var_H0 = 0.02 # 0.02 | 0.025 | 0.03 | 0.0318 | 0.035 | 0.04 # (wrong) VaR prediction

# Simulation of an iid time series under H0:
example1 = function(T,alpha){
  x = numeric(T)
  q = runif(T,0,1)
  for (i in 1:T){
    if (q[i]<= alpha){
      x[i] = runif(1,-es_H1-(es_H1-var_H0),-var_H0)
    }
    if (q[i] > alpha){
      x[i] = runif(1,-var_H0,var_H0 + 2*alpha/(1-alpha)*es_H1) # The expected value is zero
    }
  }
  return(x)
}

LS = -es_H1-(es_H1-var_H0)
RS = var_H0 + 2*alpha/(1-alpha)*es_H1

#Check
x = example1(100000,alpha)
plot(ecdf(x)); plot(density(x))
mean(x)
-es_H1
mean(x[x<=-var_H0])
# The ES is correct, the VaR is wrong predicted. 

# VaR and ES unter H0: 
var_es_H0 = matrix(0,T,2)
var_es_H0[,1] = rep(var_H0,T)
var_es_H0[,2] = var_es_H1[,2] # The ES matches

# Prepare the demoninator in for the test statistic of A/S Test 3:
denom_t_H0 = function(){
  
  # Auxilliary function for integration in the denominator
  fb = function(p){
    return(pbinom(floor(T*alpha)-1,T-1,p))
  }
  
  # Quantile function for the distribution under H0
  quantile_function = function(q){
    v = 0
    if (q<alpha){
      v = LS + q/alpha*(-var_H0-LS) 
    }
    if (q>=alpha){
      v = -var_H0 + (q-alpha)/(1-alpha)*(RS-(-var_H0))
    }
    return(v)
  } 
  
  
  # Multiply
  summe = function(x){
    return(fb(x)*quantile_function(x))
  }
  
  denom = integrate(summe,0,1)
  denom = as.numeric(denom[1]) 
  denom = -T/(floor(T*alpha)) * denom
  
  return(denom) 
  
}
denom=rep(denom_t_H0(),T) # Ignore the warnings

# Critical values for A/S 1-4 under the null hypothesis
crit_mc = crit_H0_example1(T,var_es_H0,M,alpha,p,denom)

# Input for Corbetta/Peri Backtests: 
cp_pt = numeric(T)
for (i in 1:T){
  cp_pt[i] = alpha/2
}

significance_example1(var_es_H0,H1,T,M,alpha,p,crit_mc,denom)

############# Example 1b: Turning over the distributions #############
# H0: The returns are iid t-distributed with 3 df and SD-scale of 0.01
# H1: As described in the thesis

T = 250;  # 250 | 500
H0 = c("stdt",0,0.01,3,0)
es_H1 = var_es_analytic(H0,alpha)[2]
var_H1 = 0.0436 # 0.0436 | 0.0386 | 0.0336 | 0.0318 | 0.0286 | 0.0236 
var_es_H0 = matrix(0,T,2)
var_es_H0[,1] = rep(var_es_analytic(H0,alpha)[[1]],T)
var_es_H0[,2] = rep(var_es_analytic(H0,alpha)[[2]],T)
es_H0 = var_es_analytic(H0,alpha)[[2]]
var_H0 = var_es_analytic(H0,alpha)[[1]]
var_es_H1[,1] = rep(var_H1,T)
var_es_H1[,2] = var_es_H0[,2]

# Simulation of an iid time series under H1:
example1_vv = function(T,alpha){
  x = numeric(T)
  q = runif(T,0,1)
  for (i in 1:T){
    if (q[i]<= alpha){
      x[i] = runif(1,-es_H1-(es_H1-var_H1),-var_H1)
    }
    if (q[i] > alpha){
      x[i] = runif(1,-var_H1,var_H1 + 2*alpha/(1-alpha)*es_H1) 
    }
  }
  return(x)
}

# Check
x = example1(100000,alpha)
sort(x)[floor(100000*alpha)]
mean(x)
mean(x[x<=-var_H1])

# critical values under H0. Since one has a standard distribution, there is no need for an adjustment here. 
crit_mc = crit_H0(T,H0,M,alpha,p)[[1]]

# Calculation of the significance
significance_example1_vv(var_es_H0,H0,T,M,alpha,p,crit_mc)


############# Example 2: Wrong ES, but correct VaR prediction ############# 
T = 250 # 250 | 500
H0 = c("normal",0,0.01,0,0)
var_es_analytic(H0,alpha)
crit = crit_H0(T,H0,M,alpha,p)
crit_mc = crit[[1]]
var_es_H0 = crit[[2]]
df = 20
H1 = c("stdt",0,0.01,df,0)
var_es_analytic(H1,alpha)
power_calc(T,crit_mc,var_es_H0,H0,H1,N,alpha,p) # Power in that case where both the VaR and the ES are underestimated
# Shift, such that H0 and H1 have the same VaR
shift = var_es_analytic(H1,alpha)[1] - var_es_analytic(H0,alpha)[1]
H1_shift = c("stdt",shift,0.01,df,0)
var_es_analytic(H1_shift,alpha)
power_calc(T,crit_mc,var_es_H0,H0,H1_shift,N,alpha,p) # Power in that case where only the ES is underestimated


############# An example of deliberate deception ##############
# The true distribution of the returns is given by c("stdt",0,0.01,3,0), IID.
# The H0 assumptions is described in the thesis

H1 = c("stdt",0,0.01,3,0); 
T = 250 # T=500
var_es_analytic(H1,alpha)
var_es_H1 = matrix(0,T,2)  
var_H1 = var_es_analytic(H1,alpha)[[1]]
es_H1 = var_es_analytic(H1,alpha)[[2]]
var_es_H1[,1] = rep(var_es_analytic(H1,alpha)[[1]],T)
var_es_H1[,2] = rep(var_es_analytic(H1,alpha)[[2]],T)

var_H0 = 0.04 # deliberately overestimated VaR
es_H0 = 0.045 # deliberately underestimated ES

# iid time series under H0:
example1 = function(T,alpha){
  x = numeric(T)
  q = runif(T,0,1)
  for (i in 1:T){
    if (q[i]<= alpha){
      x[i] = runif(1,-es_H0-(es_H0-var_H0),-var_H0)
    }
    if (q[i] > alpha){
      x[i] = runif(1,-var_H0,var_H0 + 2*alpha/(1-alpha)*es_H0) 
    }
  }
  return(x)
}

#Check
x = example1(100000,alpha)
mean(x)
sort(x)[floor(100000*alpha)]
mean(x[x<=-var_H0])

# VaR und ES under H0: .
var_es_H0 = matrix(0,T,2)
var_es_H0[,1] = rep(var_H0,T) 
var_es_H0[,2] = rep(es_H0,T) 

# Preparation of the denominator for AS Test 3
denom_t_H0 = function(){
  
  # The same structure as above.
  
  fb = function(p){
    return(pbinom(floor(T*alpha)-1,T-1,p))
  }
  
  quantile_function = function(x){
    
    v = 0
    if (x<alpha){
      v = -es_H0-(es_H0-var_H0) + x/alpha*(-var_H0-(-es_H0-(es_H0-var_H0)))
    }
    if (x>=alpha){
      v = -var_H0 + (x-alpha)/(1-alpha)*((var_H0 + 2*alpha/(1-alpha)*es_H0)-(-var_H0))
    }
    return(v)
  }  
  
  
  summe = function(x){
    return(fb(x)*quantile_function(x))
  }
  
  denom = integrate(summe,0,1)
  denom = as.numeric(denom[1])
  denom = -T/(floor(T*alpha)) * denom
  
  return(denom) 
  
}
denom=rep(denom_t_H0(),T) # Ignore the warnings

# Critical values for the A/S Tests under H0
# Since the same structure as in example1 is used, one can use the adjusted function for the critical values here as well.
crit_mc = crit_H0_example1(T,var_es_H0,M,alpha,p,denom)

# Preparation of the vector for Corbetta and Peri
cp_pt = numeric(T)
for (i in 1:T){
  cp_pt[i] = alpha/2
}

# Calculation of the acceptance rate
# Since the structure is the same as in example1, one can use the adjusted function for "significance" 
significance_example1(var_es_H0,H1,T,M,alpha,p,crit_mc,denom)


###########################################################
# The power of tests in a setting with time-dependent volatility
#

############# An example of underestimated SD in a non-iid-setup ##############
#
#

# Calculation of the power in an IID Setup 
T = 250 # 250 | 500
H0 = c("normal",0,0.015,0,0)
H1 = c("normal",0,1.2*0.015,0,0)
help = crit_H0(T,H0,M,alpha,p)
crit_mc = help[[1]]
var_es_H0 = help[[2]]
power_calc(T,crit_mc,var_es_H0,H0,H1,N=M,alpha,p) # Power in the IID Setup

# Preparing of Volatility setting (denominator of AS3)
fb = function(p){
  return(pbinom(floor((T)*alpha)-1,(T)-1,p))  
}

# Vola setting: 
T = 250
x = numeric(T)
as_values = matrix(0,M,4)
cp_values = matrix(0,M,2)
crit_mc = numeric(4)
power = numeric(6)
denom=numeric(T)

# Simulation under H0 and H1: 
x_sim_H0= function(T,sigma){
  for (i in 1:T){
    H0 = c("normal",0,sigma[i],0,0)
    x[i] = pl_iid_sim(1,H0)
  }
  return(x)
}
x_sim_H1= function(T,sigma){
  for (i in 1:T){
    H1 = c("normal",0,sigma[i]*1.2,0)
    x[i] = pl_iid_sim(1,H1)
  }
  return(x)
}

# 1) Critical values
for (i in 1:M){
  sigma = runif(T,1,2)
  var_es_H0 = matrix(0,T,2)
  for (j in 1:T){
    var_es_H0[j,]=var_es_analytic(c("normal",0,sigma[j],0,0),alpha)
  }
  x = x_sim_H0(T,sigma)
  as_values[i,1] = as1(x,var_es_H0[,1],var_es_H0[,2])
  as_values[i,2] = as2(x,var_es_H0[,1],var_es_H0[,2],alpha)
  as_values[i,4] = as4(x,var_es_H0[,1],var_es_H0[,2],alpha)
  
  # AS 3
  for (j in 1:(T)){
    summe = function(x){
      return(fb(x)*qnorm(x,0,sigma[j]))
    }
    denom_help = integrate(summe,0,1)
    denom[j] = denom_help[[1]]
    denom[j] = -T/(floor(T*alpha)) * denom[j]
  }
  as_values[i,3]=as3_adjusted(x,alpha,denom)
  
}
crit_mc[1] = sort(as_values[,1])[floor(M*alpha)]
crit_mc[2] = sort(as_values[,2])[floor(M*alpha)]
crit_mc[3] = sort(as_values[,3])[floor(M*alpha)]
crit_mc[4] = sort(as_values[,4])[floor(M*alpha)]
crit_mc

# Preparing (y_t) for Corbetta and Peri:
# sd = 1
# H0 = c("normal",0,sd,0,0)
# pnorm(-var_es_analytic(H0,alpha)[2],0,as.numeric(H0[3]))
# Here, cp_pt is always the same for all SD!
cp_pt = rep(pnorm(-var_es_analytic(c("normal",0,1,0,0),alpha)[2],0,1),T)

# 2) Calculation of the Power:
for (i in 1:M){
  sigma = runif(T,1,2)
  var_es_H0 = matrix(0,T,2)
  for (j in 1:T){
    var_es_H0[j,]=var_es_analytic(c("normal",0,sigma[j],0,0),alpha)
  }
  x = x_sim_H1(T,sigma)
  # A/S 1,2: 
  as_values[i,1] = as1(x,var_es_H0[,1],var_es_H0[,2])
  as_values[i,2] = as2(x,var_es_H0[,1],var_es_H0[,2],alpha)
  # AS 3:
  for (j in 1:(T)){
    summe = function(x){
      return(fb(x)*qnorm(x,0,sigma[j]))
    }
    denom_help = integrate(summe,0,1)
    denom[j] = denom_help[[1]]
    denom[j] = -T/(floor(T*alpha)) * denom[j]
  }
  as_values[i,3]=as3_adjusted(x,alpha,denom)
  # A/S 4:
  as_values[i,4] = as4(x,var_es_H0[,1],var_es_H0[,2],alpha)
  
  # C/P: 
  cp_values[i,1]= cp_adjusted(x,var_es_H0[,2],cp_pt,p)[[3]]
  cp_values[i,2]= cp_adjusted(x,var_es_H0[,2],cp_pt,p)[[4]]
}
power[1] = length(as_values[,1][as_values[,1]<= crit_mc[1]])/M
power[2] = length(as_values[,2][as_values[,2]<= crit_mc[2]])/M
power[3] = length(as_values[,3][as_values[,3]<= crit_mc[3]])/M
power[4] = length(as_values[,4][as_values[,4]<= crit_mc[4]])/M
power[5] = mean(cp_values[,1])
power[6] = mean(cp_values[,2])
power # Power in the Non-Identical Setup
#

# 
##### The analysis similar to Du/Escanciano #####
#

par(mfrow=c(1,1))
T = 250 # 250 | 500
# H0 similar to Du/Escanciano
H0 = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),  mean.model = list(armaOrder = c(1,0), include.mean=FALSE), distribution.model ="norm", fixed.pars = list(ar1 = 0.05, omega = 0.05, alpha1 = 0.1, beta1 = 0.85))

# Auxiliary function for integration (A/S Test 3)
fb = function(p){
  return(pbinom(floor((T-1)*alpha)-1,(T-1)-1,p))
}

# Results: 

T = 250 # 500
# Note that this takes some time on a custom notebook, even for a calculation on multiple cores.
# To get a rough idea, one can take 
# M=1000; N=1000


# Results: 
cl <- makePSOCKcluster(detectCores() - 1); # Calculating on multiple cores 
registerDoParallel(cl, cores = detectCores() - 1); # Calculating on multiple cores 
crit_mc = garch_crit_mc(T,M,p=0.05)
garch_rejection_rate("H0",T,N,crit_mc,alpha,p)
garch_rejection_rate("A1",T,N,crit_mc,alpha,p)
garch_rejection_rate("A2",T,N,crit_mc,alpha,p)
garch_rejection_rate("A3",T,N,crit_mc,alpha,p)
garch_rejection_rate("A4",T,N,crit_mc,alpha,p)
garch_rejection_rate("A5",T,N,crit_mc,alpha,p)
stopCluster(cl) # Return to calculation on one core


