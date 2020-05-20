#######################################################
######################## Packages #####################
#######################################################
# Please install and load the following packages:
library(skewt)
library(sn)
library(fBasics)
library(poibin) 
library(rugarch)
library(Quandl) 
library(quantmod)
library(fGarch)
library(foreach) # install.packages("doParallel")
library(iterators) # install.packages("doParallel")
library(parallel) # install.packages("doParallel")
library(doParallel) # install.packages("doParallel")
set.seed(2018)
options(scipen = 30)
rm(list=ls()) # Remove all objects
graphics.off()

#######################################################
######################## Backtests ####################
#######################################################

# Acerbi/Szekely - Test 1 (2014)
as1 = function(x,var,es){
  
  # Inputs:
  # x   | time series
  # var | VaR prediction
  # es  | ES prediction
  
  # Output: 
  # z   | test statistic
  
  z = 0
  n = length(x)
  nt= 0
  
  for (i in 1:n)
  {
    if (x[i] < - var[i]) 
    {
      nt = nt + 1
      z = z + x[i]/es[i]
    }
  }
  z = z/nt + 1
  
  return(z)
}

# Acerbi/Szekely - Test 2 (2014)
as2 = function(x,var,es,alpha){

  # Inputs:
  # x     | time series
  # var   | VaR prediction
  # es    | ES prediction
  # alpha | confidence level
  
  # Output: 
  # z     | test statistic
  
  n = length(x)
  z = 0
  for (i in 1:n){
    
    if (x[i] < - var[i]) 
    {
      z = z + (x[i])/(n*alpha*es[i])
    }
  }
  z = z+1
  return(z)
}

# Acerbi/Szekely - Test3 (2014) (IID Setup)
as3 = function(x,H0,alpha){
  
  # Inputs:
  # x     | iid time series
  # H0    | assumption 
  # alpha | confidence level
  
  # Output: 
  # z     | test statistic
  
  T = length(x)
  x = sort(x)
  h = floor(T*alpha)
  H0_dist = H0[1]
  H0_par1 = as.numeric(H0[2]) # mu
  H0_par2 = as.numeric(H0[3]) # sigma/sd 
  H0_par3 = as.numeric(H0[4]) # df / skewness
  H0_par4 = as.numeric(H0[5]) # skewness for skewt
  numerator = 0
  denominator = 0
  
  # auxiliary function for integration in the denominator
  fb = function(p){
    return(pbinom(floor(T*alpha)-1,T-1,p))
  }
  
  if (H0_dist=="normal"){
    summe = function(x){
      return(fb(x)*qnorm(x,H0_par1,H0_par2)) # 18.12. Hier geht es so
    }
  }
  if (H0_dist=="stdt"){
    summe = function(x){
      # 18.12.: return(fb(x)*H0_par2*qt((x-H0_par1),H0_par3)) # 12.12. vorher: return(fb(x)*qt((x-H0_par1)/H0_par2,H0_par3)) 
      return(fb(x)*((H0_par2*qt((x),H0_par3))+H0_par1)) # neu am 18.12.
      }
  }
  if (H0_dist=="skewn"){
    summe = function(x){
      return(fb(x)*(H0_par2*qsn((x),alpha=H0_par3)+H0_par1)) # 18.12.
    }
  }
  if (H0_dist=="skewt"){
    summe = function(x){
      return(fb(x)*(H0_par2*qskt((x),H0_par3,H0_par4)+H0_par1)) # 18.12.
    }
  }
  
  # Calculate the numerator
  for (i in 1:h){
    numerator = numerator + x[i]
  }
  numerator = -1/h * numerator
  
  # Calculate the denominator
  denominator = integrate(summe,0,1)
  denominator = as.numeric(denominator[1])
  denominator = -T/(floor(T*alpha)) * denominator
  
  z = - (numerator/denominator) + 1 
  
  return(z)
  
}

# Acerbi/Szekely - Test 3 (2014) adjustet with individual denominator
as3_adjusted = function(x,alpha,denominator){
  
  # Inputs:
  # x           | iid time series
  # alpha       | confidence level
  # denominator | the denominator as a vector
  
  # Output: 
  # z           | test statistic
  
  T = length(x)
  x = sort(x)
  h = floor(T*alpha)
  numerator = 0
  
  # Calculate numerator: 
  for (i in 1:h){
    numerator = numerator + x[i] 
  }
  
  numerator = T*numerator 
  numerator = -1/h * numerator
  
  denominator_sum = sum(denominator)
  
  z = - (numerator/denominator_sum) + 1 
  
  return(z)
  
}

# Acerbi/Szekely - Test 4 (2017) (IID Setup)
as4 = function(x,var,es,alpha){
  
  # Inputs:
  # x     | time series
  # var   | VaR prediction
  # es    | ES prediction
  # alpha | confidence level
  
  # Output: 
  # z     | test statistic
  
  T = length(x)
  z = 0
  
  for (i in 1:T){
    add = alpha*(es[i]-var[i])
    if (x[i] < - var[i]){
      add = add + x[i] + var[i]
    }
    add = add/(es[i]*alpha)
    z = z + add
  }
  z = z/T
  return(z)
}

# The simple BCBS VaR binomial test (99% VaR prediction)
var_test = function(x,var,p){
  
  # Inputs:
  # x         | time series
  # var       | VaR prediction, should be the 99 % prediction
  # p         | legal probability of error
  
  # Output: 
  # reject    | BOOLEAN
  
  T = length(x)
  reject = FALSE
  
  # Count Exceedances
  n = 0
  for (i in 1:T){
    if (x[i] <= -var[i]){
      n = n+1
    }
  }
  
  # If the p_value (1 - pbinom(n,t,0.01)) is lower than the probability of error, the tests rejects
  # Note that 0.01 is selected analogous to Basel requirements
  
  if (1 - pbinom(n,T,0.01) < p){
    reject = TRUE
  }
  
  return (reject)
}

# Corbetta/Peri - (2016) (IID Setting)
cp = function(x,rho,H0,p){
  
  # Inputs:
  # x         | time series
  # rho       | risk measure prediction
  # H0        | assumption 
  # p         | significance level
  
  # Output: 
  # z1        | test statistic of test 1
  # z2        | test statistic of test 2
  # reject1   | rejection of test 1 (as BOOLEAN)
  # reject2   | rejection of test 2 (as BOOLEAN)
  
  T = length(x)
  pb = numeric(T)
  lambda = 0
  z1 = 0
  z2 = 0
  h1 = 0 # numerator
  h2 = 0 # denominator
  H0_dist = H0[1]
  H0_par1 = as.numeric(H0[2]) # mu
  H0_par2 = as.numeric(H0[3]) # sigma/sd 
  H0_par3 = as.numeric(H0[4]) # df / skewness
  H0_par4 = as.numeric(H0[5]) # skewness for skewt
  reject1 = FALSE
  reject2 = FALSE
  
  for (i in 1:T){
    y = -rho[i] 
    
    if (H0_dist=="normal"){
      lambda = pnorm(y,H0_par1,H0_par2)
    }
    if (H0_dist=="stdt"){
      lambda = pt((y-H0_par1)/H0_par2,H0_par3) 
    }
    if (H0_dist=="skewn"){
      lambda = psn(x=(y-H0_par1)/H0_par2,alpha=H0_par3)
    }
    if (H0_dist=="skewt"){
      lambda = pskt((y-H0_par1)/H0_par2,H0_par3,H0_par4) 
    }
    pb[i] = lambda
    if (x[i]<y){
      z1 = z1 + 1
      h1 = h1 + (1 - lambda)
    }
    if (x[i]>y){
      h1 = h1 - lambda
    }
    h2 = h2 + lambda*(1-lambda)
  }
  h2 = sqrt(h2)
  z2 = h1/h2
  
  # If the p value of Z1 ~ Poiss.Bin(lambda) is lower than the significance level, the test rejects
  if((1-ppoibin(z1,pp=pb))<p){
    reject1 = TRUE
  }
  
  # If the p value of Z2 ~ N(0,1) is lower than the significance level, the test reject
  #
  # Since one does not know exactly in which direction the test statistic 
  # develops under H1, one must test in both directions.
  # 
  
  if((1-pnorm(z2))<p/2){
    reject2 = TRUE
  }
  if(pnorm(z2)<p/2){
    reject2 = TRUE
  }
  
  return(list(z1,z2,reject1,reject2))
}

# Corbetta/Peri functions, such that lambda_t = P_t(-rho_t) is delivered seperately. Non-IID and non-standard distributions possible.
cp_adjusted = function(x,rho,cp_pt,p){
  
  # Inputs:
  # x         | time series
  # rho       | risk measure prediction
  # cp_pt     | lambda_t as a vector 
  # p         | significance level
  
  # Output: 
  # z1        | test statistic of test 1
  # z2        | test statistic of test 2
  # reject1   | rejection of test 1 (as BOOLEAN)
  # reject2   | rejection of test 2 (as BOOLEAN)
  
  T = length(x)
  pb = numeric(T)
  lambda = 0 
  z1 = 0
  z2 = 0
  h1 = 0 
  h2 = 0
  reject1 = FALSE
  reject2 = FALSE
  
  for (i in 1:T){
    y = -rho[i] 
    lambda = cp_pt[i]
    pb[i] = lambda
    
    if (x[i]<y){
      z1 = z1 + 1
      h1 = h1 + (1 - lambda)
    }
    if (x[i]>y){
      h1 = h1 - lambda
    }
    h2 = h2 + lambda*(1-lambda)
  }
  h2 = sqrt(h2)
  z2 = h1/h2
  
  # If the p value of Z1 ~ Poiss.Bin(lambda) is lower than the significance level, the test rejects:
  if((1-ppoibin(z1,pp=pb))<p){
    reject1 = TRUE
  }
  
  # If the p value of Z2 ~ N(0,1) is lower than the significance level, the test rejects:
  if((1-pnorm(z2))<p/2){
    reject2 = TRUE
  }
  if(pnorm(z2)<p/2){
    reject2 = TRUE
  }
  
  return(list(z1,z2,reject1,reject2))
}

###########################################################
######################## Functions ########################
###########################################################

# Simulation of an iid-time-series
pl_iid_sim = function(T,H0){
  
  # Inputs:
  # T     | length of time series
  # H0    | assumption
  
  # H0 = c(dist,parameter1,parameter2,parameter3,parameter4)
  # If dist = normal: parameter1 = mu, parameter2 = sd
  # If dist = st: parameter1 = mu, parameter2 = sigma, parameter3 = df
  # If dist = skewn: parameter1 = mu, parameter2 = sigma, parameter3 = alpha
  # If dist = skewt: parameter1 = mu, parameter2 = sigma, parameter3 = df, parameter4 = gamma
  # Note that sigma is the SD-Shift for all distributions exept the N(mu,sd) one
  
  # Output: 
  # x     | simulated time series
  
  x = numeric(T)
  dist = H0[1]
  parameter1 = as.numeric(H0[2]) 
  parameter2 = as.numeric(H0[3])
  parameter3 = as.numeric(H0[4])
  parameter4 = as.numeric(H0[5])
  
  if (dist=="normal"){
    mu = parameter1
    sd= parameter2
    x = sd*rnorm(T) + mu
  }
  if (dist=="stdt"){
    mu = parameter1
    sigma = parameter2
    df = parameter3
    x = mu + sigma*rt(T,df) 
  }
  if (dist=="skewn"){
    mu = parameter1
    sigma = parameter2
    shape = parameter3
    x = as.numeric(rsn(n=T, xi=mu -sigma*shape/sqrt(1+(shape)^2)*sqrt(2/pi), omega=sigma, alpha=shape))
    # Since with the shape, die expected value is moved by (shape/sqrt(1+(shape)^2)*sqrt(2/pi)) 

  }
  if (dist=="skewt"){
    mu = parameter1
    sigma = parameter2
    df = parameter3
    gamma = parameter4
    x = sigma*rskt(T,df,gamma)
    x = x - mean(x) # Mean adjustment due to the skewness
    x = x + mu
  }
  return(x)
}

# Simulation of an iid-time-series with n VaR exceedances
pl_iid_sim_fixn = function(T,H0,n,alpha){
  
  # Inputs:
  # T     | length of time series
  # H0    | assumption
  # n     | number of VaR exceedances
  # alpha | confidence level
  
  # H0 = c(dist,parameter1,parameter2,parameter3,parameter4)
  # If dist = normal: parameter1 = mu, parameter2 = sd
  # If dist = st: parameter1 = mu, parameter2 = sigma, parameter3 = df
  # If dist = skewn: parameter1 = mu, parameter2 = sigma, parameter3 = alpha
  # If dist = skewt: parameter1 = mu, parameter2 = sigma, parameter3 = df, parameter4 = gamma
  # Note that sigma is the SD-Shift for all distributions exept the N(mu,sd) one
  
  # Output: 
  # x     | simulated time series
  
  x = numeric(T)
  outlier = numeric(n)
  rest = numeric(T-n)
  dist = H0[1]
  par1 = as.numeric(H0[2])
  mu = par1
  par2 = as.numeric(H0[3])
  par3 = as.numeric(H0[4])
  par4 = as.numeric(H0[5])
  
  q1 = runif(n,min=0,max=alpha)
  q2 = runif(T-n,min=alpha,max=1)
  if (dist=="normal"){
    outlier = qnorm(q1,par1,par2)
    rest = qnorm(q2,par1,par2)
  }
  if (dist=="stdt"){
    
    aux = pt(-var_es_analytic(H0,alpha)[1]/par2,par3) 
    q1 = runif(n,min=0,max=aux)
    q2 = runif(T-n,min=aux,max=1)
    
    outlier = par2*qt(q1,par3)
    rest = par2*qt(q2,par3)
  } 
  if (dist=="skewn"){  
    shape = par3
    mu= par1 -shape/sqrt(1+(shape)^2)*sqrt(2/pi) 
    outlier = qsn(q1,xi=mu,omega=par2,alpha=shape)
    rest = qsn(q2,xi=mu,omega=par2,alpha=shape)
  }
  if (dist=="skewt"){
    outlier = par2*qskt(q1,par3,par4)
    rest = par2*qskt(q2,par3,par4)
    
    # Mean Adjustment due to the skewness
    sample = par2*rskt(T,df=par3,gamma=par4)
    mean_adjust = mean(sample)
    
    outlier = outlier + mu - mean_adjust
    rest = rest + mu - mean_adjust
  }
  x = c(outlier,rest)
  return(x)
}

# Calculation of VaR und ES for standard distributions
var_es_analytic = function(H0, alpha){
  
  # Inputs:
  # H0    | assumption
  # alpha | confidence level
  
  # H0 = c(dist,parameter1,parameter2,parameter3,parameter4)
  # If dist = normal: parameter1 = mu, parameter2 = sigma
  # If dist = st: parameter1 = mu, parameter2 = sigma, parameter3 = df
  # If dist = skewn: parameter1 = mu, parameter2 = sigma, parameter3 = alpha
  # If dist = skewt: parameter1 = mu, parameter2 = sigma, parameter3 = df, parameter4 = gamma
  # Note that sigma is the SD-Shift for all distributions exept the N(mu,sd) one
  
  # Output: 
  # var   | Value-at-Risk
  # es    | Expected Shortfall 
  
  dist = H0[1]
  parameter1 = as.numeric(H0[2])
  parameter2 = as.numeric(H0[3])
  parameter3 = as.numeric(H0[4])
  parameter4 = as.numeric(H0[5])
  
  if (dist=="normal"){ 
    mu = parameter1
    sd = parameter2
    var = -sd*qnorm(alpha) - mu
    es = sd*dnorm(qnorm(1-alpha))/alpha - mu
  }
  
  if(dist == "stdt"){ 
    mu = parameter1 
    sd = parameter2
    df = parameter3
    var = -sd*qt(alpha,df) - mu
    b = beta(df/2,0.5)
    es = 2*sd/(alpha*2*sqrt(df)*b)*1/((1+((qt(alpha,df))^2/df))^((df+1)/2))*(df + (qt(alpha,df))^2)/(df-1) - mu
  }
  
  if(dist == "skewn"){
    
    sd = parameter2
    al = parameter3
    shape = parameter3
    mu = parameter1 - sd*shape/sqrt(1+(shape)^2)*sqrt(2/pi) 
    var = - sd*qsn(alpha,xi=0,omega=1, alpha=al, tau=0, dp=NULL, tol=1e-10, solver="NR") - mu
    yp= (-var-mu)/sd
    zp=sqrt(1+(al)^2)*yp 
    # Bernadi (2012)
    # es = - ((sd*sqrt(2))/(alpha*sqrt(pi))*(al*pnorm(zp) - sqrt(2*pi)*dnorm(yp)*pnorm(al*yp))) - mu
    es = - sd*as.numeric((integrate((function(x) qsn(x,omega=1,alpha=al)),0,alpha))[1])/alpha - mu 
  }
  
  if(dist == "skewt"){
    mu = parameter1
    sd = parameter2
    df = parameter3
    gamma = parameter4
    help = sd*rskt(100000,df,gamma) # Mean Adjustment
    mean_adjust = mean(help) # Mean Adjustment
    var = -sd*qskt(alpha,df,gamma) - mu + mean_adjust
    es = - sd*as.numeric((integrate((function(x) qskt(x,df,gamma)),0,alpha))[1])/alpha - mu + mean_adjust
  }
  
  return(c(var,es))
}

# Critical values for the distributions of the test statistics of A/S 1-4 for standard distributions (IID Setup)
crit_H0 = function(T,H0,M,alpha,p){
  
  # Inputs:
  # T         | length of time series
  # H0        | assumption
  # M         | number of Monte Carlo simulations
  # alpha     | significance level
  # p         | significance level
  
  # Outputs:
  # crit_mc   | vector with 4 critical values for A/S tests 1-4
  # var_es_H0 | vectors of VaR and ES for each day (iid assumption)

x = numeric(T)
z_H0_mc = matrix(0,M,4) 
crit_mc = numeric(4)
var_es_H0_help = c(0,0)
var_es_H0 = matrix(0,T,2) 

# Calculate VaR and ES under H0 and save it as a vector for every day t (iid assumption)
var_es_H0_help = var_es_analytic(H0,alpha)
var_es_H0[,1] = rep(var_es_H0_help[1],T)
var_es_H0[,2] = rep(var_es_H0_help[2],T)

# Adjustment for the skewed distributions
# The 3rd Test of A/S and the Tests of C/P calculate with values using the quantile function. 
# The re-adjustment of the mean as described in the master thesis is already transfered here. 

if(H0[[1]]=="skewn")
{
  mu = as.numeric(H0[2])
  sigma = as.numeric(H0[3])
  shape = as.numeric(H0[4])
  mean_adjust = sigma*shape/sqrt(1+(shape)^2)*sqrt(2/pi)
  H0_adj = c("skewn",mu - mean_adjust,sigma,shape,0)
}

if(H0[[1]]=="skewt")
{
  mu = as.numeric(H0[2])
  sigma = as.numeric(H0[3])
  df = as.numeric(H0[4])
  gamma = as.numeric(H0[5])
  mean_adjust = mean(sigma*rskt(100000,df,gamma))
  H0_adj = c("skewt",mu-mean_adjust,sigma,df,gamma)
}
################ 


# Calculate a Monte-Carlo sample of Z_i under H0:
for (i in 1:M){
  x = pl_iid_sim(T,H0)
  z_H0_mc[i,1]=as1(x,var_es_H0[,1],var_es_H0[,2])
  z_H0_mc[i,2]=as2(x,var_es_H0[,1],var_es_H0[,2],alpha)
  if (H0[[1]]=="skewn" || H0[[1]]=="skewt"){
    z_H0_mc[i,3]=as3(x,H0_adj,alpha) # H0 Adjusted
  }
  if (H0[[1]]=="normal" || H0[[1]]=="stdt"){
    z_H0_mc[i,3]=as3(x,H0,alpha) 
  }
  #z_H0_mc[i,3]=as3(x,H0,alpha)
  z_H0_mc[i,4]=as4(x,var_es_H0[,1],var_es_H0[,2],alpha)
}

# Calculate the critical values: 
z_H0_mc[,1]=sort(z_H0_mc[,1],na.last=TRUE) 
z_H0_mc[,2]=sort(z_H0_mc[,2])
z_H0_mc[,3]=sort(z_H0_mc[,3])
z_H0_mc[,4]=sort(z_H0_mc[,4])
crit_mc[1]=z_H0_mc[,1][floor(M*p)]
crit_mc[2]=z_H0_mc[,2][floor(M*p)]  
crit_mc[3]=z_H0_mc[,3][floor(M*p)]
crit_mc[4]=z_H0_mc[,4][floor(M*p)]

return(list(crit_mc,var_es_H0))

}

# Calculation of the power of the tests for different standard H0 and H1 distributions (IID Setup)
power_calc = function(T,crit_mc,var_es_H0,H0,H1,N,alpha,p){
  
  # Inputs:
  # T         | length of time series
  # crit_mc   | critical values of the tests 1-4 under H0
  # var_es_H0 | vector with VaR and ES for each day (iid assumption)
  # H0        | H0 assumption
  # H1        | actual assumption
  # N         | number of Monte Carlo simulations which are done to calculate the power
  # alpha     | confidence level
  # p         | significance level
  
  # Outputs:
  # power     | vector with 6 values, power of tests A/S 1-4, C/P 1,2
  
  x = numeric(T)
  z_H1_mc = matrix(0,N,4) 
  cp_tests = matrix(0,N,4) 
  var_tests = numeric(T)
  power = numeric(6)
  
  # Adjustment for the skewed distributions
  # The 3rd Test of A/S and the Tests of C/P calculate with values using the quantile function. 
  # The re-adjustment of the mean as described in the master thesis is already transfered here. 
  
  if(H0[[1]]=="skewn")
  {
    mu = as.numeric(H0[2])
    sigma = as.numeric(H0[3])
    shape = as.numeric(H0[4])
    mean_adjust = sigma*shape/sqrt(1+(shape)^2)*sqrt(2/pi)
    H0_adj = c("skewn",mu - mean_adjust,sigma,shape,0)
  }
  
  if(H0[[1]]=="skewt")
  {
    mu = as.numeric(H0[2])
    sigma = as.numeric(H0[3])
    df = as.numeric(H0[4])
    gamma = as.numeric(H0[5])
    mean_adjust = mean(sigma*rskt(100000,df,gamma))
    H0_adj = c("skewt",mu-mean_adjust,sigma,df,gamma)
  }
  ################ 
  
  for (i in 1:N){
    
    # Calculate a Monte-Carlo sample of Z under H1 (with assumed VaR and ES from H0):
    x = pl_iid_sim(T,H1)
    z_H1_mc[i,1]=as1(x,var_es_H0[,1],var_es_H0[,2])
    z_H1_mc[i,2]=as2(x,var_es_H0[,1],var_es_H0[,2],alpha)
    if (H0[[1]]=="skewn" || H0[[1]]=="skewt"){
      z_H1_mc[i,3] = as3(x,H0_adj,alpha) # H0 Adjusted!   
    }
    if (H0[[1]]=="normal" || H0[[1]]=="stdt"){
      z_H1_mc[i,3] = as3(x,H0,alpha)    
    }
    z_H1_mc[i,4]=as4(x,var_es_H0[,1],var_es_H0[,2],alpha)
    
    # Test C/P directly:
    if (H0[[1]]=="skewn" || H0[[1]]=="skewt"){
      cp_test = cp(x,var_es_H0[,2],H0_adj,p) # H0 Adjusted!   
    }
    if (H0[[1]]=="normal" || H0[[1]]=="stdt"){
      cp_test = cp(x,var_es_H0[,2],H0,p)   
    }

    cp_tests[i,1]=cp_test[[3]]
    cp_tests[i,2]=cp_test[[4]]

    
  }
  
  # Calculate the power of the tests: 
  power[1]=length(z_H1_mc[,1][z_H1_mc[,1]<crit_mc[1]])/N
  power[2]=length(z_H1_mc[,2][z_H1_mc[,2]<crit_mc[2]])/N
  power[3]=length(z_H1_mc[,3][z_H1_mc[,3]<crit_mc[3]])/N
  power[4]=length(z_H1_mc[,4][z_H1_mc[,4]<crit_mc[4]])/N
  power[5]=mean(cp_tests[,1])
  power[6]=mean(cp_tests[,2])
  
  return(power)
  
}

# Calculation of the significance of the tests (IID Setup) 
significance = function(H0,T,setting,M,alpha,p,crit_mc){ # crit_mc=numeric(4) wieder Pflicht! 
  
  # Inputs: 
  # H0            | H0 distribution H0 
  # T             | length of time series T
  # setting       | "random" or numeric, if one wants to simulate a fixed number of VaR exceedances
  # M             | number of simulations to compute the significance
  # p             | the prescribed significance level
  # crit_mc       | optional. Needed if one wants to simulate a fixed number of VaR exceedances
  
  # Output:
  # significance  | a vector with calculated significances for all tests
  
  simulations = matrix(0,6,M)
  significance = numeric(6)
  
  # Calculation of VaR and ES
  var1 = var_es_analytic(H0,alpha)[[1]]
  es1 = var_es_analytic(H0,alpha)[[2]]
  var_es = matrix(0,T,2)
  var_es[,1] = rep(var1,T)
  var_es[,2] = rep(es1,T)
  
  cp_test = 0
  
  # Adjustment for the skewed distributions
  # The 3rd Test of A/S and the Tests of C/P calculate with values using the quantile function. 
  # The re-adjustment of the mean as described in the master thesis is already transfered here. 
  
  
  if(H0[[1]]=="skewn")
  {
    mu = as.numeric(H0[2])
    sigma = as.numeric(H0[3])
    shape = as.numeric(H0[4])
    mean_adjust = sigma*shape/sqrt(1+(shape)^2)*sqrt(2/pi)
    H0_adj = c("skewn",mu - mean_adjust,sigma,shape,0)
  }
  
  if(H0[[1]]=="skewt")
  {
    mu = as.numeric(H0[2])
    sigma = as.numeric(H0[3])
    df = as.numeric(H0[4])
    gamma = as.numeric(H0[5])
    mean_adjust = mean(sigma*rskt(100000,df,gamma))
    H0_adj = c("skewt",mu-mean_adjust,sigma,df,gamma)
  }
  # End Mean Adjustment
    
    # Calculate the significance for all tests
    # Simulation of M iid p&l time series and performing of the tests
    
    for (i in 1:M){
      if (is.numeric(setting) == FALSE){
      x = pl_iid_sim(T,H0)
      }
      if (is.numeric(setting) == TRUE){
      n = setting
      x = pl_iid_sim_fixn(T,H0,n,alpha)
      }
      simulations[1,i] = as1(x,var_es[,1],var_es[,2])
      simulations[2,i] = as2(x,var_es[,1],var_es[,2],alpha)
      if (H0[[1]]=="skewn" || H0[[1]]=="skewt"){
        simulations[3,i] = as3(x,H0_adj,alpha) # H0 Adjusted
      }
      if (H0[[1]]=="normal" || H0[[1]]=="stdt"){
        simulations[3,i] = as3(x,H0,alpha)    
      }
      simulations[4,i] = as4(x,var_es[,1],var_es[,2],alpha)
      if (H0[[1]]=="skewn" || H0[[1]]=="skewt"){
        cp_test = cp(x,var_es[,2],H0_adj,p) # H0 Adjusted  
      }
      if (H0[[1]]=="normal" || H0[[1]]=="stdt"){
        cp_test = cp(x,var_es[,2],H0,p) 
      }
      simulations[5,i] = cp_test[[3]]
      simulations[6,i] = cp_test[[4]]

    }
    
    for (k in 1:4)
    {
      significance[k] = 1-(length(simulations[k,][simulations[k,]<crit_mc[k]])/M)
    }
    
    significance[5] = 1-mean(simulations[5,])
    significance[6] = 1-mean(simulations[6,])
  
  return(significance)
  
}

# Power plot for the six backtests, see section "Power"
plot_power = function(seq,power,xlab,main){
  n = length(power[1,])
  plot(seq,power[,1], type="l", col="blue", lwd=2, xlab=xlab, ylab="power", ylim = range(0:1))
  lines(seq,power[,2], col="red", lwd=2)
  lines(seq,power[,3], col="green", lwd=2)
  lines(seq,power[,4], col="orange", lwd=2)
  lines(seq,power[,5], col="black", lwd=2)
  lines(seq,power[,6], col="grey", lwd=2)
  title(main=main)
  legend("bottomright", box.lty=0, cex=0.7, c("A/S 1", "A/S 2", "A/S 3", "A/S 4", "C/P 1", "C/P 2"), col=c("blue","red","green","orange","black","grey"),  lwd=2)
}

# Power: A fixed number of exceedances. Adjustments of the functions "pl_iid_sim" and "power_calc"

# Simulation of a fixed number of VaR exceedances under H0 by a given distribution H1
pl_iid_sim_fixn_H0H1 = function(T,var_H0,H1,n,alpha_H0){
  
  # Simulate a p&l time series with n exceedances under the VaR of alpha
  # Inputs:
  # T = length of time series
  # var_H0 = the VaR under H0
  # H1 = alternative hypothesis
  # n = number of returns under VaR 
  # alpha_H0 = confidence level
  
  x = numeric(T)
  outlier = numeric(n)
  rest = numeric(T-n)
  dist = H1[1]
  par1 = as.numeric(H1[2])
  par2 = as.numeric(H1[3])
  par3 = as.numeric(H1[4])
  par4 = as.numeric(H1[5])
  
  if (dist=="normal"){
    alpha_new = pnorm(-var_H0,par1,par2)  # Adjusted
    q1 = runif(n,0,alpha_new)            # Adjusted
    q2 = runif(T-n,alpha_new,1)          # Adjusted
    outlier = qnorm(q1,par1,par2)
    rest = qnorm(q2,par1,par2)
  }
  if (dist=="stdt"){
    alpha_new = pt(-var_H0/par2,par3)    # Adjusted
    q1 = runif(n,0,alpha_new)            # Adjusted
    q2 = runif(T-n,alpha_new,1)          # Adjusted
    outlier = qt(q1,par3)*par2
    rest = qt(q2,par3)*par2
  } 
  if (dist=="skewn"){
    shape = par3
    mu= par1 -shape/sqrt(1+(shape)^2)*sqrt(2/pi) # Mean adjustment
    alpha_new = psn(-var_H0,xi=mu,omega=par2,alpha=par3) # Adjusted
    q1 = runif(n,0,alpha_new)            # Adjusted
    q2 = runif(T-n,alpha_new,1)          # Adjusted
    outlier = qsn(q1,xi=mu,omega=par2,alpha=shape)
    rest = qsn(q2,xi=mu,omega=par2,alpha=shape)
  }
  if (dist=="skewt"){
    # Mean Adjustment
    sample = par2*rskt(T,df=par3,gamma=par4)
    mean_adjust = mean(sample)
    
    outlier = qskt(q1,par3,par4)*par2 - mean_adjust
    rest = qskt(q2,par3,par4)*par2 - mean_adjust
  }
  x = c(outlier,rest)
  return(x)
}

# Calculation of the power in that case
power_calc_fixn_H0H1 = function(T,n,crit_mc,var_es_H0,H0,H1,N,alpha,p){
  
  # Everything works analogue to the power function above.
  # The only difference is that here, a simulation of a fixed number of VaR_H0 exceedances is made.
  
  x = numeric(T)
  z_H1_mc = matrix(0,N,4) # N simulations, 4 tests of A/S
  cp_tests = matrix(0,N,2) # N simulations, 2 CP Tests
  power = numeric(6)
  
  # Adjustment for the skewed distributions
  # The 3rd Test of A/S and the Tests of C/P calculate with values using the quantile function. 
  # The re-adjustment of the mean as described in the master thesis is already transfered here. 
  
  
  if(H0[[1]]=="skewn")
  {
    mu = as.numeric(H0[2])
    sigma = as.numeric(H0[3])
    shape = as.numeric(H0[4])
    mean_adjust = sigma*shape/sqrt(1+(shape)^2)*sqrt(2/pi)
    H0_adj = c("skewn",mu - mean_adjust,sigma,shape,0)
  }
  
  if(H0[[1]]=="skewt")
  {
    mu = as.numeric(H0[2])
    sigma = as.numeric(H0[3])
    df = as.numeric(H0[4])
    gamma = as.numeric(H0[5])
    mean_adjust = mean(sigma*rskt(100000,df,gamma))
    H0_adj = c("skewt",mu-mean_adjust,sigma,df,gamma)
  }
  ################ 
  
  for (i in 1:N){
    
    # Calculate a Monte-Carlo sample of Z under H1 (with Var & ES from H0):
    x = pl_iid_sim_fixn_H0H1(T,var_es_H0[1],H1,n,alpha)
    z_H1_mc[i,1]=as1(x,var_es_H0[,1],var_es_H0[,2])
    z_H1_mc[i,2]=as2(x,var_es_H0[,1],var_es_H0[,2],alpha)
    if (H0[[1]]=="skewn" || H0[[1]]=="skewt"){
      z_H1_mc[i,3] = as3(x,H0_adj,alpha) # H0 Adjusted 
    }
    if (H0[[1]]=="normal" || H0[[1]]=="stdt"){
      z_H1_mc[i,3] = as3(x,H0,alpha)    
    }
    z_H1_mc[i,4]=as4(x,var_es_H0[,1],var_es_H0[,2],alpha)
    
    # Test C/P directly:
    if (H0[[1]]=="skewn" || H0[[1]]=="skewt"){
      cp_test = cp(x,var_es_H0[,2],H0_adj,p) # H0 Adjusted  
    }
    if (H0[[1]]=="normal" || H0[[1]]=="stdt"){
      cp_test = cp(x,var_es_H0[,2],H0,p)   
    }
    cp_tests[i,1]=cp_test[[3]]
    cp_tests[i,2]=cp_test[[4]]
    
  }
  
  # Calculate the power of the tests: 
  power[1]=length(z_H1_mc[,1][z_H1_mc[,1]<crit_mc[1]])/N
  power[2]=length(z_H1_mc[,2][z_H1_mc[,2]<crit_mc[2]])/N
  power[3]=length(z_H1_mc[,3][z_H1_mc[,3]<crit_mc[3]])/N
  power[4]=length(z_H1_mc[,4][z_H1_mc[,4]<crit_mc[4]])/N
  power[5]=mean(cp_tests[,1])
  power[6]=mean(cp_tests[,2])
  
  return(power)
  
}

# Which role does VaR play? Adjustments of the functions "crit_H0" and "significance"

# Adjustement of crit_H0 for example 1
crit_H0_example1 = function(T,var_es_H0,M,alpha,p,denom){
  
  # Similar structure as the function "crit_H0", adjusted for the example1 (VaR-sensitivity)
  
  x = numeric(T)
  z_H0_mc = matrix(0,M,4) # M simulations, 4 tests
  crit_mc = numeric(4)
  
  # Calculate a Monte-Carlo sample of Z under H0:
  for (i in 1:M){
    
    x = example1(T,alpha)
    z_H0_mc[i,1]=as1(x,var_es_H0[,1],var_es_H0[,2])
    z_H0_mc[i,2]=as2(x,var_es_H0[,1],var_es_H0[,2],alpha)
    z_H0_mc[i,3]=as3_adjusted(x,alpha,denom)
    z_H0_mc[i,4]=as4(x,var_es_H0[,1],var_es_H0[,2],alpha)
  }
  
  # Calculate the critical values: 
  z_H0_mc[,1]=sort(z_H0_mc[,1],na.last=TRUE) 
  z_H0_mc[,2]=sort(z_H0_mc[,2])
  z_H0_mc[,3]=sort(z_H0_mc[,3])
  z_H0_mc[,4]=sort(z_H0_mc[,4])
  crit_mc[1]=z_H0_mc[,1][floor(M*p)]
  crit_mc[2]=z_H0_mc[,2][floor(M*p)]  
  crit_mc[3]=z_H0_mc[,3][floor(M*p)]
  crit_mc[4]=z_H0_mc[,4][floor(M*p)]
  
  return(crit_mc)
  
}

# Adjustment of the function "significance" for example 1
significance_example1 = function(var_es_H0,H1,T,M,alpha,p,crit_mc,denom){
  
  simulations = matrix(0,6,M)
  significance = numeric(6)
  
  var_es = matrix(0,T,2)
  var_es[,1] = var_es_H0[,1]
  var_es[,2] = var_es_H0[,2]
  
  as_values = matrix(0,M,4)
  
  for (i in 1:M){
    
    x = pl_iid_sim(T,H1) 
    simulations[1,i] = as1(x,var_es[,1],var_es[,2])
    simulations[2,i] = as2(x,var_es[,1],var_es[,2],alpha)
    simulations[3,i] = as3_adjusted(x,alpha,denom)
    simulations[4,i] = as4(x,var_es[,1],var_es[,2],alpha)
    cp_test = cp_adjusted(x,var_es[,2],cp_pt,p)
    simulations[5,i] = cp_test[[3]]
    simulations[6,i] = cp_test[[4]]
  }
  
  for (k in 1:4)
  {
    significance[k] = 1-(length(simulations[k,][simulations[k,]<crit_mc[k]])/M)
  }
  
  significance[5] = 1-mean(simulations[5,])
  significance[6] = 1-mean(simulations[6,])
  
  
  return(significance)
  
  
}

# Adjustment of the function "significance" for example 1b (reversed distribution)
significance_example1_vv = function(var_es_H0,H0,T,M,alpha,p,crit_mc){
  
  simulations = matrix(0,6,M)
  significance = numeric(6)
  
  var_es = matrix(0,T,2)
  var_es[,1] = var_es_H0[,1]
  var_es[,2] = var_es_H0[,2]
  
  for (i in 1:M){
    
    x = example1_vv(T,alpha) 
    simulations[1,i] = as1(x,var_es[,1],var_es[,2])
    simulations[2,i] = as2(x,var_es[,1],var_es[,2],alpha)
    simulations[3,i] = as3(x,H0,alpha) 
    simulations[4,i] = as4(x,var_es[,1],var_es[,2],alpha)
    cp_test = cp(x,var_es[,2],H0,p)
    simulations[5,i] = cp_test[[3]]
    simulations[6,i] = cp_test[[4]]
    
  }
  
  for (k in 1:4)
  {
    significance[k] = 1-(length(simulations[k,][simulations[k,]<crit_mc[k]])/M)
  }
  
  significance[5] = 1-mean(simulations[5,])
  significance[6] = 1-mean(simulations[6,])
  
  return(significance)
}

# GARCH Setup: 

# Own GARCH path simulation: 
own_garch_path <- function(T,omega,alpha,beta,gamma,model,dist,mean_model=c(),mu=c(),AR_alpha=c(),MA_beta=c(),dist_para=c()){
  
  # Inputs:
  # T          | Length of time series
  # omega      | Omega parameters of the GARCH model
  # alpha      | Alpha parameters of the GARCH model
  # beta       | Beta parameters of the GARCH model
  # gamma      | Gamma parameters of the GARCH model (EGARCH and TGARCH)
  # model      | garch/egarch/tgarch
  # dist       | Distribution for the innovations
  # mean_model | conditional, AR, MA, ARMA, TAR, GarchInMean
  # mu         | Initial values 
  # AR_alpha   | Alpha in the mean model
  # MA_beta    | Beta in the mean model
  # dist_para  | scaling parameter if dist="std"
  
  # Output
  # y          | GARCH path
  
  q = length(alpha);
  p = length(beta);
  if(is.null(mean_model)==FALSE){
    AR_p = length(AR_alpha)
    MA_q = length(MA_beta)
  }
  
  eps=c();
  garch_sig2=c();
  if(is.null(mean_model)==FALSE){
    arma_mean=c();
    y=c();
  }
  
  sig = omega;
  # Innovations:
  if (dist=="norm"){
    nu=rnorm(T,0,1);
    EW=sqrt(2/pi);
  }
  if (dist=="std"){
    if (is.null(dist_para)==TRUE){
      nu=rstd(T,0,1); # The variance is scaled to 1 (fGarch package)
    } else {
      nu=rstd(T,0,1,nu=dist_para)
    }
    EW=mean(abs(nu));
  }
  if(dist=="sstd"){
    if (is.null(dist_para)==TRUE){
      nu=rsstd(T,0,1);
    } else {
      nu=rsstd(T,0,1,nu=dist_para);
    }
    EW=mean(abs(nu));
  }
  if(dist=="mixednormal"){
    nu=(0.6*rnorm(T,1,sqrt(2))+0.4*rnorm(T,-1.5,sqrt(0.75)))/sqrt(3)
  }
  # End Innovations
  
  if (model=="egarch"){
    
    eps[1] = sqrt(exp(sig))*nu[1];
    garch_sig2[1] = exp(sig);
  }else{
    eps[1] = sqrt(sig)*nu[1];
    garch_sig2[1] = sig;
  }
  
  
  for (i in 2:T){
    sig2 = sig;
    if(is.null(mean_model)==FALSE){
      arma_mean2 = mu;
    }
    if (q>0){
      for (j in 1:q){
        if (i>j) { # Interception of the initial cases
          if (model =="garch"){sig2 = sig2 + alpha[j]*eps[i-j]^2;}
          if (model=="egarch"){sig2 = sig2 + alpha[j]*nu[i-j] + gamma[j]*(abs(nu[i-j])-EW);}
          if (model=="tgarch"){if (eps[i-j]<0) {sig2 = sig2 + (alpha[j]+gamma[j])*eps[i-j]^2 ;
          } else {sig2 = sig2 + alpha[j]*eps[i-j]^2;}}
        }
      }
    }
    if (p>0){
      for (k in 1:p){
        if (i>k) {
          if (model=="garch"|model=="tgarch"){sig2 = sig2 + beta[k]*garch_sig2[i-k];}
          if (model=="egarch"){sig2 = sig2 + beta[k]*log(garch_sig2[i-k]);}
        }
      }
    }
    
    # Set garch_sig2 as final
    if (model=="egarch"){eps[i] = sqrt(exp(sig2))*nu[i]; garch_sig2[i] = exp(sig2);}
    else{eps[i] = sqrt(sig2)*nu[i];garch_sig2[i] = sig2;}
    
    # Calculation of the Mean
    if(is.null(mean_model)==FALSE){
      arma_mean[1] = mu;
      y[1]= arma_mean[1] + eps[1];
      if (mean_model == "unconditional"){
        arma_mean2 = mu
      } else if(mean_model == "AR"){
        for (l in 1:AR_p){
          if (i>l){
            arma_mean2 = arma_mean2 + AR_alpha[l]*y[i-l];
          }
        }
      } else if (mean_model == "MA"){
        for (m in 1:MA_q){
          if (i>m){
            arma_mean2 = arma_mean2 + MA_beta[m]*eps[i-m]
          }
        }
      } else if (mean_model == "ARMA"){
        for (l in 1:AR_p){
          if (i>l){
            arma_mean2 = arma_mean2 + AR_alpha[l]*y[i-l];
          }
        }
        for (m in 1:MA_q){
          if (i>m){
            arma_mean2 = arma_mean2 + MA_beta[m]*eps[i-m]
          }
        }
      } else if (mean_model == "TAR"){
        arma_mean2=0.7*(eps[i-1]<=-2)*y[i-1]
      } else if (mean_model == "GarchInMean"){
        arma_mean2=-0.5*garch_sig2[i] 
      }
      arma_mean[i] = arma_mean2
      y[i] = arma_mean[i]+eps[i]
    }
  }
  
  return("Path"=y)
}  

# Function to calculate the critical values of A/S 1-4 under H0
garch_crit_mc = function(T,M,p){
  
  mc_as = matrix(0,M,4); crit_mc = numeric(4); denom= numeric(T-1)
  for (i in 1:M){
    # Simulation of one path under H0:
    x = own_garch_path(T,omega=0.05,alpha=0.1,beta=0.85,gamma=0,model="garch",dist="norm",mean_model="AR",mu=0,AR_alpha=c(0.05),MA_beta=c())  
    path = x
    
    # Forecasting with x under H0: 
    forecast = ugarchforecast(H0,path,n.ahead=1,n.roll=T-2,out.sample=T-1)
    # Predicted Sigma for x=2, ..., x=T
    predicted_sigma = as.numeric(sigma(forecast))
    predicted_mu = as.numeric(fitted(forecast))
    # Calculation of VaR and ES: 
    var_es = matrix(0,T-1,2)
    for (j in 1:(T-1)){
      H0_dist = c("normal",predicted_mu[j],predicted_sigma[j],0,0)
      var_es[j,]=var_es_analytic(H0_dist,alpha)
    }
    # Calculate the test statistics for A/S 1,2,4
    mc_as[i,1]=as1(path[2:T],var_es[,1],var_es[,2])
    mc_as[i,2]=as2(path[2:T],var_es[,1],var_es[,2],alpha)
    mc_as[i,4]=as4(path[2:T],var_es[,1],var_es[,2],alpha)
    
    # Preparation of the denominator for A/S 3
    for (j in 1:(T-1)){
      summe = function(x){
        return(fb(x)*qnorm(x,predicted_mu[j],predicted_sigma[j]))
      }
      denom_help = integrate(summe,0,1)
      denom[j] = denom_help[[1]]
      denom[j] = -T/(floor(T*alpha)) * denom[j]
    }
    # Calculate the test statistic for A/S 3: 
    mc_as[i,3]=as3_adjusted(path[2:T],alpha,denom)
  }
  
  # Calculate the critical values:
  crit_mc[1] = sort(mc_as[,1])[floor(M*p)]
  crit_mc[2] = sort(mc_as[,2])[floor(M*p)]
  crit_mc[3] = sort(mc_as[,3])[floor(M*p)]
  crit_mc[4] = sort(mc_as[,4])[floor(M*p)]
  
  return(crit_mc)
  
}

# Function to calculatie the rejection rate of the alternative vs H0 
garch_rejection_rate = function(alternative,T,N,crit_mc,alpha,p){
  
  power = numeric(6);
  simulations = matrix(0,N,9); 
  var_001 = numeric((T-1))
  
  
  for (i in 1:N){
    
    # Simulation of one path under each alternative:
    # If ugarchspec is used
    # H1 = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),  mean.model = list(armaOrder = c(1,0), include.mean=FALSE), distribution.model ="norm", fixed.pars = list(ar1 = 0.05, omega = 0.05, alpha1 = 0.1, beta1 = 0.85))
    # x = ugarchpath(H1,n.sim=T,n.start=0,m.sim=1)
    # path = x@path$seriesSim
    
    # Simulation of a path under the alternative distribution
    
    if (alternative=="H0"){
      x = own_garch_path(T,omega=0.05,alpha=0.1,beta=0.85,gamma=0,model="garch",dist="norm",mean_model="AR",mu=0,AR_alpha=c(0.05),MA_beta=c())  
    }
    if (alternative=="A1"){
      x = own_garch_path(T,omega=0.04,alpha=0.1,beta=0.89,gamma=0,model="garch",dist="norm",mean_model="TAR",mu=0,AR_alpha=c(),MA_beta=c())
    }
    if (alternative=="A2"){
      x = own_garch_path(T,omega=0.01,alpha=0.29,beta=0.7,gamma=0,model="garch",dist="norm",mean_model="GarchInMean",mu=0,AR_alpha=c(),MA_beta=c())    
    }
    if (alternative=="A3"){
      x = own_garch_path(T,omega=0.01,alpha=c(0.1,0.8),beta=c(),gamma=0,model="garch",dist="norm",mean_model="garch",mu=0,AR_alpha=c(0.05),MA_beta=c())
    }
    if (alternative=="A4"){
      x = own_garch_path(T,omega=0.05,alpha=0.8,beta=0.9,gamma=-0.3,model="egarch",dist="norm",mean_model="AR",mu=0,AR_alpha=c(0.05),MA_beta=c())  
    }
    if (alternative=="A5"){
      x = own_garch_path(T,omega=0.05,alpha=0.1,beta=0.85,gamma=0,model="garch",dist="mixednormal",mean_model="AR",mu=0,AR_alpha=c(0.05),MA_beta=c())  
    }
    
    path=x
    
    # forecasting with x under H0: 
    forecast = ugarchforecast(H0,path,n.ahead=1,n.roll=T-2,out.sample=T-1)
    # forecast2 = ugarchforecast(H0,path2,n.ahead=1,n.roll=T-2,out.sample=T-1)
    # Predicted Sigma for x=2, ..., x=T
    predicted_sigma = as.numeric(sigma(forecast))
    predicted_mu = as.numeric(fitted(forecast))
    
    # Calculation of VaR and ES (alpha)
    var_es = matrix(0,T-1,2)
    for (j in 1:(T-1)){
      H0_dist = c("normal",predicted_mu[j],predicted_sigma[j],0,0)
      var_es[j,]=var_es_analytic(H0_dist,alpha)
    }
    # Calculation of the test statistics Z1, .., Z4: 
    simulations[i,1]=as1(path[2:T],var_es[,1],var_es[,2])
    simulations[i,2]=as2(path[2:T],var_es[,1],var_es[,2],alpha)
    simulations[i,4]=as4(path[2:T],var_es[,1],var_es[,2],alpha)
    
    # Preparation of the denominator for A/S 3: 
    nenner=numeric(T-1)
    for (j in 1:(T-1)){
      summe = function(x){
        return(fb(x)*qnorm(x,0,predicted_sigma[j]))
      }
      nenner_help = integrate(summe,0,1)
      nenner[j] = nenner_help[[1]]
      nenner[j] = -T/(floor(T*alpha)) * nenner[j]
    }
    
    simulations[i,3]=as3_adjusted(path[2:T],alpha,nenner)
    
    # Preparation of the Corbetta and Peri Input: 
    
    cp_pt = numeric(T-1)
    # Is the same under H0 for each t! 
    for (j in 1:(T-1)){
      cp_pt[j] = pnorm(-var_es[j,2],0,predicted_sigma[j]) 
    }
    
    simulations[i,5]=cp_adjusted(path[2:T],var_es[,2],cp_pt,p)[[3]]
    simulations[i,6]=cp_adjusted(path[2:T],var_es[,2],cp_pt,p)[[4]]
    
  }  
  
  # Calculation of the rejection rate
  power[1]=length(simulations[,1][simulations[,1]<crit_mc[1]])/N
  power[2]=length(simulations[,2][simulations[,2]<crit_mc[2]])/N
  power[3]=length(simulations[,3][simulations[,3]<crit_mc[3]])/N
  power[4]=length(simulations[,4][simulations[,4]<crit_mc[4]])/N
  power[5]=mean(simulations[,5])
  power[6]=mean(simulations[,6])
  
  return(power)
  
}
