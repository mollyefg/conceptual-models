require(tseriesChaos)

### simple conceptual model of long-term jack pine budworm population dynamics
##  N: jack pine budworm population
##  Z: parasitoid population


# define functions
Nfunc = function(N, beta, Z){
  	lambda = 20
	Ntmp = lambda*N*exp(-(beta*Z)/(1+Z));
  
  return(Ntmp)
}

Zfunc = function(omega, N, beta, Z){
  Ztmp = omega*N*(1-exp(-(beta*Z)/(1+Z)))
  return(Ztmp)
}


#  set up parameters, initial conditions, and run time

beta = 4.4; #5 #650 ; #1 ; #2.0; #0.1; #0.001; 
omega = 100

n0 = .9 # 500
z0 = .8; #50

nstor = c()
zstor = c()

maxyrs = 10000
years = seq(1,maxyrs,1)

## run

for(i in 1:max(years)){
  n = Nfunc(N=n0, beta=beta, Z=z0)
    
  z = Zfunc(omega=omega, N=n0, beta=beta, Z=z0)

  n0 = n
  nstor[i] = n0
  
  z0 = z
  zstor[i] = z0
	}

## organize results and plot

count = length(nstor[nstor!=Inf])
par(mfrow=c(1,3))
fract = 0.99;

start = round(count*fract)

par(ask=TRUE)
plot(log10(nstor[start:count]),log10(zstor[start:count]),cex=0.01, xlab = "log(budworms)", ylab = "log(parasitoids")
plot(years[start:count], log10(nstor[start:count]),  pch = 20, type="l", xlab = "years", ylab = "log(budworms)")
plot(years[start:count], log10(zstor[start:count]), col = "red", pch = 20,type = "l",new="FALSE",xlab = "log(year)", ylab="log(parasitoid)")


