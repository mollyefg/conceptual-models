### conceptual long-term model of jack pine budworm population dynamics
## includes host, parasitoid, and jack pine trees

##### N = jpbw, Z = parasitoid population, T=trees

# define functions

lambfunc = function(a,T,K){
	ltmp=5.0*((a*T)/(K+T))
	return(ltmp)
	}
	

Nfunc = function(N,T, a,K,gamma, beta, Z, lambda){
  #lambda = 5.0*((a*T)/(K+T))
  Ntmp = lambda*N*exp(-(beta*Z)/(gamma+Z));
  
  return(Ntmp)
}

Zfunc = function(omega, Z, gamma,beta, N, K){
  Ztmp = omega*N*(1-exp(-(beta*Z)/(gamma+Z)))
  return(Ztmp)
}

Tfunc = function(r, T, g, d, N, K){
  Ttmp = (r*T)/(K+T) + T*exp(-g*T)*exp(-d*N)
  return(Ttmp)
}


#  set up parameters, initial conditions, and run time

# r1, r2: births of trees
# g1 = 1/av'g time in small class
# g2 = 1/av'g time in med class
# g3 = 1/av'g time in large class
# d1, d2, d3 = chance of dying in each age class each year

r = 1
g = 1/100
d = 1/100

a=2

K = 0.2
K2 = 0.5;
beta = 3.9; #6.5 ; #1 ; #2.0; #0.1; #0.001; 
gamma = 1; #0.5; #0.005; #.0005; 
omega = 1

# assign values for: N, Z, T

n0 = .41 # 500
z0 = .33; #50
t0 = 200
l0 = 5.0*((a*t0)/(K+t0))

nstor = c()
zstor = c()
tstor = c()
lambstor = c()

maxyrs = 10000; #2500
years = seq(1,maxyrs,1)

## run


for(i in 1:max(years)){
  l = lambfunc(a=a, T=t0, K=K)

  n = Nfunc(N=n0, T = t0, a=a, K=K, gamma=gamma, beta=beta, Z=z0, lambda=l0)
    
  z = Zfunc(omega=omega, Z=z0, gamma=gamma, beta = beta, N = n0, K = K)

  t = Tfunc(r = r, T=t0, g = g, d = d, N = n0, K = K2)
  
  l0 = l
  lambstor[i] = l0

  n0 = n
  nstor[i] = n0
  
  z0 = z
  zstor[i] = z0
    
  t0 = t
  tstor[i] = t0
}

## organize results and plot

count = length(nstor[nstor!=Inf])
par(mfrow=c(2,2))

fract = 0.99;
start = round(count*fract)
p=nstor[start:count]
plot(log10(nstor[start:count]),log10(zstor[start:count]), pch = 20, xlab = "log(host)", ylab = "log(parasitoid", col = rainbow(length(p))[rank(p)])
plot(years[start:count], log10(nstor[start:count]), col = "dark green", pch = 20, type="l", xlab = "log(year)", ylab = "log(host)")
plot(years[start:count], log10(zstor[start:count]), col = "red", pch = 20,type = "l",new="FALSE",xlab = "log(year)", ylab="log(parasitoid)")
plot(years[start:count], log10(tstor[start:count]), col = "red", pch = 20,type = "l",new="FALSE",xlab="log(year)",ylab = "log(trees)")



