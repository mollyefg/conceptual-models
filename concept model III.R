### conceptual long-term model of jack pine budworm population dynamics
## most complex: includes host, parasitoid, and three size classes of jack pine trees

##### N = jpbw, Z = parasitoid population, S = small trees, M = medium trees, L = large trees

# define functions

Nfunc = function(N, S, M, L, a,b,c,K,gamma, beta, Z){
  lambda = 5.0*((a*S+b*M+c*L)/(K+S+M+L))
	Ntmp = lambda*N*exp(-(beta*Z)/(gamma+Z));

	return(Ntmp)
	}

Zfunc = function(omega, Z, gamma,beta, N, a, S, b, M, c, L, K){
	Ztmp = omega*N*(1-exp(-(beta*Z)/(gamma+Z)))
	return(Ztmp)
	}

Sfunc = function(r1, M, r2, L, S, g1, d1, N, K){
	Stmp = (r1*M+r2*L)/(K+M+L) + S*exp(-g1*S)*exp(-d1*N)
	return(Stmp)
	}

Mfunc = function(M, S, g1, g2, d2, N){
	Mtmp = S*(1-exp(-g1*S)) + M*exp(-g2*M)*exp(-d2*N)
	return(Mtmp)
	}

Lfunc = function(M,g2,L,g3,d3,N){
	Ltmp = M*(1-exp(-g2*M)) + L*exp(-g3*L)*exp(-d3*N)
	return(Ltmp)
	}


#  set up parameters, initial conditions, and run time

# r1, r2: births of trees
# g1 = 1/av'g time in small class
# g2 = 1/av'g time in med class
# g3 = 1/av'g time in large class
# d1, d2, d3 = chance of dying in each age class each year

r1 = 0.1; r2 = 5
g1 = 1/10; g2 = 1/10; g3 = 1/20
d1 = 1/200; d2 = 1/10000; d3 = 1/50

d1 = d2 = d3 = 0.01;

a = 1; b = 1; c = 1;

K = 0.2
K2 = 0.5;
beta = 6.5 ; #1 ; #2.0; #0.1; #0.001; 
gamma = 1; #0.5; #0.005; #.0005; 
omega = 1

# assign values for: N, Z, S, M, L, lambda
n0 = .41 # 500
z0 = .33; #50
s0 = 100
m0 = 500
l0 = 30

nstor = c()
zstor = c()
sstor = c()
mstor = c()
lstor = c()

lambstor = c()

maxyrs = 5000
years = seq(1,maxyrs,1)


## run


for(i in 1:max(years)){
	#Nfunc = function(N, S, M, L, a,b,c,K,gamma, beta,  Z){
	n = Nfunc(N=n0, S=s0, M=m0, L=l0, a=a, b=b, c=c, K=K, gamma=gamma, beta=beta, Z=z0)
	z = Zfunc(omega=omega, Z=z0, gamma=gamma, beta = beta, N = n0, a = a, S = s0, b = b, M = m0, c = c, L = l0, K = K)
	s = Sfunc(r1=r1, M=m0, r2=r2, L=l0, S=s0, g1=g1, d1=d1, N=n0, K=K2)
	m = Mfunc(M=m0, S=s0, g1=g1, g2=g2, d2=d2, N=n0)
	l = Lfunc(M=m0, g2=g2, L=l0, g3=g3, d3=d3, N=n0)
  
	n0 = n
	nstor[i] = n0
  
	#Zfunc = function(omega, Z, gamma,beta, N, a, S, b, M, c, L, K){
	z0 = z
	zstor[i] = z0
	  
	#Sfunc = function(r1, M, r2, L, S, g1, d1, N){
	
	s0 = s
	sstor[i] = s0
	  
	#Mfunc = function(M, S, g1, g2, d2, N){
	
	m0 = m
	mstor[i] = m0
	  
	#Lfunc = function(M,g2,L,g3,d3,N){
	
	l0 = l
	lstor[i] = l0

	i = i+1
	}


## organize results and plot

count = length(nstor[nstor!=Inf])
par(mfrow=c(2,3))

fract = 0.97;
start = round(count*fract)
par(ask=TRUE);

plot(log10(nstor[start:count]),log10(zstor[start:count]),cex=2, xlab = "log(budworms)", ylab = "log(parasitoids)", pch = 20)
plot(years[start:count], log10(nstor[start:count]), col = "dark green", pch = 20, type="l", xlab = "years", ylab = "log(budworms)")
plot(years[start:count], log10(zstor[start:count]), col = "red", pch = 20,type = "l",new="FALSE", xlab = "years", ylab = "log(parasitoids)")
plot(years[start:count], log10(sstor[start:count]), col = "red", pch = 20,type = "l",new="FALSE", xlab = "years", ylab = "log(small trees)")
plot(years[start:count], log10(mstor[start:count]), col = "red", pch = 20,type = "l",new="FALSE", xlab = "years", ylab = "log(medium trees)")
plot(years[start:count], log10(lstor[start:count]), col = "red", pch = 20,type = "l",new="FALSE", xlab = "years", ylab = "log(large trees)")


## power spectrum

par(mfrow=c(1,1))
z = spectrum(log10(nstor[start:(start+round((maxyrs-start)/3))]),plot=FALSE); #Calculating the power spectrum
plot((1/z$freq),sqrt(z$spec),xlim = c(0,20),type="o",xlab="Period",ylab="Power", main = "Power spectrum");



