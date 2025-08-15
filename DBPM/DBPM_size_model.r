sizemodel<-function(params,ERSEM.det.input=F,pp,bio.det,sst,sft,U_mat,V_mat,W_mat,temp.effect=T,eps=1e-5){



with(params, {

#---------------------------------------------------------------------------------------
# Model for a dynamical ecosystem comprised of: two functionally distinct size spectra (predators and detritivores), size structured primary producers and an unstructured detritus resource pool. 
# time implicit upwind difference discretization algorithm (from Ken Andersen)
# (Numerical Recipes and Richards notes, see "C://...standardised test implicit.txt" for basic example)
# U represents abundance density (m-2) of "fish" and V abundance density (m-2) of "detritivores"
# W is the pool of detritus (expressed in biomass density, g.m-2) - not size-based
# Fish feed on other fish and benthic detritivores with size preference function
# Benthic detritivores feed on detritus produced from pelagic spectrum 
# Senescence Mortality also included.  Options to include dynamic reproduction and predator handling time (but not currently used).
#
# Code modified for global fishing mortality rate application. JLB 17/02/2014
# ---------------------------------------------------------------------------------------
# Input parameters to vary:
 		
ui0 <- 10^pp          # intercept of plankton size spectrum (estimated from biogeophysical model output or satellite data).		
W.init <- bio.det	  # initial detritus pool biomass density 
sst <- sst            # temperature in water column
sft <- sft		      # temperature near seabed



#---------------------------------------------------------------------------------------
# Functions called within sizemodel()
#---------------------------------------------------------------------------------------


#Function to build a lookup table for diet preference of all combinations of predator 
#and prey body size: diet preference (in the predator spectrum only)

phi.f=function(q){
	q=q
	q0=q0
	sd.q=sd.q

	phi=ifelse(q>0,exp(-(q-q0)*(q-q0)/(2*sd.q*sd.q))/(sd.q*sqrt(2.0*pi)),0) 
	
	return(phi)
	
	# or normalise feeding kernel to sum to 1:
	# return(phi/sum(phi))
    # this is commented out because in previous work normalising did not produce realistic growth rates
	

}

#Functions to build lookup tables for components of integration which remain constant
gphi.f = function(q1) { return( 10^(-q1)* phi.f(q1 )) }	#growth
mphi.f = function(q2) { return( 10^(alpha.u*q2)*phi.f(q2)) }	#mortality

#Function to build lookup table for components of 10^(alpha*x)
expax.f = function(x, alpha.u) { return( 10^(alpha.u*x )) }

#Function to compute convolution products: can be for growth or mortality depending on phi
convolution.f = function(phi, u) {
    res = matrix(NA,length(x),1)
    for (i in 1:length(res)) {
	res[i] = 0.0
	res[i] = sum(phi[,i] * u * dx)
    }
    return(res)
}

#--------------------------------
# Growth and Mortality Equations
#--------------------------------

#predators
#death.u<-function(u,A,expax,mphi){return(pref.pel*A*expax*(convolution.f(mphi,u)))}
death.u<-function(u,A,expax,mphi){return((pref.pel*A*expax)*(u*dx)%*%(mphi))}         #use faster matrix method instead of convolution loop function
#growth.u<-function(K0,K1,A,expax,gphi,u,v){return((pref.pel*K0*A*expax*(convolution.f(gphi,u)))+ (pref.ben*K1*A*expax*(convolution.f(gphi,v))))}
growth.u<-function(K0,K1,A,expax,gphi,u,v){return((pref.pel*K0*A*expax)*(u*dx)%*%(gphi)+ (pref.ben*K1*A*expax)*(v*dx)%*%(gphi))}
#detritivores
#death.v<-function(u,A,expax,mphi){return((pref.ben)*A*expax*(convolution.f(mphi,u)))}
death.v<-function(u,A,expax,mphi){return((pref.ben*A*expax)*(u*dx)%*%(mphi))}
growth.v<-function(K1,A,w){return((1/10^x)*(K1*A*10^(x*0.75)*w))}

#detritus output (g.m-3.yr-1)
out.w<-function(A,v,w){return(sum((A*(10^(x*0.75))*w*v)*dx))}     



#---------------------------------------------------------------------------------------
# Initialising matrices
#---------------------------------------------------------------------------------------

#q1 is a square matrix holding the log(predatorsize/preysize) for all combinations of sizes
y = x

q1  = matrix(NA, length(x), length(y))
for (i in 1:length(y)) { q1[,i] = y[i] - x}

#q2 is the reverse matrix holding the log(preysize/predatorsize) for all combinations of sizes
q2 = matrix(-q1, length(x), length(y))	

#matrix for recording the two size spectra 
V = U = array(0, c(length(x), Neq+1))

#vector to hold detrtitus biomass density (g.m-2)
W = array(0,Neq)

#matrix for keeping track of growth and reproduction from ingested food:
R.v=R.u=GG.v = GG.u   = array(0, c(length(x), Neq+1)) 

#matrix for keeping track of predation mortality
PM.v = PM.u   = array(0, c(length(x), Neq+1))   

#matrix for keeping track of catches  
Y.v = Y.u  = array(0, c(length(x), Neq+1))  

#matrix for keeping track of  total mortality (Z)
Z.v = Z.u = array(0, c(length(x), Neq+1))

#matrix for keeping track of senescence mortality and other (intrinsic) mortality
SM.v = SM.u   =OM.v = OM.u   = array(0, length(x))

#empty vector to hold fishing mortality rates at each size class
Fvec.v = Fvec.u = rep(0,length(x))



#lookup tables for terms in the integrals which remain constant over time
gphi  = gphi.f(q1)
mphi  = mphi.f(q2)

#lookup table for components of 10^(alpha.u*x)
expax = expax.f(x, alpha.u)


#---------------------------------------------------------------------------------------
# Numerical integration
#---------------------------------------------------------------------------------------
#
#INITIAL VALUES
if (Neq>2) {
   # if running to equilibrium use given initial values
   U[1:(ref-1),1]<-ui0*10^(r.plank*x)[1:(ref-1)]    #(phyto+zoo)plankton size spectrum  
   U[ref:Nx,1]<-U.init[ref:Nx]          # set initial consumer size spectrum 
   V[ref.det:Nx,1]<-sinking.rate*V.init[ref.det:Nx]  # set initial detritivore spectrum  
   W[1]<-W.init             # set initial detritus biomass density (g.m^-3) 
   

} else {
   
   # use values from last time step
   
   U[1:ref,1]<-ui0*10^(r.plank*params$x)[1:ref]   #(phyto+zoo)plankton size spectrum  
   U[ref:Nx,1]<-U_mat[ref:(Nx)]      # read in initial consumer size spectrum 
   V[ref.det:Nx,1]<-V_mat[ref.det:Nx]  # read in initial detritivore spectrum  
   W[1]<-W_mat             # read in initial detritus biomass density (g.m^-3)
   
   
}

#intrinsic natural mortality
OM.u<-mu0*10^(-0.25*x)
OM.v<-mu0*10^(-0.25*x)

#senescence mortality rate to limit large fish from building up in the system
#same function as in Law et al 2008, with chosen parameters gives similar M2 values as in Hall et al. 2006
SM.u=k.sm*10^(p.s*(x-xs))
SM.v=k.sm*10^(p.s*(x-xs))

#Fishing mortality (THESE PARAMETERS NEED TO BE ESTIMATED!)

#Fvec[Fref:Nx] = 0.09*(x[Fref:Nx]/ log10(exp(1)) ) + 0.04 # from Benoit & Rochet 2004 
 
Fvec.u[Fref.u:Nx] = rep(Fmort.u,length(Fvec.u[Fref.u:Nx]))  
Fvec.v[Fref.v:Nx] = rep(Fmort.v,length(Fvec.v[Fref.v:Nx]))  



#iteration over time, N [days]
for (i in 1:(Neq)) {

if(W[i]=="NaN"|W[i]<0)
{
	#browser()
	U[,i]<-"NaN"
	return(list(U=U[,],GG.u=GG.u[,],PM.u=PM.u[,],V=V[,],GG.v=GG.v[,],PM.v=PM.v[,],Y.u=Y.u[,],Y.v=Y.v[,],W=W[], params=params))	
}
if(i>100)
{
#browser()
if(max(abs((log(U[-c(1:91),(i-1)]))-(log(U[-c(1:91),(i-2)]))))<eps
	&max(abs((log(V[-c(1:91),(i-1)]))-(log(V[-c(1:91),(i-2)]))))<eps
	&(abs((log(W[(i-1)]))-(log(W[(i-2)]))))<eps)
{
	U[,Neq]<-U[,i-1]
	PM.v[,Neq]<-PM.v[,i-1]
	GG.v[,Neq]<-GG.v[,i-1]
	V[,Neq]<-V[,i-1]
	PM.u[,Neq]<-PM.u[,i-1]
	GG.u[,Neq]<-GG.u[,i-1]
	Y.u[,Neq]<-Y.u[,i-1]
	Y.v[,Neq]<-Y.v[,i-1]
	W[Neq]<-W[i-1]
	return(list(U=U[,],GG.u=GG.u[,],PM.u=PM.u[,],V=V[,],GG.v=GG.v[,],PM.v=PM.v[,],Y.u=Y.u[,],Y.v=Y.v[,],W=W[], params=params))	
}
}

#--------------------------------
# Calculate Growth and Mortality
#--------------------------------

if (temp.effect==T) {
pel.Tempeffect=exp(c1 - E/(k*(sst+ 273)))
ben.Tempeffect=exp(c1 - E/(k*(sft+ 273)))
}            
                        
if (temp.effect==F) {
pel.Tempeffect=1
ben.Tempeffect=1
}            
      

# feeding rates

f.pel<-pel.Tempeffect*as.vector(((A.u*10^(x*alpha.u)*pref.pel)*(U[,i]*dx)%*%(gphi))/(1+handling*(A.u*10^(x*alpha.u)*pref.pel)*(U[,i]*dx)%*%(gphi))) # yr-1

f.ben<-pel.Tempeffect*as.vector(((A.u*10^(x*alpha.u)*pref.ben)*(V[,i]*dx)%*%(gphi))/(1+handling*(A.u*10^(x*alpha.u)*pref.ben)*(V[,i]*dx)%*%(gphi))) # yr-1

f.det<-ben.Tempeffect*((1/10^x)*A.v*10^(x*alpha.v)*W[i])/(1+handling*(1/10^x)*A.v*10^(x*alpha.v)*W[i]) # yr-1 
 
# Predator growth integral 
GG.u[,i]<-(1-def.high)*K.u*(f.pel) + (1-def.low)*K.v*(f.ben)                #yr-1
      

# Reproduction
if (repro.on==1) R.u[,i]<-(1-def.high)*(1-(K.u+AM.u))*(f.pel) +  (1-def.high)*(1-(K.v+AM.v))*f.ben  #yr-1
                  
# Predator death integrals 

#Satiation level of predator for pelagic prey

sat.pel<-ifelse(f.pel>0,f.pel/((A.u*10^(x*alpha.u)*pref.pel)*(U[,i]*dx)%*%(gphi)),0)

PM.u[,i]<-as.vector((pref.pel*A.u*expax)*(U[,i]*sat.pel*dx)%*%(mphi))  #yr-1 

#PM.u[,i]<-as.vector((1-f.pel)*(A.u*10^(x*alpha.u)*pref.pel)*(U[,i]*dx)%*%(mphi))  #yr-1 


Z.u[,i]<- PM.u[,i]+ pel.Tempeffect*OM.u + SM.u + Fvec.u     #yr-1
	
# Benthos growth integral
GG.v[,i]<-(1-def.low)*K.d*f.det         #yr-1

#reproduction
if (repro.on==1) R.v[,i]<-(1-def.low)*(1-(K.d+AM.v))*(f.det) #yr-1
   
# Benthos death integral
#Satiation level of predator for benthic prey  

sat.ben<-ifelse(f.ben>0,f.ben/((A.u*10^(x*alpha.v)*pref.ben)*(V[,i]*dx)%*%(gphi)),0)

PM.v[,i]<-ifelse(sat.ben>0,as.vector((pref.ben*A.u*expax)*(U[,i]*sat.ben*dx)%*%(mphi)),0)  #yr-1

#PM.v[,i]<-as.vector((1-f.ben)*(A.u*10^(x*alpha.u)*pref.ben)*(U[,i]*dx)%*%(mphi))  #yr-1

Z.v[,i]<-PM.v[,i]+ ben.Tempeffect*OM.v + SM.v  + Fvec.u  #yr-1

#total biomass density eaten by pred (g.m-2.yr-1)

eatenbypred<-10^x*f.pel*U[,i] + 10^x*f.ben*U[,i] 

  
#detritus output (g.m-2.yr-1)

output.w<-sum(10^x*f.det*V[,i]*dx)   


#total biomass density defecated by pred (g.m-2.yr-1)
defbypred<-def.high*(f.pel)*10^x*U[,i]+ def.low*(f.ben)*10^x*U[,i]


#output fisheries catches per yr at size
Y.u[Fref.u:Nx,i]<-Fvec.u[Fref.u:Nx]*U[Fref.u:Nx,i]*10^x[Fref.u:Nx] 
#output fisheries catches per yr at size
Y.v[Fref.v:Nx,i]<-Fvec.v[Fref.v:Nx]*V[Fref.v:Nx,i]*10^x[Fref.v:Nx] 


#------------------------------------------------
# Increment values of W,U & V	for next time step  
#------------------------------------------------

#Detritus Biomass Density Pool - fluxes in and out (g.m-2.yr-1) of detritus pool and solve for detritus biomass density in next time step 

if (ERSEM.det.input==F) {
	
#considering pelagic faeces as input as well as dead bodies from both pelagic and benthic communities 
# and phytodetritus (dying sinking phytoplankton)


if (det.coupling==1.0) {
	
    # pelagic spectrum inputs (sinking dead bodies and faeces) - export ratio used for "sinking rate"
	#  + benthic spectrum inputs (dead stuff - already on/in seafloor)
	
	input.w<-(sinking.rate*(sum(defbypred[ref:Nx]*dx)
	+ sum(OM.u[1:Nx]*U[1:Nx,i]*10^(x[1:Nx])*dx) 
	+ sum(SM.u[1:Nx]*U[1:Nx,i]*10^(x[1:Nx])*dx))
	+ (sum(OM.v[1:Nx]*V[1:Nx,i]*10^(x[1:Nx])*dx) 
	+ sum(SM.v[1:Nx]*V[1:Nx,i]*10^(x[1:Nx])*dx)))
	
	}
	
	
if (det.coupling==0.0) {
	
	input.w<-sum(OM.v[1:Nx]*V[1:Nx,i]*10^(x[1:Nx])*dx) + sum(SM.v[1:Nx]*V[1:Nx,i]*10^(x[1:Nx])*dx)
	
	}

# change in detritus biomass density (g.m-2.yr-1)                     

# this one assumes additional losses to sediment?
# dW<-input.w - (output.w + W[i])     

# losses due to detritivory only:

dW<-input.w - (output.w)     

W[i+1]=W[i] + dW*delta_t        #biomass density of detritus g.m-2

}

if (ERSEM.det.input==T) {
W[i+1]=W[i]
}


#----------------------------------------------
# Pelagic Predator Density (nos.m-2)- solve for time + delta_t using implicit time Euler upwind finite difference (help from Ken Andersen and Richard Law)

# Matrix setup for implict differencing 
Ai.u<-Bi.u<-Si.u<-array(0, c(length(x), 1))   


Ai.u[idx]<-(1/log(10))*-GG.u[idx-1,i]*delta_t/dx
Bi.u[idx]<-1+(1/log(10))*GG.u[idx,i]*delta_t/dx +Z.u[idx,i]*delta_t
Si.u[idx]<-U[idx,i]

# Boundary condition at upstream end 
Ai.u[ref]<-0
Bi.u[ref]<-1
Si.u[ref]<-U[ref,i]

# Invert matrix

#recruitment at smallest consumer mass

#continuation of plankton hold constant  
U[1:ref,i+1]<-U[1:ref,i] 

# apply transfer efficency of 10% *plankton density at same size  
# if (repro.on==0)  U[ref,i+1]<-0.1*u.init.f(x,ui0,r.plank)[ref]       
# reproduction from energy allocation

if (repro.on==1)  U[ref,i+1]<-U[ref,i]+ (sum(R.u[(ref+1):Nx,i]*10^x[(ref+1):Nx]*U[(ref+1):Nx,i]*dx)*delta_t)/(dx*10^x[ref]) - (delta_t/dx)*(1/log(10))*(GG.u[ref,i])*U[ref,i] -delta_t*Z.u[ref,i]*U[ref,i]

#main loop calculation

for (j in (ref+1):(Nx)){
	
   #U[j,i+1]<-U[j,i]-(delta_t/dx)*(1/log(10))*(GG.u[j,i])*U[j,i] + (delta_t/dx)*(1/log(10))*GG.u[j-1,i]*U[j-1,i] -delta_t*Z.u[j,i]*U[j,i]
   
   U[j,i+1]<-(Si.u[j]-Ai.u[j]*U[j-1,i+1])/Bi.u[j]


}

#----------------------------------------------
# Benthic Detritivore Density (nos.m-2) 
Ai.v<-Bi.v<-Si.v<-array(0, c(length(x), 1))   
idx=(ref.det+1):Nx  #shorthand for matrix referencing

Ai.v[idx]<-(1/log(10))*-GG.v[idx-1,i]*delta_t/dx 
Bi.v[idx]<-1+(1/log(10))*GG.v[idx,i]*delta_t/dx +Z.v[idx,i]*delta_t
Si.v[idx]<-V[idx,i]

#boundary condition at upstream end
Ai.v[ref.det]<-0
Bi.v[ref.det]<-1
Si.v[ref.det]<-V[ref.det,i]  

#invert matrix
#recruitment at smallest detritivore mass  

#hold constant continution of plankton with sinking rate multiplier 
V[1:ref.det,i+1]<-V[1:ref.det,i]

# apply a very low of transfer effiency 1%* total biomass of detritus divided by minimum size
# if (repro.on==0)  V[ref.det,i+1]<-(0.01*W[i+1])/(10^x[ref]) 

if (repro.on==1)  V[ref.det,i+1]<-V[ref.det,i]+ sum(R.v[(ref.det+1):Nx,i]*10^x[(ref.det+1):Nx]*V[(ref.det+1):Nx,i]*dx)*delta_t/(dx*10^x[ref.det]) - (delta_t/dx)*(1/log(10))*(GG.v[ref.det,i])*V[ref.det,i] -delta_t*Z.v[ref.det,i]*V[ref.det,i]


#loop calculation
for (j in (ref.det+1):(Nx)){ 

   #V[j,i+1]<-V[j,i]-(delta_t/dx)*(1/log(10))*(GG.v[j,i])*V[j,i] + (delta_t/dx)*(1/log(10))*GG.v[j-1,i]*V[j-1,i]-delta_t*Z.v[j,i]*V[j,i] 

   V[j,i+1]<-(Si.v[j]-Ai.v[j]*V[j-1,i+1])/Bi.v[j]

}		
rm(j)




}#end time iteration


return(list(U=U[,],GG.u=GG.u[,],PM.u=PM.u[,],V=V[,],GG.v=GG.v[,],PM.v=PM.v[,],Y.u=Y.u[,],Y.v=Y.v[,],W=W[], params=params))


# return(list(U=U[,Neq+1],GG.u=GG.u[,Neq],PM.u=PM.u[,Neq],V=V[,Neq+1],GG.v=GG.v[,Neq],PM.v=PM.v[,Neq],Y=Y[,Neq],W=W[Neq+1], params=params))

})  
# end with(params)


}#end size-based model function


#---------------------------------------------------------------------------------------
# FUNCTION TO GET Parameters of model
#---------------------------------------------------------------------------------------


sizeparam<-function(equilibrium=F, dx=0.1,xmin=-12,xmax=6,xmin.consumer.u=-7,xmin.consumer.v=-7,tmax=10, delta_t=1/365,fmort.u =0.0,fminx.u=1, fmort.v = 0.0,fminx.v=1,depth=500,er=0.5,pp=-3,slope=-1){

#---------------------------------------------------------------------------------------
# Parameters of model
#---------------------------------------------------------------------------------------

param=list()

# plankton parameters   
param$ui0 <- 10^pp
param$r.plank= slope  # initial (fixed) slope of phyto-zooplankton spectrum

# fishing parameters 

param$Fmort.u=fmort.u      # fishing mortality rate for predators 
param$Fmort.v=fmort.v      # fishing mortality rate for detritivores
param$min.fishing.size.u= fminx.u # minimum log10 body size fished for predators
param$min.fishing.size.v= fminx.v # minimum log10 body size fished for detritivores
 
# benthic-pelagic coupling parameters

# set predator coupling to benthos, depth dependent - 0.75 above 500 m, 0.5 between 500-1800 and 0 below 	1800m (suggestions of values from Clive Trueman based on stable isotope work, and proportion of biomass, 	Rockall Trough studies)


  param$pref.ben = 0.8*exp(-1/250*depth)
  # if (depth <= 500) param$pref.ben= 0.75
  # if (depth > 500 && depth  <= 1800) param$pref.ben= 0.5
  # if (depth > 1800) param$pref.ben= 0.0

  param$pref.pel=1-param$pref.ben # preference for pelagic prey 
  param$det.coupling=1.0  # detritus coupling on? 1= yes, 0 = no
  param$sinking.rate=er # fraction of sinking detritus reaching the seafloor (from export ratio input)

  
# feeding and energy budget parameters

  param$q0=2.0      # Mean log10 predator prey mass ratio  100:1.
  param$sd.q=1.0    # 0.6 to 1.0 for lognormal prey preference function. 
  param$A.u=640     # hourly rate volume searched constant m3.yr-1 for fish. need to check this value, its quite large.
  param$A.v=param$A.u*0.1    # hourly rate volume filtered constant m3.yr-1 for benthos. this value yields believable growth curve.
              # approximatelty 10 times less than A.u
  param$alpha.u=0.82  #exponent for metabolic requirements plus swimming for predators(Ware et al 1978)
  param$alpha.v=0.75  #exponent =0.75 for whole organism basal (sedentary) metabolic rate (see growth.v) from Peters (1983) and Brown et al. (2004) for detritivores
  param$def.high=0.3  # fraction of ingested food that is defecated (Peters,1983)
  param$def.low=0.5   # low = low quality (K) food, high = high quality (K) food
  param$K.u=0.3       #net growth conversion efficiency for organisms in the "predator" spectrum form Ware (1978)
  param$K.v=0.2	      #net growth conversion efficiency for organisms in the "detritivore" spectrum
  param$K.d=param$K.v       #net growth conversion efficency for detritus
  param$AM.u=0.5      #fraction of energy required for maintenance & activity etc.
  param$AM.v=0.7
  param$handling=0    # 5.7e-7    # if handling time is > 0 Type II functional response, if = 0 linear (no predator satiation)
  param$repro.on=1  # dynamic reproduction on=1,off=0  
  param$c1=25.22 # constant used in Jennings et al. 2008 Proc B to standardize metabolism-temperature effects# for Boltzmann equation. Derived from Simon's fit to Andy Clarke's data
  param$E=0.63   # activation energy, eV
  param$k=8.62*10^-5 # Boltzmann's constant
  
   
# "other" mortality parameters

  param$mu0=0.2	      # residual natural mortality
  param$xs=3          # size at sensenscence
  param$p.s=0.3       # exponent for senescence mortality 
  param$k.sm =0.2     # constant for senescence mortality 
  
 
#---------------------------------------------------------------------------------------
# Parameters for numerical integration (size & time discretisation)
#---------------------------------------------------------------------------------------
param$dx     =  dx	  # size increment after discretization for integration (log body weight)
param$xmin   = xmin  # minimum log10 body size of plankton
param$x1     = xmin.consumer.u   # minimum log10 body size in dynamics predators
param$x1.det = xmin.consumer.v  # minimum log10 body size in dynamics benthic detritivores
param$xmax   =  xmax   # maximum log10 body size of predators
param$xmax2  =  4.0   # maximum log10 body size before senescence kicks in (departure form linearity)


#Vector with size bins
param$x    = seq(param$xmin, param$xmax, param$dx)
param$Nx  = length(param$x)

#position in x of x1
param$ref=((param$x1-param$xmin)/param$dx)+1

#position in x of x1.det
param$ref.det=((param$x1.det-param$xmin)/param$dx)+1 


#position in F vector corresponding to smallest size fished in U
param$Fref.u=((param$min.fishing.size.u-param$xmin)/param$dx)+1
#position in F vector corresponding to smallest size fished in V
param$Fref.v=((param$min.fishing.size.v-param$xmin)/param$dx)+1


#short hand for matrix indexing
param$idx=2:param$Nx


# integration parameters for equilibrium option

if (equilibrium==F) {
# only run for one step at a time (eg. for daily time-varying input)
param$tmax = 1/365 
# discretisation of year 
param$delta_t = 1/365
# number of time bins 
param$Neq=param$tmax/param$delta_t	
}


if (equilibrium==T){
# number of years to run model
param$tmax = tmax
# discretisation of year 
param$delta_t = delta_t
# number of time bins 
param$Neq=param$tmax/param$delta_t	

# # Set initial values for size spectra & detritus
	
# param$U.init=read.table("pp_0.006F_0pel_0.5ben_0.5sink_0.5_U.txt")
# param$V.init=read.table("pp_0.006F_0pel_0.5ben_0.5sink_0.5_V.txt")
# param$W.init=scan("pp_0.006F_0pel_0.5ben_0.5sink_0.5_Det.txt",quiet=TRUE)

  param$U.init<-param$ui0*10^(param$r.plank*param$x)  # (phyto+zoo)plankton + pelagic predator size spectrum  
  param$V.init<-param$sinking.rate*param$ui0*10^(param$r.plank*param$x)  # set initial detritivore spectrum  
  param$W.init<-0.00001  # abritrary initial value for detritus



}



return(param)
}#end sizeparam function



plotsizespectrum<-function(modeloutput,params, itime=1,pp=-2.2,timeaveraged=F) {

with(params, {	
#******************************PLOT**************************************
# plot changes in the two size spectra over time
ui0 <- 10^pp


U <- modeloutput[[1]][,itime]
V <- modeloutput[[4]][,itime]
W <- modeloutput[[7]][,itime]


if (timeaveraged==T) {
	
	U<-rowMeans(U)
	V<-rowMeans(V)
	W<-mean(W)

	}

maxy=max(log10(U))
miny=-20
   # par(mfrow=c(1,1))
   plot(x[ref:Nx],log10(U[ref:Nx]),type="l",col="blue",cex=1.6,ylab=c("log Abundance density [1/m3]"),xlab=c("log Body mass [g]"),xlim=c(x1.det,xmax),ylim=c(miny,maxy))
   points(x[ref.det:Nx],log10(V[ref.det:Nx]),type="l",col="red",cex=1.6,ylab=c(""),xlab=c(""))      
   # text(x1.det+abs(0.05*x1.det),miny+abs(1.3*miny),paste("day = ",itime,sep=""),pos=4)
   text(x1.det+abs(0.05*x1.det),miny+abs(0.3*miny),paste("pel.pref=",pref.pel,sep=""),pos=4)
   text(x1.det+abs(0.05*x1.det),miny+abs(0.2*miny),paste("ben.pref=",pref.ben,sep=""),pos=4)
   #text(x1.det+abs(0.05*x1.det),miny+abs(0.1*miny),paste("PP=",round(ui0,3),sep=""),pos=4)
   legend(xmax-0.48*(xmax-x1.det), maxy-0.05*maxy, c("Predators", "Detritivores"), col = c("blue","red"), lwd=2.0)

})
# end with(params) 

}	# end plot function





	
savetofile<-function(modeloutput,out_dir="C:/Ersem convert plankton/RUN MODEL MOST RECENT/",runname="test"){
   # write results of last timestep(?) to ascii files NEEDS TO BE STORED SO CAN BE READ BACK IN AS U.init,V.init and W.init.

   write.table(modeloutput[[1]],paste(out_dir,runname,"_U.txt",sep=""))
   write.table(modeloutput[[4]],paste(out_dir,runname,"_V.txt",sep=""))
   write.table(modeloutput[[7]],paste(out_dir,runname,"_Det.txt",sep=""))
 } # end of savetofile function
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

#-------------------------------------------------------------
#       Function for calculating right-hand side 
#
#       with Newton Raphson method to get steady -state
#-------------------------------------------------------------

Increment<-function(u_vals,params=params,pp=pp,bio.det= dat$BenC[ilatlon]*10,sst=dat$MLDT[ilatlon],sft=dat$BenT[ilatlon]){
with(params, {
	
ui0=10^pp        # intercept of plankton size spectrum estimated from ERSEM output.		
W.init=bio.det	 # detritus pool biomass density from ERSEM output
sst=sst          #temperature in water column
sft=sft		       #temperature near seabed
#---------------------------------------------------------------------------------------
# Initialising matrices
#---------------------------------------------------------------------------------------

#q1 is a square matrix holding the log(predatorsize/preysize) for all combinations of sizes
q1  = matrix(NA, length(x), length(y))
for (i in 1:length(y)) { q1[,i] = y[i] - x}

#q2 is the reverse matrix holding the log(preysize/predatorsize) for all combinations of sizes
q2 = matrix(-q1, length(x), length(y))	

#matrix for recording the two size spectra over time
V = U = array(0, c(length(x), 2))

#vector to hold detrtitus biomass density (g.m-2)
W = array(0,Neq)

#matrix for keeping track of growth and reproduction from ingested food:
R.v=R.u=GG.v = GG.u   = array(0, c(length(x), 2)) 

#matrix for keeping track of predation mortality
PM.v = PM.u   = array(0, c(length(x), 2))   

#matrix for keeping track of  total mortality (Z)
Z.v = Z.u = array(0, c(length(x), 2))

#matrix for keeping track of senescence mortality and other (intrinsic) mortality
SM.v = SM.u   =OM.v = OM.u   = array(0, length(x))

#empty vector to hold fishing mortality rates at each size class
Fvec = rep(0,length(x))

#lookup tables for terms in the integrals which remain constant over time
gphi  = gphi.f(q1)
mphi  = mphi.f(q2)

#lookup table for components of 10^(alpha.u*x)
expax = expax.f(x, alpha.u)


#------------------------
# Initial values
#------------------------


U[1:(ref-1),1]<-u.init.f(x,ui0,r.plank)[1:(ref-1)]
U[ref:Nx,1]<-u_vals[1:(Nx-ref+1)]
V[ref.det:Nx,1]<-u_vals[(Nx-ref+2):(Nx-ref.det+Nx-ref+2)]
W[1]<-W.init

# Predator Feeding integral
                    
#intrinsic natural mortality
OM.u<-mu0*10^(-0.25*x)
OM.v<-mu0*10^(-0.25*x)

#senescence mortality rate to limit large fish from building up in the system
#same function as in Law et al 2008, with chosen parameters gives similar M2 values as in Hall et al. 2006
SM.u=k.sm*10^(p.s*(x-xs))
SM.v=k.sm*10^(p.s*(x-xs))

#fishing mortality
Fvec[Fref:Nx] = 0.09*(x[Fref:Nx]/ log10(exp(1)) ) + 0.04  





#--------------------------------
# Calculate Growth and Mortality
#--------------------------------
             
                        
pel.Tempeffect=exp(c1 - E/(k*(sst+ 273)))
ben.Tempeffect=exp(c1 - E/(k*(sft+ 273)))


# feeding rates
f.pel<-pel.Tempeffect*as.vector(((A.u*10^(x*alpha.u)*pref.pel)*(U[,1]*dx)%*%(gphi))/(1+handling*(A.u*10^(x*alpha.u)*pref.pel)*(U[,1]*dx)%*%(gphi))) # yr-1
f.ben<-pel.Tempeffect*as.vector(((A.u*10^(x*alpha.u)*pref.ben)*(V[,1]*dx)%*%(gphi))/(1+handling*(A.u*10^(x*alpha.u)*pref.ben)*(V[,1]*dx)%*%(gphi))) # yr-1
f.det<-ben.Tempeffect*((1/10^x)*A.v*10^(x*alpha.v)*W[1])/(1+handling*(1/10^x)*A.v*10^(x*alpha.v)*W[1]) # yr-1 
 
# Predator growth integral 
GG.u[,1]<-(1-def.high)*K.u*(f.pel) + (1-def.low)*K.v*(f.ben)                #yr-1
   
# Reproduction
if (repro.on==1) R.u[,1]<-(1-def.high)*(1-(K.u+AM.u))*(f.pel) +  (1-def.high)*(1-(K.v+AM.v))*f.ben  #yr-1
                  
# Predator death integrals 

#Satiation level of predator for pelagic prey
sat.pel<-ifelse(f.pel>0,f.pel/((A.u*10^(x*alpha.u)*pref.pel)*(U[,1]*dx)%*%(gphi)),0)
PM.u[,1]<-as.vector((pref.pel*A.u*expax)*(U[,1]*sat.pel*dx)%*%(mphi))  #yr-1 

Z.u[,1]<- PM.u[,1]+ pel.Tempeffect*OM.u + SM.u + Fvec*Fmult     #yr-1
	
# Benthos growth integral
GG.v[,1]<-(1-def.low)*K.d*f.det         #yr-1

#reproduction
if (repro.on==1) R.v[,1]<-(1-def.low)*(1-(K.d+AM.v))*(f.det) #yr-1
   
# Benthos death integral
#Satiation level of predator for benthic prey  
sat.ben<-ifelse(f.ben>0,f.ben/((A.u*10^(x*alpha.v)*pref.ben)*(V[,1]*dx)%*%(gphi)),0)
PM.v[,1]<-ifelse(sat.ben>0,as.vector((pref.ben*A.u*expax)*(U[,1]*sat.ben*dx)%*%(mphi)),0)  #yr-1

Z.v[,1]<-PM.v[,1]+ ben.Tempeffect*OM.v + SM.v    #yr-1

#total biomass density eaten by pred (g.m-2.yr-1)
eatenbypred<-10^x*f.pel*U[,1] + 10^x*f.ben*U[,1] 

#detritus output (g.m-2.yr-1)
output.w<-sum(10^x*f.det*V[,1]*dx)   

#total biomass density defecated by pred (g.m-2.yr-1)
defbypred<-def.high*f.pel*10^x*U[,1]+ def.low*f.ben*10^x*U[,1]

#------------------------------------------------
# Increment values of W,U & V	for next time step  
#------------------------------------------------

#Detritus Biomass Density Pool - fluxes in and out (g.m-2.yr-1) of detritus pool and solve for detritus biomass density in next time step 
#if (ERSEM.det.input==F) {
##considering pelagic faeces as input as well as dead bodies from both pelagic and benthic communities 
## and phytodetritus (dying sinking phytoplankton)
#if (det.coupling==1.0) input.w<-sinking.rate*(sum(defbypred[ref:Nx]*dx)+ sum(OM.u[1:Nx]*U[1:Nx,1]*10^(x[1:Nx])*dx) + sum(SM.u[1:Nx]*U[1:Nx,1]*10^(x[1:Nx])*dx)) + sum(OM.v[1:Nx]*V[1:Nx,1]*10^(x[1:Nx])*dx) + sum(SM.v[1:Nx]*V[1:Nx,1]*10^(x[1:Nx])*dx)
#if (det.coupling==0.0) input.w<-sum(OM.v[1:Nx]*V[1:Nx,1]*10^(x[1:Nx])*dx) + sum(SM.v[1:Nx]*V[1:Nx,1]*10^(x[1:Nx])*dx)
#
## change in detritus biomass density (g.m-2.yr-1)                     
#dW<-input.w - (output.w + W[1])     
#W[i+1]=W[1] + dW*delta_t        #biomass density of detritus g.m-2
#
#}
#
#if (ERSEM.det.input==T) {
             
W[i+1]=W[1]   #hold constant as ERSEM input
#}



#hold plankton constant thru time on left hand boundary
U[1:(ref-1),2]<-u.init.f(x,ui0,r.plank)[1:(ref-1)]

#recruitment at smallest consumer mass  1.-hold constant OR 2.- scale as function of food consumed
# 3. combination of background plankton + eggs produced

#old one
#U[ref,2]<-((sum(gsi*eatenbypred[ref:Nx]*dx) + u.init.f(x,ui0,r.plank)[ref]) - U[ref,1])/delta_t

#new renewal as a rate

U[ref,2]<-((U[ref,1]
  + (sum(R.u[(ref+1):Nx,1]*10^x[(ref+1):Nx]*U[(ref+1):Nx,1]*dx)*delta_t)/(dx*10^x[ref])
  - (delta_t/dx)*(1/log(10))*(GG.u[ref,1])*U[ref,1] 
  -delta_t*Z.u[ref,1]*U[ref,1])
  - U[ref,1])/delta_t
 
 
#loop calculation
for (j in (ref+1):Nx){

U[j,2]<--(Z.u[j,1]*U[j,1])-(1/log(10))*((GG.u[j,1]*U[j,1]-GG.u[j-1,1]*U[j-1,1])/dx)
}

  
  #Option 1:hold micro and meiofauna constant on left hand boundary

  #V[1:(ref-1),i+1]=V[1:(ref-1),1]

  # Option2: start benthos at smallest size (or could start at ref +1 (x = x1)
        #or could let recruitment of detritivores determined by food density
 #old one: 
 #V[ref.det,2]<- (((1/K.d*10^x[ref.det])*sum(gsi*10^x[(ref.det+1):Nx]*GG.v[(ref.det+1):Nx,1]*V[(ref.det+1):Nx,1]*dx) + u.init.f(x,vi0,r.v)[ref.det])- V[ref.det,1])/delta_t
 
 #new renewal as a rate
V[ref.det,2]<-((V[ref.det,1]
+ sum(R.v[(ref.det+1):Nx,1]*10^x[(ref.det+1):Nx]*V[(ref.det+1):Nx,1]*dx)*delta_t/(dx*10^x[ref.det]) 
- (delta_t/dx)*(1/log(10))*(GG.v[ref.det,1])*V[ref.det,1] 
- delta_t*Z.v[ref.det,1]*V[ref.det,1])
-V[ref.det,1])/delta_t

  #loop calculation

  for (j in (ref.det+1):Nx){ 
    V[j,2]<--(Z.v[j,1]*V[j,1])-(1/log(10))*((GG.v[j,1]*V[j,1]-GG.v[j-1,1]*V[j-1,1])/dx)
  }		
	

return(c(U[ref:Nx,2],V[ref.det:Nx,2]))

} )
#end with
 
}#end Increment function








#----------------------------------------------------
#
#
#     ers2siz
#
#     FUNCTION TO CONVERT ERSEM OUTPUT INTO PHYTOPLANKTON ISZE SPECTRUM INTERCEPT
#
#
#     this needs to be called before any of the above
#
#-------------------------------------------------------


ers2siz<-function(Diatoms,Flagellates, Picophyto,nanoflagellates, microzooplankton,plot_results=F,zooponly=T,MLD){
	
	
# Convert plankton functional group biomass density estimates (g C m^-3) from GETM-BFM output
# into abundance-body mass distributions for size spectrum model input

# NOTE : I think we need to use the total biomass density integrated in the water column  (g C m^-2) and then divide by depth to get (g C m^-3)
# Otherwise, won't we underestimate the amount of plankton biomass?   Maybe compare both for the time being.


#Phytoplankton gC m^-2 (if depth intergrated, then make sure divide by mixed layer depth later, outside of function)
#Diatoms<-8          #2-200 um  P1 ERSEM guide ESD 20-200 mu m
#Flagellates<-8      #2-200 um  P2 ERSEM guide (nanoflagellates) ESD 2-20 mu m 
#Picophyto<-10       #0.2-2 um  P3 ERSEM guide ESD 0.2-2 mu m
                    # Inedibles (e.g. dinoflagellates) ERSEM guide ESD > 100 mu m

#Nanomicrophyto<-Diatoms+Flagellates   # P1+P2

#Zooplankton gC m^-2:
#nanoflagellates<-15      #2-20 um Z6
#microzooplankton<-15     #20-200 um Z5

#INPUT SIZE RANGE = 0.2 to 200 um  (pico,nano,micro phyto and zooplankton, exclude meso zooplankton (200-2000 um) )

#######################################################################  
#convert gC M^-2 to biomass (grams wet weight M^-2)
#from coversions in Boudreau & Dickie 1991
#biomass
#1 g C = 10kcal = 10 g
# 1 mgC/m3=1E-3 gC/m3 = 10*1E-3 kcal/m3=1E-2 g/m3
fbiomass<-function(val){val=val*10}

#######################################################################
#size
#convert um to grams wet weight from Boudreau & Dickie
#convert diameter in um to radius in centimetres and then to volume (cubic centimetre)
#assume 1 cc= 1 mL =1 g = 1 kcal
fsize<-function(size){(4.0/3.0)*pi*((0.5*0.0001*size)^3)}

#######################################################################


#set up size ranges min,max converting from um to grams wet weight and create bins of log10 grams=0.1
# --> create volume weigth ranges based on size ranges
picosize<-seq(round(log10(fsize(c(0.2))),1),round(log10(fsize(c(2))),1)-0.1,0.1)
nanosize<-seq(round(log10(fsize(c(2))),1),round(log10(fsize(c(20))),1)-0.1,0.1)
microsize<-seq(round(log10(fsize(c(20))),1),round(log10(fsize(c(200))),1)-0.1,0.1)
nanomicrosize<-seq(round(log10(fsize(c(2))),1),round(log10(fsize(c(200))),1)-0.1,0.1)
#fullsize<-c(picosize,nanosize,microsize)
fullsize<-c(nanosize,microsize)
fullN<-rep(0,length(fullsize)) # create zero vector of equal length as fullsize

sizespec<-function(totalbiomass,sizeclasses,dx=0.1,r=1){rep(totalbiomass/(dx*length(sizeclasses)),length(sizeclasses))*10^(-r*sizeclasses)}
#assumes equal biomass across size classes, B=M^0  and numerical density scales as N=M^-1
#holds if including zooplankton or across more than one trophic level
#note that Simon (Proceedings paper) assumed scaling of B=M^1/4 and N=M^-3/4 for primary producers only

#returns numbers.m^-2

##phyto<-c(sizespec(fbiomass(Picophyto),fullsize[fullsize==picosize],r=1),sizespec(fbiomass(Nanomicrophyto),fullsize[fullsize>-11.5],r=1))



phyto<-c(sizespec(fbiomass(Nanomicrophyto),fullsize[fullsize>-11.5],r=1))


zoo<-c(sizespec(fbiomass(nanoflagellates),fullsize[fullsize==nanosize],r=1),sizespec(fbiomass(microzooplankton),fullsize[fullsize==microsize],r=1))


#fullN[fullsize==picosize]<-phyto[fullsize==picosize]
fullN[fullsize==nanosize]<-phyto[fullsize==nanosize] + zoo[1:length(nanosize)]
fullN[fullsize==microsize]<-phyto[fullsize==microsize]+zoo[(length(nanosize)+1):length(zoo)]

if (plot_results) {
 # plot(fullsize,log10(phyto),col=1)                  #phytoplankton size spectrum
 # points(fullsize[fullsize>-11.5],log10(zoo),col=2)  #zooplankton size spectrum
  plot(fullsize,log10(fullN),col=4,type="p")                #phyto + zoo whole plankton size spectrum, bg=4,pch=21

  #maxy=max( max(log10(phyto)),max(log10(zoo)) )
  #legend(max(fullsize)-2.5, floor(maxy), c("phyto", "zoo", "full"), col = c(1,2,4), lwd=2.0)
  #rm(maxy)
}

#check slope & intercept
#lm(log10(fullN[fullN>0])~fullsize[fullN>0],na.action=na.exclude)





#output intercept to be used as input for size spectrum model - slope coefficient not used.

if (sum(fullN)==0) pp<-NA

if (sum(fullN)>0) {

     	if(zooponly==F) pp<-coef(lm(log10(fullN[fullN>0])~fullsize[fullN>0],na.action=na.exclude))[1]

		if(zooponly==T) pp<-coef(lm(log10(zoo[zoo>0])~fullsize[zoo>0],na.action=na.exclude))[1]

}

# pp<-cbind(fullsize,fullN) #could instead just output estimated densities-at-size (not used)



return(pp)

}#end ers2siz function



#Method to calculate intercept from appendix of Barnes et al. 2010 JPR


GetPPIntSlope<-function(Picophyto,Nanomicrophyto,mmin=10^-14.25, mmid=10^-10.184,mmax=10^-5.25) {
	
#mmin=10^-14.25 gww, approx 0.2 ESD
#mmid=10^-10.184 gww, approax 5 ESD - division between small and large phytoplankton in GFDL-Topaz
#mmax=10^-5.25 gww, approx 200 ESD

	
#from Appendix of Barnes 2010 JPR
	
#the scaling of biomass with body mass can be described as B=aM^b

#the exponent b (also equivalent to the slope b in a log B vs log M relationship) can be assumed:
# 0.25 (results in N ~ M^-3/4) or 0 (results in N ~ M^-1)
# most studies seem to suggest N~M^-1, so can assume that and test sensitivity of our results to this assumption. 

#b=0

# or can be calculated from the difference in biomass of pico and nanophytoplankton at thier assumed midpoint sizes

###------Approach One------###
#Original

#midpico<-(0.2+5)/2
#midnano<-(5+200)/2
#
#b<-(Picophyto-Nanomicrophyto)/(midpico-midnano)
#
##if totalbiomass is known within each group, with a known ( or assumed ) size range and a known within-group b value, a can be calculated as follows:
#
#totalbiomass=(Picophyto+Nanomicrophyto)*10 # first convert g C to gww
#
#a=(b+1)*totalbiomass/(mmax^(b+1) - mmin^(b+1)) # and use body mass range in gww
#nb<-b-1 #Exponent/slope for numerical density versus M
#
#pp<-c(log10(a),nb) # log10(a) will equal to the log10 density of particles at mass = 1 g and at log10 (mass)=0
#
###------Approach Two------###
#Change to calculate a and b in log B (log10 abundance) vs. log M (log10 gww)

midpico<-log10((mmin+mmid)/2) #in log10 gww
midnano<-log10((mmid+mmax)/2) #in log10 gww

pico<-log10((Picophyto*10)/10^midpico)       #convert to log10 (gww/size class median size) for log10 abundance
nano<-log10((Nanomicrophyto*10)/10^midnano)  #convert to log10 (gww/size class median size) for log10 abundance

b<-(pico-nano)/(midpico-midnano)

a = nano - b*midnano  #a is really log10(a), same a when pico, midpico are used
pp<-c(a,b) # log10(a) will equal to the log10 density of particles at mass = 1 g and at log10 (mass)=0
			 # a could be used directly to replace 10^pp in sizemodel()
}


# Added 1 May, 2013
# Methods to calculate:

# 1. primary producer size distribution (as described by MB10, 50, 90)
# 2. primary prodcuer slope and intercept (from Barnes et al. and Simon's notes)
# 3.detritus input dependent on primary production and depth - proportion of sinking phytoplankton for detritus input
# 4. minimum consumer size
# 5. input for preferred predator-prey mass ratio function - beta and sigma

# These are calculated from measurements of temperature, primary production 

GetPPSpectrum<-function(pp = 0.77, sst = 20, chla=NA) {
	
# List of parameters
# sst= temperature (celsius) - INPUT, in degrees C
# P= primary production Note: needs to be in g C m-2 d-1  - INPUT
# chla = chl a biomass Note: mg Chl m-3, if used instead of pp input !=NA


# OUTPUTS:
# a= intercept phytoplankton size spectrum (orginal units are (x) cell mass pg C and (y) log10 pg C m-3 converted to g vs. number density per m^3)
# b= slope phytoplankton size spectrum 
# mb50= as Barnes et al

# /* phytoplankton size spectrum parameters */

 log10pp=log10(pp*1000)

# predict log 10 (pg C) median phytoplankton cell size

 mb50log=(-0.081*sst)+(0.875*log10pp)-0.661

# predict phytoplankton ss intercept from PP as g C m-2 d-1
# spectrum units are log10 pg C m-3 (y) v cell mass log10 pg C (x)

 mbint=(0.554*log10pp)+8.087

# predict phytoplankton ss slope from PP as g C m-2 d-1 
# spectrum units are log10 pg C m-3 (y) v cell mass log10 pg C (x)

 mbslope=(0.182*log10pp)-1.715
 
 
# fixed log predator=prey mass ratio to calculate minimum consumer size in the model

beta = 2.0
 
# need to convert log 10 pg C to g C (10^-12) then g C to g wet weight (10^-1; 1 g C = 10 g ww)

mb50a = 10^(mb50log)*10^-12*10^-1

# mb50 now in log 10 (grams wet weight)

mb50 = log10(mb50a)

return(list(

# OUTPUTS: 
# a is intercept - log 10 (pg C m^-3) - convert to log10 (g ww m^3) -> log10 (density m^3) 
a=log10((10^mbint)*10^-13),
# pp slope
b = mbslope,
# log 10 (median pp size, g ww),
med.x = mb50,
# estimate min consumer size from PPMR and median phytoplankton size ( see also Woodworth et al 2012)
min.cons.x = log10(10^beta*10^mb50)

))
 
}




#------------------------------------------------------------------------------------- 
#
# Function to Optimise yield , based on a single parameter, (F-even, a single F across all sizes) 
# 
#-------------------------------------------------------------------------------------


## ------------------------------------------------------------------
##  Do a single-parameter optimization of  yield
## ------------------------------------------------------------------


optimizeYield <- function(params,com,pp=pp, bio.det= W, sst=sst, sft=sft, U_mat=U, V_mat=V,sizeindex=(1:params$Nx)) {

  
  # Helper function to return the quantity to be optimized based on the
  # fishing mortality :
  #
    
    
    
  totalYield <- function(newF,params, com,sizeindex) {
    
    print(newF)
            
    params$Fmort = newF
    
    
    params$tstepdays = 1.0
    params$tmaxyears = 10   #change to 10 years( not running all the way to equilib within this)
    params$Neq=as.integer(365.0*params$tmaxyears/params$tstepdays)	# total number of steps in a period
    params$delta_t=params$tmaxyears/params$Neq               		# time step (rates per year)  
    params$siz_time=1:params$Neq

    
    params$U.init[,1]=com$U
    params$V.init[,1]=com$V
    params$W.init[1]=com$W 
    
    com <- sizemodel(params=params, pp, bio.det, sst, sft, U_mat, V_mat, temp.effect=T) 
           
           
    # mean yield through time for specific size range ( or entire fished range)    
    
    yield = com$Y[sizeindex]
 

    totalYield <- sum(yield)*1e-6
           
    
    -totalYield
  }
 
    
  result <- optimize(f=totalYield,interval=c(0,1.5),maximum=F,params=params, com=com)
  
  # result <- nlm(f=totalYield,p=0.8,params=params, com=com)
  
  # result <- genoud(fn=totalYield, nvars=1,Domains=cbind(0,2), max=FALSE,params=params, com=com)
  
  # result<-GenSA(par=0.1, lower=0, upper=3, fn=totalYield,params=params,com=com,control=list(maxit=100))
  
  return(list(result,com))
  
}


