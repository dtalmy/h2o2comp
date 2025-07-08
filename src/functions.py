#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 13:57:23 2024

@author: dkm
"""

#Base parameter set  (leaky or no H are special cases) 

# allomtry for half saturation constants
kss = 0.17*(4/3*3.142*0.5**3)*0.27 # (mumol N l-1) 
ksp = 0.17*(4/3*3.142*0.3**3)*0.27 # (mumol N l-1) 

k1p =  0.02*100000     #Pro alpha
k1s =  0.01*100000    #Syn alpha 
k2p =  0.63    # moore and chisholm 1999
k2s =  0.75    # moore, goericke, chisholm, 1995
dp = 0.1   #pro delta
ds =  0.1   #syn delta
kdam = 0.005555   #hooh mediated damage rate of Pro  
deltah = 0.001       # background decay rate
phi = 1.7e-6    #0007  #detoxification-based decay of HOOH via Syn in this case
rho =  0.2
Sh = 400
SN = 8000
Qnp = 9.4e-15/14.0*1e+6  #Nitrogen Quota for Pro from Bertilison in mumol N cell-1
Qns = 20.0e-15/14.0*1e+6 # and for syn 

params = [ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh]

##############
#make param for heterotroph detox value and heterotroph to be unchanging with current N and also to be callable into leaky function
#############

#functions for model 

def leak(y,t,params):
    ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55 = params
    P,S,N,H = max(y[0],1e-30),max(y[1],1e-30),max(y[2],1e-30),max(y[3],1e-30)
    dPdt = (k2p * N /( (ksp) + N) )*P - (dp *P) - kdam*H*P
    dSdt =(k2s * N /( (kss) + N))*S - (ds *S)     
    dNdt =  SN - Qnp*(k2p * N /(ksp + N)*P) - Qns*(k2s * N /(kss + N)*S) - rho*N    
    dHdt = Sh - deltah*H  - phi*S*H - 13e-6*EZ55*H
    return [dPdt,dSdt,dNdt,dHdt]

#################################################
#equilibrium solutions for leaky 
################################################

def death_all(params):
    ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55 = params
    Nstar = SN/(rho)
    Pstar = 0
    Sstar = 0
    Hstar = Sh/(deltah+13e-6*EZ55)
    return  Nstar, Pstar, Sstar, Hstar 

def Pwins (params): 
    ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55 = params
    Hstar = Sh/(deltah+13e-6*EZ55)
    Nstar =  ((ds +kdam*Hstar)*ksp)/(k2p - ((ds+kdam*Hstar)))
    Pstar = ((SN-rho*Nstar)*(Nstar+ksp)) / (k2p*Nstar*Qnp)
    Sstar = 0
    return  Nstar, Pstar, Sstar, Hstar 

def Swins (params): 
    ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55 = params
    Nstar = (ds*kss)/(k2s-ds)
    Pstar = 0
    Sstar = ((SN-rho*Nstar)*(Nstar+kss)) / (k2s*Nstar*Qns)
    Hstar = Sh/(phi*Sstar + deltah+13e-6*EZ55)
    return  Nstar, Pstar, Sstar, Hstar 

def Coexist (params): 
    ksp,kss,k2p,k2s,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55 = params
    Nstar = (ds*kss)/(k2s-ds)
    Hstar = 1/kdam*(k2p*Nstar/(ksp+Nstar)-ds)
    Sstar = (Sh/Hstar-deltah-13e-6*EZ55) /(phi)
    Pstar = (SN - ((k2s*Nstar*Qns)/(Nstar + (kss)))*Sstar  - rho*Nstar)  /((Qnp*k2p*Nstar)/ ((Nstar +( ksp))))
    return  Nstar, Pstar, Sstar, Hstar 

