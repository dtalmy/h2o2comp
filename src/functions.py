#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 13:57:23 2024

@author: dkm
"""

#Base parameter set
kss = 0.17*(4/3*3.142*0.5**3)*0.27 # (mumol N l-1) 
ksp = 0.17*(4/3*3.142*0.3**3)*0.27 # (mumol N l-1) 
mumaxp =  0.63    # moore and chisholm 1999
mumaxs =  0.75    # moore, goericke, chisholm, 1995
dp = 0.1   #pro delta
ds =  0.1   #syn delta
kdam = 0.005555   #hooh mediated damage rate of Pro mccullough 2025 
deltah = 0.001       # background decay rate mccullough 2025
phi = 1.7e-6    #0007  #syn detox decay of HOOH mccullough 2025
rho =  0.2 # assumed
Qnp = 9.4e-15/14.0*1e+6  #Nitrogen Quota for Pro from Bertilison in mumol N cell-1
Qns = 20.0e-15/14.0*1e+6 # and for syn 

##############
#make param for heterotroph detox value and heterotroph to be unchanging with current N and also to be callable into leaky function
#############

# tendencies
def leak(y,t,params):
    ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55 = params
    P,S,N,H = max(y[0],1e-30),max(y[1],1e-30),max(y[2],1e-30),max(y[3],1e-30)
    dPdt = (mumaxp * N /( (ksp) + N) )*P - (dp *P) - kdam*H*P
    dSdt =(mumaxs * N /( (kss) + N))*S - (ds *S)     
    dNdt =  SN - Qnp*(mumaxp * N /(ksp + N)*P) - Qns*(mumaxs * N /(kss + N)*S) - rho*N    
    dHdt = Sh - deltah*H  - phi*S*H - 13e-6*EZ55*H
    return [dPdt,dSdt,dNdt,dHdt]

#################################################
# equilibrium solutions
################################################

def death_all(params):
    ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55 = params
    Nstar = SN/(rho)
    Pstar = 0
    Sstar = 0
    Hstar = Sh/(deltah+13e-6*EZ55)
    return  Nstar, Pstar, Sstar, Hstar 

def Pwins (params): 
    ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55 = params
    Hstar = Sh/(deltah+13e-6*EZ55)
    Nstar =  ((ds +kdam*Hstar)*ksp)/(mumaxp - ((ds+kdam*Hstar)))
    Pstar = ((SN-rho*Nstar)*(Nstar+ksp)) / (mumaxp*Nstar*Qnp)
    Sstar = 0
    return  Nstar, Pstar, Sstar, Hstar 

def Swins (params): 
    ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55 = params
    Nstar = (ds*kss)/(mumaxs-ds)
    Pstar = 0
    Sstar = ((SN-rho*Nstar)*(Nstar+kss)) / (mumaxs*Nstar*Qns)
    Hstar = Sh/(phi*Sstar + deltah+13e-6*EZ55)
    return  Nstar, Pstar, Sstar, Hstar 

def Coexist (params): 
    ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55 = params
    Nstar = (ds*kss)/(mumaxs-ds)
    Hstar = 1/kdam*(mumaxp*Nstar/(ksp+Nstar)-ds)
    Sstar = (Sh/Hstar-deltah-13e-6*EZ55) /(phi)
    Pstar = (SN - ((mumaxs*Nstar*Qns)/(Nstar + (kss)))*Sstar  - rho*Nstar)  /((Qnp*mumaxp*Nstar)/ ((Nstar +( ksp))))
    return  Nstar, Pstar, Sstar, Hstar 

