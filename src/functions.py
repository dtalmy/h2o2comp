#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 13:57:23 2024

@author: dkm
"""

# Base parameter set
# units: N in mM, cells in cells mL-1, H2O2 in pmol mL-1
kss = 0.17*(4/3*3.142*0.5**3)*0.27 * 1e-3 # mM
ksp = 0.17*(4/3*3.142*0.3**3)*0.27 * 1e-3 # mM
mumaxp =  0.63    # day-1, moore and chisholm 1999
mumaxs =  0.75    # day-1, moore, goericke, chisholm, 1995
dp = 0.1   # pro mortality, day-1
ds =  0.1  # syn mortality, day-1
kdam = 0.002      # mL pmol-1 day-1
deltah = 0.001    # background H2O2 decay, day-1, mccullough 2025
phi = 1.7e-6      # syn detox, mL cell-1 day-1, mccullough 2025
rho =  0.1        # ammonium loss rate, day-1
Qnp = 9.4e-15/14.0*1e+6  # Pro N quota, umol N cell-1, Bertilsson
Qns = 20.0e-15/14.0*1e+6 # Syn N quota, umol N cell-1

# tendencies
def leak(y,t,params):
    ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55,flag = params
    P,S,N,H = max(y[0],1e-30),max(y[1],1e-30),max(y[2],1e-30),max(y[3],1e-30)
    dPdt = (mumaxp * N /( (ksp) + N) )*P - (dp *P) - kdam*H*P
    dSdt =(mumaxs * N /( (kss) + N))*S - (ds *S)     
    dNdt =  SN - Qnp*(mumaxp * N /(ksp + N)*P) - Qns*(mumaxs * N /(kss + N)*S) - rho*N    
    dHdt = Sh - deltah*H  - phi*S*H - 13e-6*EZ55*H
    if flag == True:
        print(P,S,N,H)
        print(dPdt,dSdt,dNdt,dHdt)
    return [dPdt,dSdt,dNdt,dHdt]

#################################################
# equilibrium solutions
################################################

def death_all(params):
    ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55,flag = params
    Nstar = SN/(rho)
    Pstar = 0
    Sstar = 0
    Hstar = Sh/(deltah+13e-6*EZ55)
    return  Nstar, Pstar, Sstar, Hstar 

def Pwins (params): 
    ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55,flag = params
    Hstar = Sh/(deltah+13e-6*EZ55)
    Nstar =  ((ds +kdam*Hstar)*ksp)/(mumaxp - ((ds+kdam*Hstar)))
    Pstar = ((SN-rho*Nstar)*(Nstar+ksp)) / (mumaxp*Nstar*Qnp)
    Sstar = 0
    return  Nstar, Pstar, Sstar, Hstar 

def Swins (params): 
    ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55,flag = params
    Nstar = (ds*kss)/(mumaxs-ds)
    Pstar = 0
    Sstar = ((SN-rho*Nstar)*(Nstar+kss)) / (mumaxs*Nstar*Qns)
    Hstar = Sh/(phi*Sstar + deltah+13e-6*EZ55)
    return  Nstar, Pstar, Sstar, Hstar 

def Coexist (params): 
    ksp,kss,mumaxp,mumaxs,dp,ds,kdam,deltah,phi,rho,SN,Sh,Qnp,Qns,EZ55,flag = params
    Nstar = (ds*kss)/(mumaxs-ds)
    Hstar = 1/kdam*(mumaxp*Nstar/(ksp+Nstar)-ds)
    Sstar = (Sh/Hstar-deltah-13e-6*EZ55) /(phi)
    Pstar = (SN - ((mumaxs*Nstar*Qns)/(Nstar + (kss)))*Sstar  - rho*Nstar)  /((Qnp*mumaxp*Nstar)/ ((Nstar +( ksp))))
    return  Nstar, Pstar, Sstar, Hstar 

