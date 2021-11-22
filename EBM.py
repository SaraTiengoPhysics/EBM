#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import math

cp = 1.05*1e23/(4*np.pi*(6370*1e3)**2) # global heat capacity averaged and divided by Earth surface, in [J/K/m2]
dt = 60*60*24*100 # time step in s  [100 days] 
So = 1340.0       # W/m2, default solar constant
S2 = -0.477     # constant in polynomial expr of solar lat-distribution
alpha_ice = 0.62
alpha_no_ice = 0.3
Tice = 258. #[K]
Tnoice = 288. #[K]
A = -367.3 
B = 2.09


def SW_flux(Ts):
    #alb = albedo(Ts)
    S = (So/4) * (1-albedo(Ts))
    return S
    
def LW_flux(Ts):
    F = A+B*Ts
    return F

def albedo(Ts):
    Xs = 100
    if Ts > Tnoice:
        Xs = 1
    elif Ts < Tice:
        Xs = 0
    else:
        Xs = 1+((Ts-Tnoice)/(Tnoice-Tice))
        
    albedo = alpha_ice + (alpha_no_ice - alpha_ice)*((1-S2/2)*Xs+S2/2*Xs**3)
    return albedo

def find_zero_crossing(y,x):
    i = 0
    solutions = np.zeros(10)
    #print(y[0])
    dx = y[0]
    for j in np.arange(1,len(y)):
        dx1 = y[j]
        if(dx1*dx>0):
            dx = dx1
        else:
            solutions[i] = np.interp((abs(dx)),(0,np.abs(dx1-dx)),x[j-1:j+1])
            i = i+1
            dx = dx1
    return solutions[0:i]
            
def runEBM0d(Tso,nt):
    #Ts0: initial condition
    #nt: number of steps that I want to calculate
    TsSerie = np.zeros(int(nt))
    Ti = Tso #initial value of T
    for i in range(nt):
        
        dT = ((So/4*(1-albedo(Ti))-(A+B*Ti))*(dt/cp)) #devo cambiare alpha
        Tf = Ti + dT #value of T after a period dt 
        TsSerie[i] = Tf
        Ti = Tf
    return TsSerie
