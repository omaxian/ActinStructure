#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 14:35:34 2026

@author: sophiachainani
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def ode(t, y, k1, k2):
    A = y[0]
    B = y[1]
    AB = y[2]
    
    dA = -k1*A*B + k2*AB
    dB = -k1*A*B + k2 * AB 
    dAB = k1*A*B - k2 * AB
    return [dA, dB, dAB]

def microscopic(r1, r2, NA, NB, NAB, tmax):
    t = 0
    times = [t]
    A = [NA]
    B = [NB]
    AB = [NAB]
    while t<tmax:
        T1 = -np.log(np.random.rand())/(r1*NA*NB)
        T2 = -np.log(np.random.rand())/(r2*NAB)
        if (T1<T2):
            NA = NA - 1
            NB = NB - 1 
            NAB = NAB + 1
            t = t+ T1 
        else: 
            NA = NA + 1 
            NB = NB + 1 
            NAB = NAB - 1 
            t = t + T2
        times.append(t)
        A.append(NA)
        B.append(NB)
        AB.append(NAB)
    return times, A, B, AB

k1 = 2 
k2 = 1 
Ai = 1
Bi = 1 
ABi = 0 
tmax = 3

tgraph = np.linspace(0, tmax, 500)
res = solve_ivp(ode, [0, tmax], [Ai, Bi, ABi],t_eval = tgraph, args = (k1, k2))

plt.plot(res.t, res.y[2])
plt.xlabel("Time")
plt.ylabel("[AB]")
plt.show()

vols = [20000]

plt.plot(res.t, res.y[2], label = "Macroscopic ODE")

N = 5
for V in vols: 
    r1 = k1/V
    r2 = k2 
    
    NAi = int(Ai * V)
    NBi = int(Bi * V)
    NABi = int(ABi * V)
    
    everyAB = []
    
    for j in range(N):   
        times, A, B, AB = microscopic(r1, r2, NAi, NBi, NABi, tmax)
        print(j)
        print(times[:5])
        print(AB[:5])
        concAB = np.array(AB)
        
        collection_AB = []
        for time in tgraph: 
            place = np.searchsorted(times, time, side="right") - 1
            collection_AB.append(concAB[place])
        everyAB.append(collection_AB)
            
        plt.plot(times, concAB) # make these very thin
    
    everyAB = np.array(everyAB)
    mean_AB = np.mean(everyAB, axis = 0)
    #axis = 0 means it is going down the rows and is averaging each column 
    error = 2*np.std(everyAB, axis = 0)/np.sqrt(N)
    
    plt.plot(tgraph, mean_AB) # make this a thick line
    plt.fill_between(tgraph, mean_AB - error, mean_AB + error, alpha = 0.4)

    
    print("V =", V)
    print("AB:", AB[:15])
    #plt.step(times, concAB, label = f"Microscopic ODE Simulationfor Volume = {V}")
plt.xlabel("Time")
plt.ylabel("[AB]")
plt.legend()
plt.show()

# times 1, 2, 3, 4...1,000

actin = [1, 2, 5, 10]
ratios = np.linspace(0, 2, 50)

#need to pick an actin concentration, pick a ratio, determine the amount of profilin, run the simulaion until you reach equilibirum, and then calculate the fraction of actin that is still free, document it and go to the next part of the loop 
for conc in actin: 
    finalconcs = []
    for ratio in ratios: 
        profilin = ratio * conc 
        comp = 0
        res = solve_ivp(ode, [0, 3], [conc, profilin, comp], args = (2,1))
        finact = res.y[0, -1]
        fracfreeact = finact/conc
        finalconcs.append(fracfreeact)
    plt.plot(ratios, finalconcs, label = f"{conc}μM actin")
x=np.linspace(0, 1, 100)
y = 1 - x
plt.plot(x, y, ls='dotted', label=f"y=1-x")
plt.xlabel("Profilin:actin ratio")
plt.ylabel("Fraction free actin")
plt.legend()
plt.show()

#remember how it is important that the binding equilibirum depends on concentration (not only on the ratio)

        
        

            