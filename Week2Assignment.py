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
tmax = 15

tgraph = np.linspace(0, tmax, 500)
res = solve_ivp(ode, [0, tmax], [Ai, Bi, ABi],t_eval = tgraph, args = (k1, k2))

plt.plot(res.t, res.y[2])
plt.xlabel("Time")
plt.ylabel("[AB]")
plt.show()

vols = [20, 200, 2000, 20000]

plt.plot(res.t, res.y[2], label = "Macroscopic ODE")

for V in vols: 
    r1 = k1/V
    r2 = k2 
    
    NAi = int(Ai * V)
    NBi = int(Bi * V)
    NABi = int(ABi * V)
    
    times, A, B, AB = microscopic(r1, r2, NAi, NBi, NABi, tmax)
    concAB = np.array(AB)/V
    print("V =", V)
    print("AB:", AB[:15])
    plt.step(times, concAB, label = f"Microscopic ODE Simulation for Volume = {V}")
plt.xlabel("Time")
plt.ylabel("[AB]")
plt.legend()
plt.show()

            