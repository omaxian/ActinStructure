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

# NEW FILE for stuff below here
Conc = 5; # in uM
ConcArp = 0;
Tf = 28800
LBox = 5
SeedConc = 0
#for seed in range(nTrials):
    #FileName = 'Tf'+str(Tf)+'_Box'+str(LBox)+'_Actin'+str(Conc)+'uM_Seed'+str(SeedConc)+'uM_Arp'+str(int(ConcArp*1000))+'nM_'+str(seed)+'.txt';


run_1 = np.loadtxt("/Users/sophiachainani/Documents/ActinStructure/Python-Cpp/Actin5run_1/NumFibsTf28800_Box5_Actin5uM_Seed_0_KProf1_Prof0uM_Arp0nM_Formin0em4uM_1.txt")
run_2 = np.loadtxt("/Users/sophiachainani/Documents/ActinStructure/Python-Cpp/Actin5run_2/NumFibsTf28800_Box5_Actin5uM_Seed_0_KProf1_Prof0uM_Arp0nM_Formin0em4uM_1.txt")
run_3 = np.loadtxt("/Users/sophiachainani/Documents/ActinStructure/Python-Cpp/Actin5run_3/NumFibsTf28800_Box5_Actin5uM_Seed_0_KProf1_Prof0uM_Arp0nM_Formin0em4uM_1.txt")

total_runs = np.array([run_1, run_2, run_3])
print(total_runs.shape)

mean_numfibs = np.mean(total_runs, axis = 0)
N = 3
error = 2 * np.std(total_runs, axis = 0)/np.sqrt(N)
time = np.arange(2880)*10

plt.plot(time, mean_numfibs, label = "Mean number of fibers over time")
plt.fill_between(time, mean_numfibs - error, mean_numfibs + error, alpha = 0.2)
plt.xlim(0,1000)

plt.xlabel("Time")
plt.ylabel("Number of fibers")
plt.legend()

print(mean_numfibs[:10])
print(error[:10])

#going to do the stats part here but just doing random branches and lengths rn will import code later 
branches = [0, 2, 1, 1, 0, 2, 1, 1, 1, 0]
lengths = []
count = 0
branch_counts = []
current_lengths = []
in_struc = False 
for val, length in zip(branches, lengths): 
    if val == 2:
        in_struc = True 
        #mother/start of structure
        count = 0
    elif val == 1: 
       #daughter branch
        count +=1
        current_lengths.append(length)
    elif val == 0 and in_struc == True:
        print(f"structure had {count} branches")
        branch_counts.append(count)
        count = 0
    else: 
        continue 
print(branch_counts)

        
        


