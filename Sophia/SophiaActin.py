#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 08:06:58 2026

@author: sophiachainani
"""
Conc = 5; # in uM
ConcArp = 80 * 10**(-3);
Tf = 2880
LBox = 5
SeedConc = 0
nTrials = 5

import sys

# Option A: Append to the end of the search path
sys.path.append("/Users/sophiachainani/Documents/ActinStructure/Python-Cpp")
sys.path.append("/Users/ondrejmaxian/Documents/ACTIN/ActinStructure/Python-Cpp")

import numpy as np 
import matplotlib.pyplot as plt
import InitialActinSim as sim


#for seed in range(nTrials):
   #FileName = 'Tf'+str(Tf)+'_Box'+str(LBox)+'_Actin'+str(Conc)+'uM_Seed'+str(SeedConc)+'uM_Arp'+str(int(ConcArp*1000))+'nM_'+str(seed)+'.txt';
#will incorporate this later on - having trouble because couldn't change the Arp concentration successfully in the old Actin Sim file so had to make a copy, but need help changing the labeling  



num_fibers = sim.NumFibers
number_per_fiber = sim.NumberPerFiber
branched_or_linear = sim.BranchedOrLinear
free_monomers = sim.FreeMonomers
NArp23 = sim.NArp23
total_actin = sim.Nmon

print(num_fibers.shape)
print(len(number_per_fiber))
print(len(branched_or_linear))

avg_length_over_time = []
branch_count_over_time = []
num_struc_over_time = []
avg_branches_per_struc_over_time = []
branch_density_over_time = []
bound_arp_over_time = []
count = 0
num_monomers_per_len = 2000
branched_monomers_over_time = []
#normally, branches_here = branched_or_linear[count:count+n] and lengths_here = number_per_fiber[count:count+n]
# Example 1 had branches_here = np.array([2, 1, 1, 0, 0, 2, 1, 0]) and lengths_here = np.array([1000, 4000, 800, 700, 500, 3000, 900, 600])
# Example 2 had branches_here = np.array([0, 2, 1, 0, 2, 1, 1, 1]) and lengths_here = np.array([500, 200, 1000, 800, 3000, 400, 600, 100])
for i in range(len(num_fibers)):
    n = num_fibers[i]
    branches_here = branched_or_linear[count:count+n]
    lengths_here = number_per_fiber[count:count+n]
    count += n
    total = 0
    branch_totals = []
    current_lengths = []
    avg_branch_lengths = []
    curr_struc_num_monomers = []
    branch_densities = []
    branched_monomers_this_time = 0
    
    #total is our number of branches (daughters) in our current structure 
    for j in range(len(lengths_here)):
        val = branches_here[j]
        length = lengths_here[j]
        if val == 2 or val==0:
            total = 0
            current_lengths = [] 
            curr_struc_num_monomers = [length]
        elif val == 1: 
            total += 1 
            current_lengths.append(length)
            curr_struc_num_monomers.append(length)
        if (j == len(branches_here) - 1) or (branches_here[j+1] !=1):
            branch_totals.append(total)
            if len(current_lengths) > 0:
                avg_branch_lengths.append(np.mean(current_lengths))
            struc_length = np.sum(curr_struc_num_monomers)/num_monomers_per_len
            if struc_length > 0:
                density = total/struc_length
                branch_densities.append(density)
            if total>0:
                branched_monomers_this_time += np.sum(curr_struc_num_monomers)
            total = 0 
            current_lengths = []
            curr_struc_num_monomers = []
    if len(avg_branch_lengths) > 0:
        avg_length_over_time.append(np.mean(avg_branch_lengths))
    else:
        avg_length_over_time.append(0)
    if len(branch_totals)>0:
        num_struc_over_time.append(len(branch_totals))
        branch_count_over_time.append(np.sum(branch_totals))
        avg_branches_per_struc_over_time.append(np.mean(branch_totals))
        bound_arp_over_time.append(np.sum(branch_totals))
    else:
        num_struc_over_time.append(0)
        branch_count_over_time.append(0)
        avg_branches_per_struc_over_time.append(0)
        bound_arp_over_time.append(0)
    if len(branch_densities)>0:
       branch_density_over_time.append(np.mean(branch_densities)) 
    else:
        branch_density_over_time.append(0)
    branched_monomers_over_time.append(branched_monomers_this_time)
    
free_arp_over_time = NArp23 - np.array(bound_arp_over_time)
# num of branched structures = len(branch_totals), avg per struc = np.mean(branch_totals), total num of branches = np.sum(branch_totals)

        
time = np.arange(len(num_fibers)) * 10

plt.figure()
plt.plot(time, num_struc_over_time)
plt.xlabel("Time")
plt.ylabel("Number of branched structures")
plt.show()

plt.figure()
plt.plot(time, avg_branches_per_struc_over_time)
plt.xlabel("Time")
plt.ylabel("Average number of branches per structure")
plt.show()

plt.figure()
plt.plot(time, avg_length_over_time)
plt.xlabel("Time")
plt.ylabel("Average branch length")
plt.show()
    
plt.figure()
plt.plot(time, branch_density_over_time)
plt.xlabel("Time")
plt.ylabel("Branch density")   
plt.show()

plt.figure()
plt.plot(time, free_monomers)
plt.xlabel("Time")
plt.ylabel("Free monomers")
plt.show()

plt.figure()
plt.plot(time, free_arp_over_time)
plt.xlabel("Time")
plt.ylabel("Free Arp 2/3")
plt.show()


print(free_monomers)
print(num_fibers)
print(time)


fraction_bound_arp = np.array(bound_arp_over_time)/NArp23

polymerized_actin_over_time = total_actin - free_monomers

fraction_branched_actin = np.divide(np.array(branched_monomers_over_time), polymerized_actin_over_time, out = np.zeros(len(polymerized_actin_over_time)), where=polymerized_actin_over_time !=0)

plt.figure()
plt.plot(time, fraction_bound_arp)
plt.xlabel("Time")
plt.ylabel("Overall branch density")
plt.show()

plt.figure()
plt.plot(time, fraction_branched_actin)
plt.xlabel("Time")
plt.ylabel("Fraction of polymerized actin in branched structures")
plt.show()

#checking to be sure that the length of all arrays is the same 
print("Test to check length of each array is the same")
print(len(time))
print(len(num_struc_over_time))
print(len(avg_branches_per_struc_over_time))
print(len(avg_length_over_time))
print(len(branch_density_over_time))
print(len(free_monomers))
print(len(free_arp_over_time))
print(len(fraction_bound_arp))
print(len(fraction_branched_actin))


#trying to get all of the plots in one panel
fig, axes = plt.subplots(2, 4, figsize = (16, 8))

axes[0,0].plot(time, num_struc_over_time)
axes[0,0].set_xlabel("Time")
axes[0,0].set_ylabel("Number of structures")

axes[0,1].plot(time, avg_branches_per_struc_over_time)
axes[0,1].set_title("Average branches per structure")
axes[0,1].set_xlabel("Time")

axes[0,2].plot(time, avg_length_over_time)
axes[0,2].set_title("Average branch length")
axes[0,2].set_xlabel("Time")

axes[0,3].plot(time, branch_density_over_time)
axes[0,3].set_title("Branch density")
axes[0,3].set_xlabel("Time")

axes[1,0].plot(time, free_monomers)
axes[1,0].set_title("Free monomers")
axes[1,0].set_xlabel("Time")

axes[1,1].plot(time, free_arp_over_time)
axes[1,1].set_title("Free Arp2/3")
axes[1,1].set_xlabel("Time")

axes[1,2].plot(time, fraction_bound_arp)
axes[1,2].set_title("Fraction of Arp bound")
axes[1,2].set_xlabel("Time")

axes[1,3].plot(time, fraction_branched_actin)
axes[1,3].set_title("Fraction of polymerized actin in branched structures")
axes[1,3].set_xlabel("Time")

plt.tight_layout()
plt.show()



















# used to check my handwritten examples
#print("Testing")
#print(branch_totals)
#print(len(branch_totals))
#print(np.mean(branch_totals))
#print(branch_densities)







#OLD STUFF
# mean_numfibs = np.mean(total_runs, axis = 0)
   # N = 3
    #error = 2 * np.std(total_runs, axis = 0)/np.sqrt(N)
    #time = np.arange(2880)*10
    
    #plt.plot(time, mean_numfibs, label = "Mean number of fibers over time")
    #plt.fill_between(time, mean_numfibs - error, mean_numfibs + error, alpha = 0.2)
    #plt.xlim(0,1000)
    
    #plt.xlabel("Time")
    #plt.legend()
    
    #print(mean_numfibs[:10])
    #print(error[:10])


#Old filtering code 
#  for val, length in zip(branches_here, lengths_here): 
    #if val == 2:
       # curr_struc_num_monomers = [length]
        #need to finish it 
       # if in_struc == True: 
           # branch_totals.append(total)
            #if total>0:
             #   avg_branch_lengths.append(np.mean(current_lengths))
        #mother/start of structure
       # in_struc = True
        #total = 0
        #current_lengths = [] # FIx: mother should be part of current lengths
    #elif val == 1 and in_struc == True: 
       #daughter branch
       #total +=1
        #curr_struc_num_monomers.append(length)
        #current_lengths.append(length)
    #elif val == 0 and in_struc == True:
        #end of current branch
        #density = np.sum(curr_struc_num_monomers) / 2000
        #branch_densities.append(density)
        #branch_totals.append(total)
        #if total>0:
        #    avg_branch_lengths.append(np.mean(current_lengths))
        #total = 0
        #current_lengths = []
        #in_struc = False 
        
        


