import numpy as np
import sys
from AllTheActin import AllTheActin
from scipy.spatial import cKDTree
from math import pi

# Parameters
a = 4e-3;
kbT = 4.1e-3;
mu = 0.01;
RxnRadius = 2*a;

# Parameters from Kovar & Pollard paper
kp1 = 3.5e-6; # uM^(-1)*s^(-1) 
TwoMonOffRate = 0.041; #s^(-1)
kpB = 11.6; # uM^(-1)*s^(-1) 
BarbedUnbindingRate = 1.4; #s^(-1)
kpP = 1.3; #uM^(-1)*s^(-1)
PointedUnbindRate = 0.8; #s^(-1)

# Convert to microscopic assuming well-mixed system
uMInvToMicron3 = 1.0e15/(6.022e17);
TwoMonRate = kp1*uMInvToMicron3/(2*pi/3*RxnRadius**3);
print('Two mon rate %f' %TwoMonRate);
BarbedBindRate = kpB*uMInvToMicron3/(4*pi/3*RxnRadius**3);
print('Barbed bind rate %f' %BarbedBindRate);
PointedBindRate = kpP*uMInvToMicron3/(4*pi/3*RxnRadius**3);
print('Pointed bind rate %f' %PointedBindRate);
print('Diffusion coefficient %f' %(kbT/(6*pi*mu*a)))


Nmon = 1000;
LBox = 1;
Lens=np.array([LBox,LBox,LBox]);
seed = int(sys.argv[1]);
nInFibers = [];#np.concatenate((3*np.ones(1,dtype=np.int64),\
    #4*np.ones(5,dtype=np.int64),5*np.ones(39,dtype=np.int64))); # initialize nucleates
print(nInFibers)


# Reaction rates
#[_BindingRadius  _TwoMonRate _TwoMonOffRate _BarbedBindingRate  _BarbedUnbindingRate   _PointedBindingRate   _PointedUnbindingRate]
RxnRates = [RxnRadius, TwoMonRate, TwoMonOffRate, BarbedBindRate, BarbedUnbindingRate,PointedBindRate,PointedUnbindRate];

nThr=8;
AllActin = AllTheActin(Nmon,nInFibers,Lens,a,kbT,mu, RxnRates, seed,nThr);

Tf = 100;
logdt = 4;
dt = 10**(-logdt);
nSteps = int(Tf/dt+1e-6);
nSaves = 1000;
AllX = np.zeros(((nSaves+1)*Nmon,3));
AllIDs = np.zeros((nSaves+1)*Nmon,dtype=np.int64);
nFibs = np.zeros((nSteps+1),dtype=np.int64);

saveEvery = nSteps//nSaves;

for i in range(nSteps):
    if (i % saveEvery==0):
        print(str(i/nSteps*100)+'% done')
        saveIndex = i//saveEvery;
        AllX[saveIndex*Nmon:(saveIndex+1)*Nmon,:] = AllActin.getX();
        AllIDs[saveIndex*Nmon:(saveIndex+1)*Nmon] = AllActin.getStructureIDs();
        print('Number fibers %d' %nFibs[i])
    nFibs[i+1]=AllActin.React(dt);
    AllActin.Diffuse(dt);
    
AllX[nSaves*Nmon:(nSaves+1)*Nmon,:] = AllActin.getX();  
AllIDs[nSaves*Nmon:(nSaves+1)*Nmon] = AllActin.getStructureIDs(); 
np.savetxt('Dt'+str(logdt)+'X_'+str(seed)+'.txt',AllX);
np.savetxt('Dt'+str(logdt)+'IDs_'+str(seed)+'.txt',AllIDs);
np.savetxt('Dt'+str(logdt)+'Fibers_'+str(seed)+'.txt',nFibs);
