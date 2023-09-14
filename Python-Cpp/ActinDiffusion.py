import numpy as np
import sys
from AllTheActin import AllTheActin
from scipy.spatial import cKDTree

a = 4e-3;
kbT = 4.1e-3;
mu = 0.068;
Nmon = 1000;
LBox = 1;
Lens=np.array([LBox,LBox,LBox]);
seed = int(sys.argv[1]);
np.random.seed(seed);
#x = 0.5+np.array([np.arange(Nmon)*2*a]);
#y = 0.5+np.zeros((Nmon,1));
#z = y;
#X = np.concatenate((x.T,y,z),axis=1)
X = np.random.rand(Nmon,3)*Lens;

# Reaction rates
#[_BindingRadius  _TwoMonRate _TwoMonOffRate _BarbedBindingRate  _BarbedUnbindingRate   _PointedBindingRate   _PointedUnbindingRate]
RxnRates = [2*a, 12, 5, 2, 1.5,4,2.5];

AllActin = AllTheActin(X,Lens,a,kbT,mu, RxnRates, seed);

Tf = 5;
dt = 0.01;
nSteps = int(Tf/dt+1e-6);
nSaves = 500;
AllX = np.zeros(((nSaves+1)*Nmon,3));
AllIDs = np.zeros((nSaves+1)*Nmon,dtype=np.int64);
nFibs = np.zeros((nSteps+1),dtype=np.int64);

#rcut = 4*a;
#KDTree = cKDTree(X,boxsize=Lens);
#Neighbors = KDTree.query_pairs(rcut,output_type='ndarray');
#np.savetxt('NeighborsKD.txt',Neighbors)
saveEvery = nSteps//nSaves;

for i in range(nSteps):
    if (i % saveEvery==0):
        print(str(i/nSteps*100)+'% done')
        saveIndex = i//saveEvery;
        AllX[saveIndex*Nmon:(saveIndex+1)*Nmon,:] = AllActin.getX();
        AllIDs[saveIndex*Nmon:(saveIndex+1)*Nmon] = AllActin.getStructureIDs();
    nFibs[i+1]=AllActin.React(dt);
    AllActin.Diffuse(dt);
    
AllX[nSaves*Nmon:(nSaves+1)*Nmon,:] = AllActin.getX();  
AllIDs[nSaves*Nmon:(nSaves+1)*Nmon] = AllActin.getStructureIDs(); 
np.savetxt('X_'+str(seed)+'.txt',AllX);
np.savetxt('IDs_'+str(seed)+'.txt',AllIDs);
np.savetxt('Fibers_'+str(seed)+'.txt',nFibs);
