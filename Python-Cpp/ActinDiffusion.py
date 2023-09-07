import numpy as np
import sys
from AllTheActin import AllTheActin
from scipy.spatial import cKDTree

a = 4e-3;
kbT = 4.1e-3;
mu = 0.068*100;
Nmon = 1000;
LBox = 1;
Lens=np.array([LBox,LBox,LBox]);
seed = int(sys.argv[1]);
np.random.seed(seed);
#x = 0.5+np.array([np.arange(Nmon)*2.5*a]);
#y = 0.5+np.zeros((Nmon,1));
#z = y;
#X = np.concatenate((x.T,y,z),axis=1)
X = np.random.rand(Nmon,3)*Lens;

# Reaction rates
#[_BindingRadius  _TwoMonRate  _BarbedBindingRate  _BarbedUnbindingRate   _PointedBindingRate   _PointedUnbindingRate]
RxnRates = [10*a, 10, 0, 5, 0, 7];

AllActin = AllTheActin(X,Lens,a,kbT,mu, RxnRates, seed);

Tf = 1;
dt = 0.01;
nSteps = int(Tf/dt+1e-6);
#AllX = np.zeros(((nSteps+1)*Nmon,3));
#AllIDs = np.zeros((nSteps+1)*Nmon,dtype=np.int64);
nMonomers = np.zeros((nSteps+1),dtype=np.int64);
nMonomers[0]=Nmon;

#rcut = 4*a;
#KDTree = cKDTree(X,boxsize=Lens);
#Neighbors = KDTree.query_pairs(rcut,output_type='ndarray');
#np.savetxt('NeighborsKD.txt',Neighbors)

for i in range(nSteps):
    #AllX[i*Nmon:(i+1)*Nmon,:] = AllActin.getX();
    #AllIDs[i*Nmon:(i+1)*Nmon] = AllActin.getStructureIDs();
    if (i % (nSteps/100) == 0):
        print(str(i/nSteps*100)+'% done')
    AllActin.Diffuse(dt);
    nMonomers[i+1]=AllActin.React(dt);
    
#AllX[nSteps*Nmon:(nSteps+1)*Nmon,:] = AllActin.getX();  
#AllIDs[nSteps*Nmon:(nSteps+1)*Nmon] = AllActin.getStructureIDs(); 
#np.savetxt('X_'+str(seed)+'.txt',AllX);
np.savetxt('R1Dt2nFreeMon_'+str(seed)+'.txt',nMonomers);
