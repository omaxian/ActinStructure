import numpy as np
import sys
from AllTheActin import AllTheActin
from scipy.spatial import cKDTree

a = 4e-3;
kbT = 4.1e-3;
mu = 1;
Nmon = 1000;
LBox = 1;
Lens=np.array([LBox,LBox,LBox]);
seed = int(sys.argv[1]);
np.random.seed(seed);
#x = np.array([np.arange(Nmon)*2.5*a]);
#y = np.zeros((Nmon,1));
#z = y;
#X = np.concatenate((x.T,y,z),axis=1)
X = np.random.rand(Nmon,3)*Lens;

AllActin = AllTheActin(X,Lens,a,kbT,mu, seed);

dt = 1e-4;
nSteps = 1000;
AllX = np.zeros(((nSteps+1)*Nmon,3));
AllIDs = np.zeros((nSteps+1)*Nmon,dtype=np.int64);

#rcut = 4*a;
#KDTree = cKDTree(X,boxsize=Lens);
#Neighbors = KDTree.query_pairs(rcut,output_type='ndarray');
#np.savetxt('NeighborsKD.txt',Neighbors)

for i in range(nSteps):
    AllX[i*Nmon:(i+1)*Nmon,:] = AllActin.getX();
    AllIDs[i*Nmon:(i+1)*Nmon] = AllActin.getStructureIDs();
    AllActin.React(dt);
    AllActin.Diffuse(dt);
    
AllX[nSteps*Nmon:(nSteps+1)*Nmon,:] = AllActin.getX();  
AllIDs[nSteps*Nmon:(nSteps+1)*Nmon] = AllActin.getStructureIDs(); 
np.savetxt('X_'+str(seed)+'.txt',AllX);
np.savetxt('IDs_'+str(seed)+'.txt',AllIDs);
