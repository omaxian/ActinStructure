import numpy as np
import sys
from AllTheActin import AllTheActin

a = 4e-3;
kbT = 4.1e-3;
mu = 1;
Nmon = 1000;
LBox = 1;
seed = int(sys.argv[1]);
np.random.seed(seed);
X = np.random.rand(Nmon,3);
AllActin = AllTheActin(X,a,kbT,mu);

dt = 1e-4;
nSteps = 1000;
X = np.zeros(((nSteps+1)*Nmon,3));

for i in range(nSteps):
    X[i*Nmon:(i+1)*Nmon,:] = AllActin.getX();
    W = np.random.randn(Nmon*3);
    AllActin.Diffuse(dt,W)
    
X[nSteps*Nmon:(nSteps+1)*Nmon,:] = AllActin.getX();   
np.savetxt('X_'+str(seed)+'.txt',X);
