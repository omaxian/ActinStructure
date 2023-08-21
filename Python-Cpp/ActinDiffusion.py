import numpy as np
import sys
from AllTheActin import AllTheActin

Nmon = 28;
AllActin = AllTheActin();


seed = int(sys.argv[1]);
np.random.seed(seed);
dt = 1e-4;
nSteps = 1;
X = np.zeros(((nSteps+1)*Nmon,3));

for i in range(nSteps):
    X[i*Nmon:(i+1)*Nmon,:] = AllActin.getX();
    W = np.random.randn(Nmon*3);
    np.savetxt('W1.txt',W)
    Wdr = np.random.randn(12);
    drift = False;
    np.savetxt('WDrift.txt',Wdr)
    AllActin.Diffuse(dt,W,Wdr,drift)
    
X[nSteps*Nmon:(nSteps+1)*Nmon,:] = AllActin.getX();   
np.savetxt('X_'+str(seed)+'.txt',X);
