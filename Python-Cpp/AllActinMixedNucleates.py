import numpy as np
import sys
from ActinMixedNucleates import ActinMixedNucleates

# Parameters
a = 4e-3;
kbT = 4.1e-3;
mu = 0.01;
LBox = 5; # in um
Conc = 5; # in uM
ConcFormin = 0.00012; # in uM

# Parameters from Kovar & Pollard paper for actin alone
# RETURN TO THE ACTUAL PARAMS LATER!
kplusDimer = 3.5e-6; # uM^(-1)*s^(-1) 
kminusDimer = 0.041; #s^(-1)
kplusTrimer = 13e-5; # uM^(-1)*s^(-1) 
kminusTrimer = 22; #s^(-1)
kplusBarbed = 11.6; # uM^(-1)*s^(-1) 
kminusBarbed = 1.4; #s^(-1)
kplusPointed = 1.3; #uM^(-1)*s^(-1)
kminusPointed = 0.8; #s^(-1)

# Formin rates
kForNuc = 2e-4; # uM^(-2)*s^(-1)
kplusFor = 29.1; # uM^(-1)*s^(-1)
kminusFor = 8.1e-5; # s^(-1)
ForminEnhance = 1;

# Convert to microscopic assuming well-mixed system
Volume = LBox**3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; # everything will be in s^(-1)
RxnRates=[kplusDimer*ConversionFactor, kminusDimer, kplusTrimer*ConversionFactor, kminusTrimer, \
    kplusBarbed*ConversionFactor, kminusBarbed, kplusPointed*ConversionFactor, kminusPointed];
RxnRatesFormin = [kForNuc*ConversionFactor**2, kplusFor*ConversionFactor, kminusFor, ForminEnhance];

Nmon = int(Conc*Volume/uMInvToMicron3);
print('Number of monomers %d' %Nmon)
NFormin = int(ConcFormin*Volume/uMInvToMicron3);
print('Number of formins %d' %NFormin)

Lens=np.array([LBox,LBox,LBox]);
seed = int(sys.argv[1]);
nInFibers = [];#np.concatenate((3*np.ones(1,dtype=np.int64),\
    #4*np.ones(5,dtype=np.int64),5*np.ones(39,dtype=np.int64))); # initialize nucleates

nThr=1;
AllActin = ActinMixedNucleates(Nmon,nInFibers,Lens,a,kbT,mu, RxnRates, seed,nThr);
if (ConcFormin > 0):
    AllActin.InitializeFormins(NFormin,RxnRatesFormin);

Tf = 2000;
dt = 1;
nSteps = int(Tf/dt+1e-6);
StructInfo = np.array([],dtype=np.int64);
BoundFormins = np.array([],dtype=np.int64);
NumFibers = np.zeros(nSteps,dtype=np.int64);

for i in range(nSteps):
    NumFibers[i]=AllActin.React(dt);
    ThisStruct = AllActin.getStructureInfo()
    StructInfo = np.append(StructInfo,ThisStruct);
    BoundFormins = np.append(BoundFormins,AllActin.getBoundFormins());
    print('Time %f, Percent free %f' %((i+1)*dt, ThisStruct[0]/Nmon))
    
#np.savetxt('NumFibs.txt',NumFibers);
#np.savetxt('StructInfo.txt',StructInfo);
#np.savetxt('BoundFormins.txt',BoundFormins);
np.savetxt('NumFibs'+str(Conc)+'uM_Formin'+str((ConcFormin*1000))+'nM_Alpha'+str(ForminEnhance)+'_'+str(seed)+'.txt',NumFibers);
np.savetxt('StructInfo'+str(Conc)+'uM_Formin'+str((ConcFormin*1000))+'nM_Alpha'+str(ForminEnhance)+'_'+str(seed)+'.txt',StructInfo);
np.savetxt('BoundFormin'+str(Conc)+'uM_Formin'+str((ConcFormin*1000))+'nM_Alpha'+str(ForminEnhance)+'_'+str(seed)+'.txt',BoundFormins);

