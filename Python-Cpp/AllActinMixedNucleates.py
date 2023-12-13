import numpy as np
import sys
from ActinMixedNucleates import ActinMixedNucleates

# Parameters
a = 4e-3;
kbT = 4.1e-3;
mu = 0.01;
LBox = 3; # in um
Conc = 2; # in uM
ConcFormin = 0; # in uM
ConcArp23 = 0.5;
 # in uM

# Parameters from Kovar & Pollard paper for actin alone
# RETURN TO THE ACTUAL PARAMS LATER!
kplusDimer = 3.5e-3; # uM^(-1)*s^(-1) 
kminusDimer = 0.041; #s^(-1)
kplusTrimer = 13e-1; # uM^(-1)*s^(-1) 
kminusTrimer = 22; #s^(-1)
kplusBarbed = 1.6#11.6; # uM^(-1)*s^(-1) 
kminusBarbed = 1.4; #s^(-1)
kplusPointed = 1.3; #uM^(-1)*s^(-1)
kminusPointed = 0; #s^(-1)

# Formin rates
kForNuc = 2e-3; # uM^(-2)*s^(-1)
kplusFor = 29.1; # uM^(-1)*s^(-1)
kminusFor = 8.1e-2; # s^(-1)
ForminEnhance = 2;

# Arp 2/3 rates
kplusARF = 5.2e-2;
kMinusARF = 0.1#*3.4e-3;

# Convert to microscopic assuming well-mixed system
Volume = LBox**3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; # everything will be in s^(-1)
RxnRates=[kplusDimer*ConversionFactor, kminusDimer, kplusTrimer*ConversionFactor, kminusTrimer, \
    kplusBarbed*ConversionFactor, kminusBarbed, kplusPointed*ConversionFactor, kminusPointed];
RxnRatesFormin = [kForNuc*ConversionFactor**2, kplusFor*ConversionFactor, kminusFor, ForminEnhance];
RxnRatesArp23 = [kplusARF*ConversionFactor**2, kMinusARF];

Nmon = int(Conc*Volume/uMInvToMicron3);
print('Number of monomers %d' %Nmon)
NFormin = int(ConcFormin*Volume/uMInvToMicron3);
print('Number of formins %d' %NFormin)
NArp23 = int(ConcArp23*Volume/uMInvToMicron3);
print('Number of Arp 2/3 %d' %NArp23)

Lens=np.array([LBox,LBox,LBox]);
seed = int(sys.argv[1]);
nInFibers = [];#np.concatenate((3*np.ones(1,dtype=np.int64),\
    #4*np.ones(5,dtype=np.int64),5*np.ones(39,dtype=np.int64))); # initialize nucleates

nThr=1;
AllActin = ActinMixedNucleates(Nmon,nInFibers,Lens,a,kbT,mu, RxnRates, seed,nThr);
if (ConcFormin > 0):
    AllActin.InitializeFormins(NFormin,RxnRatesFormin);
if (ConcArp23 > 0):
    AllActin.InitializeArp(NArp23,RxnRatesArp23);

Tf = 40;
dt = 0.5;
nSteps = int(Tf/dt+1e-6);

NumFibers = np.zeros(nSteps,dtype=np.int64);
NumberPerFiber = np.array([],dtype=np.int64);
BranchedOrLinear = np.array([],dtype=bool);
BoundFormins = np.array([],dtype=bool);


for i in range(nSteps):
    AllActin.React(dt);
    NumOnEach = AllActin.NumMonOnEachFiber();
    NumberPerFiber = np.append(NumberPerFiber,NumOnEach)
    NumFibers[i]=len(NumOnEach)-3;
    BranchedOrLinear = np.append(BranchedOrLinear,AllActin.BranchedOrLinear(False))
    BoundFormins = np.append(BoundFormins,AllActin.BoundFormins())
    if (i==0):
        AllX = AllActin.getX();
    else:
        AllX = np.append(AllX,AllActin.getX(),axis=0);
    print('Time %f, Percent free %f' %((i+1)*dt, NumOnEach[0]/Nmon))
    
np.savetxt('NumFibs'+str(seed)+'.txt',NumFibers);
np.savetxt('StructInfo'+str(seed)+'.txt',NumberPerFiber);
np.savetxt('BoundFormins'+str(seed)+'.txt',BoundFormins);
np.savetxt('BranchedOrLinear'+str(seed)+'.txt',BranchedOrLinear);
#np.savetxt('AllX.txt',AllX);
#np.savetxt('NumFibs'+str(Conc)+'uM_Formin'+str((ConcFormin*1000))+'nM_Alpha'+str(ForminEnhance)+'_'+str(seed)+'.txt',NumFibers);
#np.savetxt('StructInfo'+str(Conc)+'uM_Formin'+str((ConcFormin*1000))+'nM_Alpha'+str(ForminEnhance)+'_'+str(seed)+'.txt',StructInfo);
#np.savetxt('BoundFormin'+str(Conc)+'uM_Formin'+str((ConcFormin*1000))+'nM_Alpha'+str(ForminEnhance)+'_'+str(seed)+'.txt',BoundFormins);

