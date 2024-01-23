import numpy as np
import sys
from ActinMixedNucleates import ActinMixedNucleates

# Parameters
a = 4e-3;
spacing = 0.5; #units of a
kbT = 4.1e-3;
mu = 0.01;
LBox = 3; # in um
Conc = 2; # in uM
ConcFormin = 0.1; # in uM
ConcArp23 = 0;
 # in uM

# Parameters from Kovar & Pollard paper for actin alone
# RETURN TO THE ACTUAL PARAMS LATER!
kplusDimer = 3.5e-3; # uM^(-1)*s^(-1) 
kminusDimer = 0.041; #s^(-1)
kplusTrimer = 13e-2; # uM^(-1)*s^(-1) 
kminusTrimer = 22; #s^(-1)
kplusBarbed = 1.6#11.6; # uM^(-1)*s^(-1) 
kminusBarbed = 1.4; #s^(-1)
kplusPointed = 1.3; #uM^(-1)*s^(-1)
kminusPointed = 0.8; #s^(-1)

# Formin rates
kForNuc = 2e-3; # uM^(-2)*s^(-1)
kplusFor = 0*29.1; # uM^(-1)*s^(-1)
kminusFor = 0*8.1e-2; # s^(-1)
ForminEnhance = 1.0;

# Arp 2/3 rates
kplusARF = 5.2e-2;
kMinusARF = 0.5;

# Convert to microscopic assuming well-mixed system
Volume = LBox**3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; # everything will be in s^(-1)
SpontaneousRxnRates=[kplusDimer*ConversionFactor, kminusDimer, kplusTrimer*ConversionFactor, kminusTrimer, \
    kplusBarbed*ConversionFactor, kminusBarbed, kplusPointed*ConversionFactor, kminusPointed];

RxnRatesArp23 = [kplusARF*ConversionFactor**2, kMinusARF];

Nmon = int(Conc*Volume/uMInvToMicron3);
print('Number of monomers %d' %Nmon)
NFormin = int(ConcFormin*Volume/uMInvToMicron3);
print('Number of formins %d' %NFormin)
NArp23 = int(ConcArp23*Volume/uMInvToMicron3);
print('Number of Arp 2/3 %d' %NArp23)

Lens=np.array([LBox,LBox,LBox]);
seed = int(sys.argv[1]);


nThr=1;
AllActin = ActinMixedNucleates(Nmon,Lens,SpontaneousRxnRates,a,spacing,kbT,mu, seed,nThr);
if (ConcFormin > 0):
    NBarbed = [NFormin];
    BarbedOnOff = [kplusFor*ConversionFactor, kminusFor];
    AlphaDimersMinus = [1,0];
    AlphaDimersPlus = [1,kForNuc*ConversionFactor/kplusDimer];
    AlphaTrimersMinus = [1,0];
    AlphaTrimersPlus = [1, 10]#(kplusBarbed+kplusPointed)/kplusTrimer];
    print(AlphaTrimersPlus)
    AlphaBarbed = [1, ForminEnhance];
    AllActin.InitializeBarbedBinders(NBarbed,BarbedOnOff,AlphaDimersMinus,AlphaTrimersMinus);
    AllActin.InitializeRateMatrices(AlphaDimersPlus,AlphaTrimersPlus,AlphaBarbed);
if (ConcArp23 > 0):
    AllActin.InitializeBranchers(NArp23,RxnRatesArp23);

Tf = 40;
dt = 0.5;
nSteps = int(Tf/dt+1e-6);

NumFibers = np.zeros(nSteps,dtype=np.int64);
FreeMonomers = np.zeros(nSteps,dtype=np.int64);
NumberPerFiber = np.array([],dtype=np.int64);
BranchedOrLinear = np.array([],dtype=bool);
BoundProteins = np.array([],dtype=bool);

for i in range(nSteps):
    AllActin.React(dt);
    NumOnEach = AllActin.NumMonOnEachFiber();
    NumberPerFiber = np.append(NumberPerFiber,NumOnEach)
    NumFibers[i] = len(NumOnEach);
    FreeMonomers[i] = AllActin.nFreeMonomers();
    BranchedOrLinear = np.append(BranchedOrLinear,AllActin.BranchedOrLinear(False))
    BoundProteins = np.append(BoundProteins,AllActin.BoundBarbedStates())
    print('Time %f, Percent free %f' %((i+1)*dt, FreeMonomers[i]/Nmon))
    
#np.savetxt('ForminArpNoUnbind_NumFibs'+str(seed)+'.txt',NumFibers);
#np.savetxt('ForminArpNoUnbind_StructInfo'+str(seed)+'.txt',NumberPerFiber);
#np.savetxt('ForminArpNoUnbind_BoundFormins'+str(seed)+'.txt',BoundFormins);
#np.savetxt('ForminArpNoUnbind_BranchedOrLinear'+str(seed)+'.txt',BranchedOrLinear);
np.savetxt('ForminOnly_FreeMons'+str(seed)+'.txt',FreeMonomers);
np.savetxt('ForminOnly_NumFibs'+str(seed)+'.txt',NumFibers);
np.savetxt('ForminOnly_StructInfo'+str(seed)+'.txt',NumberPerFiber);
np.savetxt('ForminOnly_BoundFormins'+str(seed)+'.txt',BoundProteins);
np.savetxt('ForminOnly_BranchedOrLinear'+str(seed)+'.txt',BranchedOrLinear);

