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
ConcArp23 = 0.1;
ConcProf = 0.05;
 # in uM

# Parameters from Kovar & Pollard paper for actin alone
# RETURN TO THE ACTUAL PARAMS LATER!
kplusDimer = 3.5e-3; # uM^(-1)*s^(-1) 
kminusDimer = 0.041; #s^(-1)
kplusTrimer = 13e-2; # uM^(-1)*s^(-1) 
kminusTrimer = 22; #s^(-1)
kplusBarbed = 1.6#11.6; # uM^(-1)*s^(-1) 
kminusBarbed = 0*1.4; #s^(-1)
kplusPointed = 1.3; #uM^(-1)*s^(-1)
kminusPointed = 0*0.8; #s^(-1)

# Formin rates
kForNuc = 2e-3; # uM^(-2)*s^(-1)
kplusFor = 5; # uM^(-1)*s^(-1)
kminusFor = 8.1e-2; # s^(-1)

# Arp 2/3 rates
kplusARF = 10e-2;
kMinusARF = 0.5;

# Profilin rate
ProfEq = 5;

# Convert to microscopic assuming well-mixed system
Volume = LBox**3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; # everything will be in s^(-1)
SpontaneousRxnRates=[kplusDimer*ConversionFactor, kminusDimer, kplusTrimer*ConversionFactor, kminusTrimer, \
    kplusBarbed*ConversionFactor, kminusBarbed, kplusPointed*ConversionFactor, kminusPointed];
ForFac = kForNuc*ConversionFactor/kplusDimer;

Nmon = int(Conc*Volume/uMInvToMicron3);
print('Number of monomers %d' %Nmon)
NFormin = int(ConcFormin*Volume/uMInvToMicron3);
print('Number of formins %d' %NFormin)
NArp23 = int(ConcArp23*Volume/uMInvToMicron3);
print('Number of Arp 2/3 %d' %NArp23)
NProf = int(ConcProf*Volume/uMInvToMicron3);
print('Number of profilin %d' %NProf)

Lens=np.array([LBox,LBox,LBox]);
seed = int(sys.argv[1]);


nThr=1;
AllActin = ActinMixedNucleates(Nmon,Lens,SpontaneousRxnRates,a,spacing,kbT,mu, seed,nThr);
if (ConcFormin > 0):
    NBarbed = [NFormin, NFormin-300];
    BarbedOnOff = [kplusFor*ConversionFactor, kminusFor, 0.6*kplusFor*ConversionFactor, 0.2*kminusFor];
    AlphaDimersMinus = [1,0.2,0.7];
    AlphaTrimersMinus = [1,0.4,1.3];
    AlphaBarbedMinus = [1, 0.1,1.2];
    AllActin.InitializeBarbedBinders(NBarbed,BarbedOnOff,AlphaDimersMinus,AlphaTrimersMinus,AlphaBarbedMinus);
    #AlphaDimersPlus = [1,];
    #AlphaTrimersPlus = [1, 10]
    #AlphaBarbedPlus = [1, ForminEnhance];
    #AllActin.InitializeRateMatrices(AlphaDimersPlus,AlphaTrimersPlus,AlphaBarbedPlus);
if (ConcProf > 0):
    NMonProts = [NProf, NProf-200]
    KsMon = [ProfEq*ConversionFactor, 0.8*ConversionFactor];
    AlphasPointed = [1,0.5, 0.1];
    AllActin.InitializeMonomerBinders(NMonProts,KsMon,AlphasPointed);
    # Re-initialize the rate matrices with all rates
    #if (ConcFormin == 0):
        #AlphaDimersPlus = [1,1.7];
        #AlphaTrimersPlus = [1, 0.4]#(kplusBarbed+kplusPointed)/kplusTrimer];
        #AlphaBarbedPlus = [1, 2];
    #else:
    AlphaDimersPlus = [1,ForFac,1.8*ForFac, 1.7, 0.5*ForFac, 0.9*ForFac,  0.6, 3*ForFac, 1.2*ForFac];
    AlphaTrimersPlus = [1, 10,3,0.4, 2, 1.9, 1.5, 2, 0.5]
    AlphaBarbedPlus = [1, 2, 1.5, 2, 2.5, 1.25, 0.1, 0.7, 1.1];
    AllActin.InitializeRateMatrices(AlphaDimersPlus,AlphaTrimersPlus,AlphaBarbedPlus);
if (ConcArp23 > 0):
    RxnRatesArp23 = [kplusARF*ConversionFactor**2, kMinusARF];
    if (ConcProf > 0):
        AlphasArp23 = [1, 1.3, 0.6];
    else:
        AlphasArp23 = [1];
    AllActin.InitializeBranchers(NArp23,RxnRatesArp23,AlphasArp23);   
    

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

