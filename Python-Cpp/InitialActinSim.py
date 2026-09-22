import numpy as np
import sys
from ActinMixedNucleates import ActinMixedNucleates

# Parameters
Conc = 5; # in uM
ConcProf = 0;
ConcArp = 80e-3;
SeedConc = 0;
ConcFormin = 0;

a = 4e-3;
kbT = 4.1e-3;
spacing = 0.5; # units of a
mu = 0.01;
LBox = 10;
if (ConcProf<1.5):
    LBox = 5;

# Parameters from Kovar & Pollard paper for actin alone
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
kminusFor = 5e-4; # s^(-1)
ForminAlphaNoProf = 0.5;

# Arp 2/3 rates
kplusARF = 2e-4*spacing*a; # This is now a rate in uM^2 per monomer of mother (default 5e-3)
kMinusARF = 3.4e-3;

# Profilin equilibrium constant
ProfEq = 1; # uM^(-1) (DEFAULT IS 4.5 but Aidan's data is 1)
AlphaWithProf = 1;
ForminAlphaWithProf = 3;
AlphaPtdProf = 0; # Profilin blocks interaction at ptd end.

# Convert to microscopic assuming well-mixed system
Volume = LBox**3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; # everything will be in s^(-1)
SpontaneousRxnRates=[kplusDimer*ConversionFactor, kminusDimer, kplusTrimer*ConversionFactor, kminusTrimer, \
    kplusBarbed*ConversionFactor, kminusBarbed, kplusPointed*ConversionFactor, kminusPointed];

Nmon = int(Conc*Volume/uMInvToMicron3);
NmonSeeds = int(SeedConc*Volume/uMInvToMicron3);
print('Number of monomers %d' %Nmon)
NArp23 = int(ConcArp*Volume/uMInvToMicron3);
print('Number of Arp 2/3 %d' %NArp23)

Lens=np.array([LBox,LBox,LBox]);
seed = 1;

nThr=1;
AllActin = ActinMixedNucleates(Nmon,Lens,SpontaneousRxnRates,a,spacing,kbT,mu, seed,nThr);
if (ConcFormin > 0):
    NFormin = int(ConcFormin*Volume/uMInvToMicron3);
    print('Number of formins %d' %NFormin)
    NBarbed = [NFormin];
    BarbedOnOff = [kplusFor*ConversionFactor, kminusFor];
    AlphaDimersMinus = [1,0];
    AlphaTrimersMinus = [1,0];
    AlphaBarbedMinus = [1, 1];
    AllActin.InitializeBarbedBinders(NBarbed,BarbedOnOff,AlphaDimersMinus,AlphaTrimersMinus,AlphaBarbedMinus);
    AlphaDimersPlus = [1,kForNuc*ConversionFactor/kplusDimer];
    AlphaTrimersPlus = [1, (ForminAlphaNoProf*kplusBarbed+kplusPointed)/kplusTrimer]
    AlphaBarbedPlus = [1, ForminAlphaNoProf];
    AllActin.InitializeRateMatrices(AlphaDimersPlus,AlphaTrimersPlus,AlphaBarbedPlus);
if (ConcProf > 0):
    NProf = int(ConcProf*Volume/uMInvToMicron3);
    print('Number profilin %d' %NProf)
    NMonProts = [NProf]
    KsMon = [ProfEq*ConversionFactor];
    AlphasPointed = [1, AlphaPtdProf];
    AllActin.InitializeMonomerBinders(NMonProts,KsMon,AlphasPointed);
    # Re-initialize the rate matrices with all rates
    if (ConcFormin == 0):
        # Rates with profilin
        AlphaDimersPlus = [1,0];
        AlphaTrimersPlus = [1,0]
        AlphaBarbedPlus = [1,AlphaWithProf];
    else:
        AlphaDimersPlus = [1,kForNuc*ConversionFactor/kplusDimer, 0,0];
        AlphaTrimersPlus = [1, (ForminAlphaNoProf*kplusBarbed+kplusPointed)/kplusTrimer,\
                            0, (ForminAlphaWithProf*kplusBarbed+kplusPointed)/kplusTrimer];
        AlphaBarbedPlus = [1, ForminAlphaNoProf, AlphaWithProf, ForminAlphaWithProf];
    AllActin.InitializeRateMatrices(AlphaDimersPlus,AlphaTrimersPlus,AlphaBarbedPlus);
if (ConcArp > 0):
    NArp23 = int(ConcArp*Volume/uMInvToMicron3);
    BranchRates = [kplusARF*ConversionFactor**2, kMinusARF];
    AlphaBranch = [1];
    if (ConcProf > 0):
        AlphaBranch = [1,0];
    AllActin.InitializeBranchers(NArp23,BranchRates,AlphaBranch);


Tf = 28800;
dt = 10;
nSteps = int(Tf/dt+1e-6);

NumFibers = np.zeros(nSteps,dtype=np.int64);
FreeMonomers = np.zeros(nSteps,dtype=np.int64);
NumberPerFiber = np.array([],dtype=np.int64);
BranchedOrLinear = np.array([],dtype=bool);
BoundProteins = np.array([],dtype=bool);
# Initialize seeds (HAS TO BE DONE AFTER INITIALIZING RATES!!!)
AllActin.InitSeeds(NmonSeeds,100);

for i in range(nSteps):
    AllActin.React(dt);
    NumOnEach = AllActin.NumMonOnEachFiber();
    # Eliminate dimers and trimers (these don't count as actual fibers)
    NumDiTri = sum(NumOnEach<=3);
    NumOnEach = NumOnEach[NumDiTri:];
    NumberPerFiber = np.append(NumberPerFiber,NumOnEach)
    NumFibers[i] = len(NumOnEach);
    FreeMonomers[i] = AllActin.nFreeMonomers();
    BranchedState = AllActin.BranchedOrLinear(False); # 2 for mother, 1 for in branched structure, 0 for linear
    BranchedState=BranchedState[NumDiTri:];
    BranchedOrLinear = np.append(BranchedOrLinear,BranchedState)
    BoundProteins = np.append(BoundProteins,AllActin.BoundBarbedStates())
    print('Time %f, Percent free %f, number fibs %d' %((i+1)*dt, FreeMonomers[i]/Nmon, NumFibers[i]))

#FileName = 'Tf'+str(Tf)+'_Box'+str(LBox)+'_Actin'+str(Conc)+'uM_Seed_'+str(SeedConc)+'_KProf'+str(ProfEq)+'_Prof'+ \
    str(ConcProf)+'uM_Arp'+str(int(ConcArp*1000))+'nM_Formin'+str(int(ConcFormin*1e4)) \
    +'em4uM_'+str(seed)+'.txt';
# Write the output
FileName = 'Tf'+str(Tf)+'_Box'+str(LBox)+'_Actin'+str(Conc)+'uM_Seed'+str(SeedConc)+'uM_Arp'+str(int(ConcArp*1000))+'nM.txt';
# FreeMonConc is a nT array that gives you the concentration of free monomers in uM (used to measure when polymerization is complete)
np.savetxt('FreeMonConc'+FileName,FreeMonomers/Nmon*Conc);
# NumFibs is a nT array that gives you the number of fibers at each time step (use it to parse the next two arrays)
np.savetxt('NumFibs'+FileName,NumFibers);
# NumberPerFiber is a #fibers array that gives you the number of monomers on each fiber. Here I have already removed the number of monomers bound to dimers/trimers (because of that this number won't add exactly to Nmon-FreeMon)
np.savetxt('NumberPerFiber'+FileName,NumberPerFiber);
#np.savetxt('BoundProteins'+FileName,BoundProteins);
# BranchedOrLinear is a #fibers array that gives you the branching state of each filament. 2 = mother filament. 1 = daughter filament on a branched structure. 0 = linear filament. 
np.savetxt('BranchedOrLinear'+FileName,BranchedOrLinear);
