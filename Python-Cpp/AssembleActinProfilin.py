import numpy as np
import sys
from ActinMixedNucleates import ActinMixedNucleates

# Parameters for actin and profilin which you can change
ActinConc = 1.5; # Actin concentration in uM
ProfConc = 1.5; # Profilin concentration in uM
SeedConc = 0; # concentration of seeds 
LBox = 10; # size of domain being simulated in um (one side)
ProfEq = 1; # Profilin equilibrium constant in uM^(-1)
AlphaDimerProf = 0; # Relative rate that profilin-actin forms dimers
AlphaTrimerProf = 0; # Relative rate that profilin-actin adds to dimers to form trimers
AlphaBarbedProf = 1; # Relative ratio for profilin adding to the barbed end 
AlphaPtdProf = 0; # Relative ratio for profilin adding to the pointed end 
Tf = 28800;     # Simulation time in s
dt = 10;        # Save interval

               
# Parameters from Kovar & Pollard paper for actin alone
kplusDimer = 3.5e-6; # uM^(-1)*s^(-1) 
kminusDimer = 0.041; #s^(-1)
kplusTrimer = 13e-5; # uM^(-1)*s^(-1) 
kminusTrimer = 22; #s^(-1)
kplusBarbed = 11.6; # uM^(-1)*s^(-1) 
kminusBarbed = 1.4; #s^(-1)
kplusPointed = 1.3; #uM^(-1)*s^(-1)
kminusPointed = 0.8; #s^(-1)

# Parameters for diffusion (these aren't involved in the code you're running)
a = 4e-3;
kbT = 4.1e-3;
spacing = 0.5; # units of a
mu = 0.01;
                
# Convert to microscopic assuming well-mixed system
Volume = LBox**3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; # everything will be in s^(-1)
SpontaneousRxnRates=[kplusDimer*ConversionFactor, kminusDimer, kplusTrimer*ConversionFactor, kminusTrimer, \
    kplusBarbed*ConversionFactor, kminusBarbed, kplusPointed*ConversionFactor, kminusPointed];
Nmon = int(ActinConc/ConversionFactor);
NmonSeeds = int(SeedConc/ConversionFactor);
print('Number of monomers %d' %Nmon)
Lens=np.array([LBox,LBox,LBox]);
seed = int(sys.argv[1]);

nThr=1;
AllActin = ActinMixedNucleates(Nmon,Lens,SpontaneousRxnRates,a,spacing,kbT,mu, seed,nThr);
               
if (ProfConc > 0):
    NProf = int(ProfConc*Volume/uMInvToMicron3);
    print('Number profilin %d' %NProf)
    NMonProts = [NProf]
    KsMon = [ProfEq*ConversionFactor];
    AlphasPointed = [1, AlphaPtdProf];
    AllActin.InitializeMonomerBinders(NMonProts,KsMon,AlphasPointed);
    # Initialize rates with profilin
    AlphaDimersPlus = [1,AlphaDimerProf];
    AlphaTrimersPlus = [1,AlphaTrimerProf]
    AlphaBarbedPlus = [1,AlphaBarbedProf];
    AllActin.InitializeRateMatrices(AlphaDimersPlus,AlphaTrimersPlus,AlphaBarbedPlus);
    
nSteps = int(Tf/dt+1e-6);
NumFibers = np.zeros(nSteps+1,dtype=np.int64);
FreeMonomers = np.zeros(nSteps+1,dtype=np.int64);
ProfBdMonomers = np.zeros(nSteps+1,dtype=np.int64);
NumberPerFiber = np.array([],dtype=np.int64);
# Initialize seeds (HAS TO BE DONE AFTER INITIALIZING RATES!!!)
nMonPerSeed=100;
AllActin.InitSeeds(NmonSeeds,nMonPerSeed);
FreeMonomers[0] = Nmon-NmonSeeds;

# Do the simulation
for i in range(nSteps):
    AllActin.React(dt);
    FreeAndProfBd = AllActin.FreeAndProfBoundMonomers();
    FreeMonomers[i+1] = FreeAndProfBd[0];
    ProfBdMonomers[i+1] = FreeAndProfBd[1];
    NumOnEach = AllActin.NumMonOnEachFiber();
    NumFibers[i+1] = len(NumOnEach);
    NumberPerFiber = np.append(NumberPerFiber,NumOnEach)
    nFibFour=np.sum(NumOnEach>3);
    print('Time %f, Percent free %f, number fibs %d' %((i+1)*dt, FreeMonomers[i]/Nmon, nFibFour))
                
FileName = 'Actin'+str(ActinConc)+'uM_KProf'+str(ProfEq)+'_Prof'+ \
    str(ProfConc)+'uM_'+str(seed)+'.txt';
np.savetxt('FreeMonConc'+FileName,FreeMonomers);    
np.savetxt('ProfBoundConc'+FileName,ProfBdMonomers);   
np.savetxt('NumFibs'+FileName,NumFibers);
np.savetxt('NumPerFib'+FileName,NumberPerFiber);

