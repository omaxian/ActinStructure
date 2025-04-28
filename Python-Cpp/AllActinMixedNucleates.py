import numpy as np
import sys
from ActinMixedNucleates import ActinMixedNucleates

WritePos = False;

for ConcProf in [0,0.75,1.5,6]:
    for ConcArp in [80e-3, 160e-3]:
        for SeedConc in [0,0.001,0.002]:
            for ConcFormin in [0]: # in uM
                # Parameters
                Conc = 1.5;
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
                seed = int(sys.argv[1]);

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
                    NumberPerFiber = np.append(NumberPerFiber,NumOnEach)
                    NumFibers[i] = len(NumOnEach);
                    nFibFour=np.sum(NumOnEach>3);
                    FreeMonomers[i] = AllActin.nFreeMonomers();
                    BranchedOrLinear = np.append(BranchedOrLinear,AllActin.BranchedOrLinear(False))
                    BoundProteins = np.append(BoundProteins,AllActin.BoundBarbedStates())
                    print('Time %f, Percent free %f, number fibs %d' %((i+1)*dt, FreeMonomers[i]/Nmon, nFibFour))
                    if (i==0 and WritePos):
                        AllX0 = AllActin.AllX0();
                        AllTau = AllActin.AllTaus();
                    elif (WritePos):
                        AllX0 = np.append(AllX0,AllActin.AllX0(),axis=0);
                        AllTau = np.append(AllTau,AllActin.AllTaus(),axis=0);
                
                FileName = 'Tf'+str(Tf)+'_Box'+str(LBox)+'_Actin'+str(Conc)+'uM_Seed_'+str(SeedConc)+'_KProf'+str(ProfEq)+'_Prof'+ \
                    str(ConcProf)+'uM_Arp'+str(int(ConcArp*1000))+'nM_Formin'+str(int(ConcFormin*1e4)) \
                    +'em4uM_'+str(seed)+'.txt';
                np.savetxt('FreeMons'+FileName,FreeMonomers);    
                np.savetxt('NumFibs'+FileName,NumFibers);
                np.savetxt('StructInfo'+FileName,NumberPerFiber);
                np.savetxt('BoundProteins'+FileName,BoundProteins);
                np.savetxt('BranchedOrLinear'+FileName,BranchedOrLinear);
                if (WritePos):
                    np.savetxt('AllX0'+str(Conc)+FileName,AllX0);
                    np.savetxt('AllTau'+str(Conc)+FileName,AllTau);
