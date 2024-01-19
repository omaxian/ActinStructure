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
 # in uM

# Parameters from Kovar & Pollard paper for actin alone
# RETURN TO THE ACTUAL PARAMS LATER!
kplusDimer = 0; # uM^(-1)*s^(-1) 
kminusDimer = 0; #s^(-1)
kplusTrimer = 0; # uM^(-1)*s^(-1) 
kminusTrimer = 0; #s^(-1)
kplusBarbed = 0; # uM^(-1)*s^(-1) 
kminusBarbed = 1.4; #s^(-1)
kplusPointed = 0; #uM^(-1)*s^(-1)
kminusPointed = 0.8; #s^(-1)

# Formin rates
kForNuc = 2e-3; # uM^(-2)*s^(-1)
kplusFor = 29.1; # uM^(-1)*s^(-1)
kminusFor = 8.1e-1; # s^(-1)
ForminEnhance = 2;

# Arp 2/3 rates
kplusARF = 0;
kMinusARF = 0.5#*3.4e-3;

NPerBranch = np.array([4,6,8,10],dtype=np.int64);
Mothers = np.array([0,0,0,0],dtype=np.int64);
AttachPts = np.array([0,3,3,3],dtype=np.int64);
BoundBarbed = np.zeros(len(NPerBranch),dtype=np.int64);

# Convert to microscopic assuming well-mixed system
Volume = LBox**3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; # everything will be in s^(-1)
RxnRates=[kplusDimer*ConversionFactor, kminusDimer, kplusTrimer*ConversionFactor, kminusTrimer, \
    kplusBarbed*ConversionFactor, kminusBarbed, kplusPointed*ConversionFactor, kminusPointed];
RxnRatesFormin = [kForNuc*ConversionFactor**2, kplusFor*ConversionFactor, kminusFor, ForminEnhance];
RxnRatesArp23 = [kplusARF*ConversionFactor**2, kMinusARF];

Nmon =  np.sum(NPerBranch);
print('Number of monomers %d' %Nmon)
NArp23 = len(NPerBranch)-1;
print('Number of Arp 2/3 %d' %NArp23)

Lens=np.array([LBox,LBox,LBox]);
nError=10;
BrkTime=np.zeros(nError);
nTrial=1000;
TrialPerSeed=nTrial//nError
for er in range(nError):
    for seed in range(TrialPerSeed):
        nInFibers = [];
        seed+=er*TrialPerSeed;
        nThr=1;
        AllActin = ActinMixedNucleates(Nmon,Lens,RxnRates,a,spacing,kbT,mu, seed,nThr);
        NumOnEach = AllActin.NumMonOnEachFiber();
        AllActin.InitializeBranchers(NArp23,RxnRatesArp23);
        AllActin.InitBranchedStructure(NPerBranch,Mothers,AttachPts,BoundBarbed);

        Tf = 100;
        dt = 0.1;
        nSteps = int(Tf/dt+1e-6);
        NumStructures = np.zeros(nSteps+1,dtype=np.int64);

        AllX = AllActin.getX();
        np.savetxt('AllX.txt',AllX);
        NumberPerStruct = AllActin.NumMonOnEachFiber();
        BoundProteins = AllActin.BoundBarbedStates();
        nStruct=len(NPerBranch);
        NumStructures[0]=len(NPerBranch);
        Intact = True
        i=0;

        while (Intact):
            AllActin.React(dt);
            NumOnEach = AllActin.NumMonOnEachFiber();
            NumberPerStruct = np.append(NumberPerStruct,NumOnEach)
            BoundProteins = np.append(BoundProteins,AllActin.BoundBarbedStates())
            try:
                NumStructures[i+1]=len(NumOnEach)-3;
            except:
                print('No depoly - exiting')
                np.savetxt('AllX'+str(seed)+'.txt',AllX);
                np.savetxt('NumStruct'+str(seed)+'.txt',NumStructures);
                np.savetxt('StructInfo'+str(seed)+'.txt',NumberPerStruct);
                np.savetxt('BoundProteins'+str(seed)+'.txt',BoundProteins);
                import sys
                sys.exit()
            nStruct=len(NumOnEach)-3
            try:
                Intact = nStruct > 1 or NumOnEach[3] > 4;
            except:
                Intact = False;
            i+=1;
            #print(NumOnEach)
            #if (not Intact):
            #    print(NumOnEach)
            NmonNow = NumOnEach[0]+2*NumOnEach[1]+3*NumOnEach[2]+np.sum(NumOnEach[3:]);
            if (NmonNow != np.sum(NPerBranch)):
                print(NumOnEach)
                import sys
                sys.exit()
            AllX = np.append(AllX,AllActin.getX(),axis=0);
            #print('Time %f, Percent free %f' %((i+1)*dt, NumOnEach[0]/Nmon))
        BrkTime[er]+= i*dt/TrialPerSeed;
        #print('Break up time %f' %(i*dt))
        
        del AllActin;
    print('Mean break up time %f' %BrkTime[er])
print(BrkTime)
print(np.mean(BrkTime))
print(2*np.std(BrkTime)/np.sqrt(nError))
np.savetxt('AllX'+str(seed)+'.txt',AllX);
np.savetxt('NumStruct'+str(seed)+'.txt',NumStructures);
np.savetxt('StructInfo'+str(seed)+'.txt',NumberPerStruct);
np.savetxt('BoundProteins'+str(seed)+'.txt',BoundProteins);


#np.savetxt('NumFibs'+str(Conc)+'uM_Formin'+str((ConcFormin*1000))+'nM_Alpha'+str(ForminEnhance)+'_'+str(seed)+'.txt',NumFibers);
#np.savetxt('StructInfo'+str(Conc)+'uM_Formin'+str((ConcFormin*1000))+'nM_Alpha'+str(ForminEnhance)+'_'+str(seed)+'.txt',StructInfo);
#np.savetxt('BoundFormin'+str(Conc)+'uM_Formin'+str((ConcFormin*1000))+'nM_Alpha'+str(ForminEnhance)+'_'+str(seed)+'.txt',BoundFormins);

