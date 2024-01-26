#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ActinStructure.cpp"

namespace py=pybind11;
typedef py::array_t<int, py::array::c_style | py::array::forcecast> npInt; 
typedef py::array_t<double, py::array::c_style | py::array::forcecast> npDoub; 

class ActinMixedNucleates {
    
    public:
          
    ActinMixedNucleates(const int &NMon, const vec3 &Lengths, const vec &SpontaneousRxnRates,
         double a, double spacingInTermsOfA, double kbT, double mu, int seed, int nThr){
        /*
        Initialization
        Nmon = total number of monomers, including those attached to fibers
        NperFiber = integer vector of number of monomers on each of the nascent fibers
        Lengths = 3-vector of (x,y,z) periodic domain lengths
        a = hydrodynamic radius
        kbT = thermal energy
        mu = fluid viscosity
        RxnRates = vector with the rates (see below for ordering)
        */
        
        // Initialize random number generators
        for (int d=0; d<3;d++){
            _Lens[d]=Lengths[d];
        }
        
        // Diffusion variables
        _a = a;
        _kbT = kbT;
        _mu = mu;
        _spacing = spacingInTermsOfA*a;
        
        // Reaction variables
        // We assume all reaction rates are given in units of s^(-1)
        _DimerPlus = vec(1,SpontaneousRxnRates[0]);
        _DimerMinus = vec(1,SpontaneousRxnRates[1]);
        _TrimerPlus = vec(1,SpontaneousRxnRates[2]);
        _TrimerMinus = vec(1,SpontaneousRxnRates[3]);
        
        // Fiber binding rates
        _BarbedPlus = vec(1,SpontaneousRxnRates[4]); // Barbed binding
        _BarbedMinus = vec(1,SpontaneousRxnRates[5]); // Barbed unbinding
        _PointedPlus = vec(1,SpontaneousRxnRates[6]); // Pointed binding
        _PointedMinus = SpontaneousRxnRates[7]; // Pointed unbinding
        
        // Initialize variables relating to barbed-end-bound and monomer bound proteins
        _nFreeMon = NMon;
        _nDimers = intvec(1,0);
        _nTrimers = intvec(1,0);
        _nBarbedProts = 0;
        _nBarbedBinders = intvec(0);
        _BarbedBindersRates = vec(0);
        _FreeBranchers = 0;
        _BranchOnRates = vec(0);
        _BranchOffRate = 0;
        _nMonProts = 0;
        _nMonBinders = intvec(0);
        _MonBindKEqs = vec(0);                
        
        rngu.seed(seed);
        unifdist = std::uniform_real_distribution<double>(0.0,1.0);
        rng.seed(seed);
        normaldist = std::normal_distribution<double>(0.0,1.0);
    }
    
    void InitializeMonomerBinders(const intvec &nMonBinders, const vec &EquilConsts, const vec &PointedAlphas){
        _nMonProts = nMonBinders.size();
        _nMonBinders = intvec(_nMonProts);
        std::memcpy(_nMonBinders.data(),nMonBinders.data(),nMonBinders.size()*sizeof(int));
        _MonBindKEqs = vec(_nMonProts);
        std::memcpy(_MonBindKEqs.data(),EquilConsts.data(),EquilConsts.size()*sizeof(double));
        _PointedPlus.resize(_nMonProts+1);
        int pAlphSize = PointedAlphas.size();
        if (_nMonProts+1 != pAlphSize){
            throw std::runtime_error("Size mismatch in pointed alpha array");
        }
        for (int iProt=1; iProt <= _nMonProts; iProt++){
            _PointedPlus[iProt] = _PointedPlus[0]*PointedAlphas[iProt];
        }
    } 

    void InitializeBarbedBinders(const intvec &nBarbedBinders,const vec &BarbedBinderRates, 
        const vec &AlphaDimer, const vec &AlphaTrimer, const vec &AlphaBarbedMinus){
        
        _nBarbedProts = nBarbedBinders.size();
        _nBarbedBinders = intvec(_nBarbedProts);
        std::memcpy(_nBarbedBinders.data(),nBarbedBinders.data(),nBarbedBinders.size()*sizeof(int));
        _BarbedBindersRates = vec(2*_nBarbedProts);
        std::memcpy(_BarbedBindersRates.data(),BarbedBinderRates.data(),BarbedBinderRates.size()*sizeof(double));
        
        _DimerMinus.resize(_nBarbedProts+1);
        _TrimerMinus.resize(_nBarbedProts+1);
        _BarbedMinus.resize(_nBarbedProts+1);
        int nDimSize = AlphaDimer.size();
        if (nDimSize != _nBarbedProts+1){
            throw std::runtime_error("Size mismatch in barbed alpha array (must have leading 1)");
        }
        for (int iProt=1; iProt <= _nBarbedProts; iProt++){
            _DimerMinus[iProt] = _DimerMinus[0]*AlphaDimer[iProt];
            _TrimerMinus[iProt]= _TrimerMinus[0]*AlphaTrimer[iProt];
            _BarbedMinus[iProt] = _BarbedMinus[0]*AlphaBarbedMinus[iProt];
        }

        _nDimers.resize(_nBarbedProts+1);
        _nTrimers.resize(_nBarbedProts+1);
         // The rates with proteins attached to monomers or at barbed ends
    }
    
    void InitializeRateMatrices(const vec &DimerAlphas, const vec &TrimerAlphas, const vec &BarbedEndAlphas){
        // Matrices of rates (i=monomer bound protein, j = barbed bound protein)
        int nBarbedDim = _nBarbedProts+1;
        int nMonDim = _nMonProts+1;
        int ExpectedDim = nBarbedDim*nMonDim;
        int nDimSize = DimerAlphas.size();
        int nTriSize = TrimerAlphas.size();
        int nBAlphSize = BarbedEndAlphas.size();
        if (nDimSize != ExpectedDim || nTriSize != ExpectedDim || nBAlphSize != ExpectedDim){
            throw std::runtime_error("Size mismatch in matrix of alphas (must have leading 1)");
        }
        _DimerPlus.resize(ExpectedDim);
        _TrimerPlus.resize(ExpectedDim);
        _BarbedPlus.resize(ExpectedDim);
        // Each row is a monomer bound protein, each column is a barbed protein
        for (int iMB=0; iMB < nMonDim; iMB++){
            for (int iBB=0; iBB < nBarbedDim; iBB++){
                int index = MatrixIndex(iMB,iBB);
                if (index > 0){
                    _DimerPlus[index]=_DimerPlus[0]*DimerAlphas[index];
                    _TrimerPlus[index]=_TrimerPlus[0]*TrimerAlphas[index];
                    _BarbedPlus[index]=_BarbedPlus[0]*BarbedEndAlphas[index];
                }
            }
        }
        for (uint i =0; i <  _DimerPlus.size(); i++){
            std::cout << "Dimer plus = " << _DimerPlus[i] << std::endl;
        }
    }
        
   
    void InitializeBranchers(int nBranchers, const vec &BranchRates, const vec &AlphaBranchOn){
        _FreeBranchers = nBranchers;
        _BranchOffRate = BranchRates[1];
        // The on rate should be in units of 1/(monomers on the mother * time)
        _BranchOnRates = vec(AlphaBranchOn.size());
        int nDimSize = AlphaBranchOn.size();
        if (nDimSize != _nMonProts+1){
            throw std::runtime_error("Size mismatch in matrix of branching alphas (must have leading 1)");
        }
        for (int iMB=0; iMB < _nMonProts+1; iMB++){
            _BranchOnRates[iMB] = AlphaBranchOn[iMB]*BranchRates[0];
        }
    }
    
    void InitBranchedStructure(const intvec &nMonsPerBranch, const intvec &Mothers,const intvec AttachPts, const intvec &ProteinsOnBarbed){
        /* 
        Initialize a branched filament (for debugging mostly). The intvector 
        nMonsPerBranch gives the number of monomers on each branch, and the 
        vector ProteinsOnBarbed tells us if there is a formin on each branch 
        */
        int nMonsIn = std::accumulate(nMonsPerBranch.begin(), nMonsPerBranch.end(),0);
        _nFreeMon-=nMonsIn;
        uint nBranches = nMonsPerBranch.size()-1;
        _FreeBranchers-=nBranches;
        // Initialize fiber
        _Fibers.push_back(std::make_shared<Fiber>(_Lens, nMonsPerBranch[0], _nMonProts, _nBarbedProts, 
                    _BarbedBindersRates, _PointedPlus, _BarbedPlus, _PointedMinus, _BarbedMinus, _BranchOnRates,_BranchOffRate,
                    _spacing, _a, _mu, _kbT,-1));
        // Branches
        uint ind = _Fibers.size()-1;
        for (uint iFib=0; iFib <= nBranches; iFib++){
            if (iFib==1){
                _Fibers[ind] =  std::make_shared<BranchedFiber>(_Fibers[ind].get());
            } else if (iFib > 1) {
                _Fibers[ind]->addBranch(Mothers[iFib]);
            } 
            _Fibers[ind]->SetNumMons(iFib,nMonsPerBranch[iFib]);
            _Fibers[ind]->SetAttachPt(iFib,AttachPts[iFib]);
            _Fibers[ind]->SetBoundBarbed(iFib,ProteinsOnBarbed[iFib]);
            if (ProteinsOnBarbed[iFib] > 0){
                int iProt = ProteinsOnBarbed[iFib]-1;
                _nBarbedBinders[iProt]--;
            }
        }
    }
        
    
    void Diffuse(double dt){
        // Diffuse the fibers
        for(auto&& FibObj: _Fibers){
            FibObj->Diffuse(dt);
        }
    }
    
    intvec NumBoundToEachProtein(int nFreeMon){
        intvec NToEach(_nMonProts+1,nFreeMon);
        if (_nMonProts==0){
            return NToEach;
        }
        // Now we are in the case when there are multiple proteins
        double denom=1;
        for (int j=0; j < _nMonProts; j++){
            denom+=_MonBindKEqs[j]*_nMonBinders[j];
        }
        int TotalBound=0;
        for (int j=1; j <= _nMonProts; j++){
            NToEach[j]=1.0/denom*nFreeMon*_nMonBinders[j-1]*_MonBindKEqs[j-1];
            TotalBound+=NToEach[j];
        }
        NToEach[0] = nFreeMon-TotalBound;
        return NToEach;
    }
    
    void React(double dt){
        /*
        Event-driven simulation with well-mixed monomers and nucleates.
        */       
        double t = 0;
        while (t < dt){
            int index = 0;
            //std::cout << "System time " << t << std::endl;
            // Delta t is the minimum time, or time for next reaction. It gets modified 
            // in each method
            intvec NumOnEach = NumBoundToEachProtein(_nFreeMon);
            double deltaT = TimeNucleationReactions(index,NumOnEach);
            deltaT = TimeNucleateBarbedBindReactions(index, deltaT);
            int numNucRxns = 5*(_nBarbedProts+1);
            int numBindNucRxns = 4*_nBarbedProts;
            //std::cout << "Nuc time and index " << deltaT << " , " << index << std::endl;
            uint nFib = _Fibers.size();
            //std::cout << "Number of fibers " << nFib << std::endl;
            if (nFib > 0){
                vec FiberEventTimes(nFib);
                // TODO: make this loop parallel
                for (uint iFib = 0; iFib < nFib; iFib++){
                    // Rxn 0: add to pointed end
                    FiberEventTimes[iFib] =  _Fibers[iFib]->NextEventTime(NumOnEach, _nBarbedBinders, _FreeBranchers);
                    //std::cout << "Event time " << FiberEventTimes[iFib] << std::endl;
                }
                // Find the minimum 
                auto it = std::min_element(std::begin(FiberEventTimes), std::end(FiberEventTimes));
                int MinFibIndex = std::distance(std::begin(FiberEventTimes), it);
                if (FiberEventTimes[MinFibIndex] < deltaT){
                    deltaT = FiberEventTimes[MinFibIndex];
                    index = numNucRxns+numBindNucRxns+MinFibIndex;
                }
            }
            //
            if (t+deltaT > dt || deltaT < 0){
                return;
            }
            if (index < numNucRxns){
                ProcessNucleateReaction(index);
            } else if (index < numNucRxns+numBindNucRxns){
                ProcessNucleateBarbedBindReaction(index);
            } else {
                int iFib = index - (numNucRxns+numBindNucRxns);
                int StructChange = _Fibers[iFib]->NeedsStructureChange();
                if (StructChange==0){
                    _Fibers[iFib]->ReactNextEvent(_nFreeMon,_nBarbedBinders,_FreeBranchers);
                } else if (StructChange == -1){
                    // Delete structure or turn branched fiber into linear fiber
                    uint nFibers = _Fibers[iFib]->nFibers();
                    int BarbedProtOnDepol = _Fibers[iFib]->getBoundBarbed(nFibers-1);
                    if (nFibers == 1){
                          // Find out what's on the end
                        _Fibers.erase(_Fibers.begin() + iFib);   
                        _nTrimers[BarbedProtOnDepol]++;
                        _nFreeMon++;
                    } else {
                        // Convert branched filament to linear 
                        vec3 X0, tau;
                        int nMonomers = _Fibers[iFib]->PrepareForLinearSwitch(X0,tau);
                        int BarbedOnRemain = _Fibers[iFib]->getBoundBarbed(0);
                        if (nMonomers < 4){
                            std::cout << "Fiber should have more than 4 monomers! " << std::endl;
                        }
                        _Fibers[iFib] =  std::make_shared<Fiber>(_Lens, nMonomers, _nMonProts, _nBarbedProts, 
                            _BarbedBindersRates, _PointedPlus, _BarbedPlus, _PointedMinus, _BarbedMinus, _BranchOnRates,_BranchOffRate,
                            _spacing, _a, _mu, _kbT,-1);
                        _Fibers[iFib]-> SetX0AndTau(X0,tau);
                        _Fibers[iFib]-> SetBoundBarbed(0,BarbedOnRemain);
                        _nFreeMon++;
                        _FreeBranchers++;
                        if (BarbedProtOnDepol > 0){
                            _nBarbedBinders[BarbedProtOnDepol-1]++;
                        }
                    }
                } else if (StructChange==1){ // Converting linear filament to branched  
                    _Fibers[iFib] =  std::make_shared<BranchedFiber>(_Fibers[iFib].get());
                    _nFreeMon--;
                    _FreeBranchers--;
                } 
            }
            // Check conservation of monomers
            t+=deltaT;
        } // end while loop       
    }
    
    npDoub getX(){
        vec AllX(0);
        for(auto&& FibObj: _Fibers){
            vec XFib = FibObj->getX();
            AllX.insert(AllX.end(), XFib.begin(), XFib.end());
        }
        return makePyDoubleArray(AllX);
    }
    
    npDoub AllX0(){
        vec AllX0(0);
        for(auto&& FibObj: _Fibers){
            vec X0Fib = FibObj->getX0();
            AllX0.insert(AllX0.end(), X0Fib.begin(), X0Fib.end());
        }
        return makePyDoubleArray(AllX0);
    }
    
    npDoub AllTaus(){
        vec AllTau(0);
        for(auto&& FibObj: _Fibers){
            vec TauFib = FibObj->getTau();
            AllTau.insert(AllTau.end(), TauFib.begin(), TauFib.end());
        }
        return makePyDoubleArray(AllTau);
    }
    
    int nFreeMonomers(){
        return _nFreeMon;
    }
    
    npInt NumMonOnEachFiber(){
        intvec nPerFiber(0);
        // Dimers
        for (int iB=0; iB <= _nBarbedProts; iB++){
            for (int iD = 0; iD < _nDimers[iB]; iD++){
                nPerFiber.push_back(2);
            }
        }
        // Trimers
        for (int iB=0; iB <= _nBarbedProts; iB++){
            for (int iD = 0; iD < _nTrimers[iB]; iD++){
                nPerFiber.push_back(3);
            }
        }
        // Fibers
        for(auto&& FibObj: _Fibers){
            int nFibThis = FibObj->nFibers();
            for (int jFib=0; jFib < nFibThis; jFib++){
                nPerFiber.push_back(FibObj->NumMonomers(jFib));
            }
        }    
        return makePyArray(nPerFiber);
    }
    
    npInt BoundBarbedStates(){
        intvec BoundStates(0);
        // Dimers
        for (int iB=0; iB <= _nBarbedProts; iB++){
            for (int iD = 0; iD < _nDimers[iB]; iD++){
                BoundStates.push_back(iB);
            }
        }
        // Trimers
        for (int iB=0; iB <= _nBarbedProts; iB++){
            for (int iD = 0; iD < _nTrimers[iB]; iD++){
                BoundStates.push_back(iB);
            }
        }
        // Fibers
        for(auto&& FibObj: _Fibers){
            int nFibThis = FibObj->nFibers();
            for (int jFib=0; jFib < nFibThis; jFib++){
                BoundStates.push_back(FibObj->getBoundBarbed(jFib));
            }
        }    
        return makePyArray(BoundStates);
    }
    
    npInt BranchedOrLinear(bool MothersAreBranched){
        intvec BranchOrLin(0);
        // Dimers
        for (int iB=0; iB <= _nBarbedProts; iB++){
            for (int iD = 0; iD < _nDimers[iB]; iD++){
                BranchOrLin.push_back(0);
            }
        }
        // Trimers
        for (int iB=0; iB <= _nBarbedProts; iB++){
            for (int iD = 0; iD < _nTrimers[iB]; iD++){
                BranchOrLin.push_back(0);
            }
        }
        for(auto&& FibObj: _Fibers){
            int nFibThis = FibObj->nFibers();
            for (int jFib=0; jFib < nFibThis; jFib++){
                if (jFib==0 && nFibThis > 1){
                    BranchOrLin.push_back(2); // Mother of branched
                } else if (nFibThis > 1){
                    BranchOrLin.push_back(1); // Branched, but not a mother
                } else {
                    BranchOrLin.push_back(0);
                }    
            }
        }    
        return makePyArray(BranchOrLin);
    }      
        
    private:
        // Base parameters 
        std::vector <std::shared_ptr<Fiber>> _Fibers;
        int _nFreeMon;
        intvec _nDimers, _nTrimers;
        double _a, _kbT, _mu, _spacing;
        vec3 _Lens;
        vec _X;  
        std::normal_distribution<double> normaldist;
        std::mt19937_64 rng;
        std::uniform_real_distribution<double> unifdist;
        std::mt19937_64 rngu;
        
        // The rates with proteins attached to monomers or at barbed ends
        int _nMonProts, _nBarbedProts;
        intvec _nMonBinders, _nBarbedBinders;
        vec _MonBindKEqs, _BarbedBindersRates; // rates to come on and off
        vec _DimerPlus, _TrimerPlus, _BarbedPlus; // matrices of rates (i=monomer bound protein, j = barbed bound protein)
        vec _DimerMinus, _TrimerMinus, _BarbedMinus; // vector of rates only depending on barbed protein
        vec _PointedPlus; // vector of rates only depending on monomer protein
        double _PointedMinus;
        
        // Branchers
        int _FreeBranchers;
        vec _BranchOnRates;
        double _BranchOffRate;
        
        int TotalMonomers(){
            int TotMon = _nFreeMon;
            TotMon+= 2*std::accumulate(_nDimers.begin(), _nDimers.end(),0);
            TotMon+= 3*std::accumulate(_nTrimers.begin(), _nTrimers.end(),0);
            for(auto&& FibObj: _Fibers){
                TotMon+= FibObj->TotalMonomers();
            }
            return TotMon;
        }
        
        uint nTotalFibers(){
            uint nFib = 0;
            for(auto&& FibObj: _Fibers){
                nFib+=FibObj->nFibers();
            }
            return nFib;     
        }    
               
        double TimeNucleationReactions(int &index, const intvec &NumBoundToEach){
            int NumPossBarbed = _nBarbedProts+1;
            int NumPossMon = _nMonProts+1;
            double deltaT = 1.0/0.0;
            for (int iBp = 0; iBp < NumPossBarbed; iBp++){
                // Reaction 0: formation of dimers
                for (int iMp = 0; iMp < NumPossMon; iMp++){
                    int matindex = MatrixIndex(iMp, iBp);
                    double RateDimerForm = _DimerPlus[matindex]*NumBoundToEach[iMp]*NumBoundToEach[iMp];
                    if (iBp > 0){
                        RateDimerForm*= _nBarbedBinders[iBp-1]; 
                    }
                    if (NumBoundToEach[iMp] < 2){
                        RateDimerForm = 0;
                    }    
                    double TryDeltaT = logrand()/RateDimerForm;
                    if (TryDeltaT < deltaT){
                        deltaT = TryDeltaT;
                        index = iBp;
                    }    
                }
                // Reaction 1: break-up of dimers
                double RateDimerBreakup = _DimerMinus[iBp]*_nDimers[iBp];
                double TryDeltaT = logrand()/RateDimerBreakup; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = NumPossBarbed+iBp;
                }
                // Reaction 2: trimer formation
                for (int iMp = 0; iMp < NumPossMon; iMp++){
                    int matindex = MatrixIndex(iMp, iBp);
                    double RateTrimerForm = _TrimerPlus[matindex]*NumBoundToEach[iMp]*_nDimers[iBp];
                    TryDeltaT = logrand()/RateTrimerForm; 
                    if (TryDeltaT < deltaT){
                        deltaT = TryDeltaT;
                        index = 2*NumPossBarbed+iBp;
                    }
                }
                // Reaction 3: trimer breakup
                double RateTrimerBreakup = _TrimerMinus[iBp]*_nTrimers[iBp];
                TryDeltaT = logrand()/RateTrimerBreakup; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = 3*NumPossBarbed+iBp;
                }
                // Reaction 4: tetramer formation
                for (int iMp = 0; iMp < NumPossMon; iMp++){
                    int matindex = MatrixIndex(iMp, iBp);
                    double TetOnRate = _BarbedPlus[matindex]+_PointedPlus[iMp];
                    double RateTetramerForm = TetOnRate*_nTrimers[iBp]*NumBoundToEach[iMp];
                    TryDeltaT = logrand()/RateTetramerForm; 
                    if (TryDeltaT < deltaT){
                        deltaT = TryDeltaT;
                        index = 4*NumPossBarbed+iBp;
                    }
                }
            }
            return deltaT;
        }
        
        double TimeNucleateBarbedBindReactions(int &index, double deltaT){
            if (_nBarbedProts==0){
                return deltaT;
            }
            int FirstIndex = 5*(_nBarbedProts+1);
            // There is something which can bind barbed ends
            for (int iBp = 0; iBp < _nBarbedProts; iBp++){ 
                // Binding empty dimers
                double RateBindToDimer = _BarbedBindersRates[2*iBp]*_nBarbedBinders[iBp]*_nDimers[0];
                double TryDeltaT = logrand()/RateBindToDimer; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = FirstIndex+iBp;
                }
                double RateUnbindFromDimer = _BarbedBindersRates[2*iBp+1]*_nDimers[iBp+1];
                TryDeltaT = logrand()/RateUnbindFromDimer; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = FirstIndex+_nBarbedProts+iBp;
                }
                double RateBindToTrimer = _BarbedBindersRates[2*iBp]*_nBarbedBinders[iBp]*_nTrimers[0];
                TryDeltaT = logrand()/RateBindToTrimer; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = FirstIndex+2*_nBarbedProts+iBp;
                }
                double RateUnbindTrimer = _BarbedBindersRates[2*iBp+1]*_nTrimers[iBp+1];
                TryDeltaT = logrand()/RateUnbindTrimer; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index =  FirstIndex+3*_nBarbedProts+iBp;
                }             
            }
            return deltaT;
        }
        
        void ProcessNucleateReaction(const int &index){
            int NumPossBarbed = _nBarbedProts+1;
            int ProtOnBarbed = index % NumPossBarbed;
            if (index < NumPossBarbed){
                // Dimer formation
                _nFreeMon-=2;
                _nDimers[ProtOnBarbed]++;
                if (ProtOnBarbed > 0){
                    _nBarbedBinders[ProtOnBarbed-1]--;
                }
            } else if (index < 2*NumPossBarbed){
                // Dimer break
                _nFreeMon+=2;
                _nDimers[ProtOnBarbed]--;
                if (ProtOnBarbed > 0){
                    _nBarbedBinders[ProtOnBarbed-1]++;
                }
            } else if (index < 3*NumPossBarbed){
                // Trimer formation
                _nFreeMon--;
                _nDimers[ProtOnBarbed]--;
                _nTrimers[ProtOnBarbed]++;
            } else if (index < 4*NumPossBarbed){
                // Trimer breakup
                _nFreeMon++;
                _nDimers[ProtOnBarbed]++;
                _nTrimers[ProtOnBarbed]--;
            } else if (index < 5*NumPossBarbed){
                // Forming tetramer - establish object for it
                    _nFreeMon--;
                    _nTrimers[ProtOnBarbed]--;
                    _Fibers.push_back(std::make_shared<Fiber>(_Lens, 4, _nMonProts, _nBarbedProts, 
                            _BarbedBindersRates, _PointedPlus, _BarbedPlus, _PointedMinus, _BarbedMinus, _BranchOnRates,_BranchOffRate,
                            _spacing, _a, _mu, _kbT,-1));
                    int iB = index- 4*NumPossBarbed;
                    _Fibers[_Fibers.size()-1]->SetBoundBarbed(0,iB);
            } else {
                std::cout << "Should not be here! " << std::endl;
            }
        }  
        
        void ProcessNucleateBarbedBindReaction(const int &index){
            int FirstIndex = 5*(_nBarbedProts+1);
            int RelIndex  = index-FirstIndex;
            int iBp = RelIndex % _nBarbedProts;
            // There is something which can bind barbed ends
            if (RelIndex < _nBarbedProts){
                // Binding to dimer
                _nDimers[0]--;
                _nDimers[iBp+1]++;
                _nBarbedBinders[iBp]--;
            } else if (RelIndex < 2*_nBarbedProts){
                // Unbinding from dimer
                _nDimers[0]++;
                _nDimers[iBp+1]--;
                _nBarbedBinders[iBp]++;
            } else if (RelIndex < 3*_nBarbedProts){
                // Binding to trimer
                _nTrimers[0]--;
                _nTrimers[iBp+1]++;
                _nBarbedBinders[iBp]--;
            } else { 
                // Unbinding from trimer
                _nTrimers[0]++;
                _nTrimers[iBp+1]--;
                _nBarbedBinders[iBp]++;
            }
        }
                        
        int MatrixIndex(int iMB, int iBB){
            return (_nBarbedProts+1)*iMB+iBB;     
        }  

        double logrand(){
            return -log(1.0-unifdist(rngu));
        }


        npDoub makePyDoubleArray(vec &cppvec){
            ssize_t              ndim    = 2;
            std::vector<ssize_t> shape   = { (long) cppvec.size()/3 , 3 };
            std::vector<ssize_t> strides = { sizeof(double)*3 , sizeof(double) };

            // return 2-D NumPy array
            return py::array(py::buffer_info(
            cppvec.data(),                       /* data as contiguous array  */
            sizeof(double),                          /* size of one scalar        */
            py::format_descriptor<double>::format(), /* data type                 */
            ndim,                                    /* number of dimensions      */
            shape,                                   /* shape of the matrix       */
            strides                                  /* strides for each axis     */
            ));
        } 
                
        npInt makePyArray(intvec &cppvec){
            // allocate py::array (to pass the result of the C++ function to Python)
            auto pyArray = py::array_t<int>(cppvec.size());
            auto result_buffer = pyArray.request();
            int *result_ptr    = (int *) result_buffer.ptr;
            // copy std::vector -> py::array
            std::memcpy(result_ptr,cppvec.data(),cppvec.size()*sizeof(int));
            return pyArray;
        }    
};

PYBIND11_MODULE(ActinMixedNucleates, m) {
    py::class_<ActinMixedNucleates>(m, "ActinMixedNucleates")
        .def(py::init<uint, vec3, vec, double, double, double, double,int,int>())
        .def("InitializeMonomerBinders",&ActinMixedNucleates::InitializeMonomerBinders)
        .def("InitializeBarbedBinders",&ActinMixedNucleates::InitializeBarbedBinders)
        .def("InitializeRateMatrices",&ActinMixedNucleates::InitializeRateMatrices)
        .def("InitializeBranchers",&ActinMixedNucleates::InitializeBranchers)
        .def("InitBranchedStructure",&ActinMixedNucleates::InitBranchedStructure)
        .def("Diffuse",&ActinMixedNucleates::Diffuse)
        .def("React",&ActinMixedNucleates::React)
        .def("getX", &ActinMixedNucleates::getX)
        .def("AllX0",&ActinMixedNucleates::AllX0)
        .def("AllTaus",&ActinMixedNucleates::AllTaus)
        .def("nFreeMonomers", &ActinMixedNucleates::nFreeMonomers)
        .def("NumMonOnEachFiber",&ActinMixedNucleates::NumMonOnEachFiber)
        .def("BranchedOrLinear", &ActinMixedNucleates::BranchedOrLinear)
        .def("BoundBarbedStates",&ActinMixedNucleates::BoundBarbedStates);
}    
