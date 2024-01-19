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
        _BarbedMinus = SpontaneousRxnRates[5]; // Barbed unbinding
        _PointedPlus = vec(1,SpontaneousRxnRates[6]); // Pointed binding
        _PointedMinus = SpontaneousRxnRates[7]; // Pointed unbinding
        
        // Initialize variables relating to barbed-end-bound and monomer bound proteins
        _nFreeMon = Nmon;
        _nDimers = intvec(1,0);
        _nTrimers = intvec(1,0);
        _nBarbedProts = 0;
        _nBarbedBinders = intvec(0);
        _BarbedBindersRates = vec(0);
        _FreeBranchers = 0;
        _BranchRates = vec(0);
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
        std::memcpy(_MonBindKEqs.data(),EquilConsts.data(),EquilConsts.size()*sizeof(int));
        _PointedPlus.resize(_nMonProts+1);
        if (_nMonProts != PointedAlphas.size()){
            throw std::runtime_error("Size mismatch in pointed alpha array");
        }
        for (int iProt=1; iProt <= _nMonProts; iProt++){
            _PointedPlus[iProt] = _PointedPlus[0]*PointedAlphas[iProt-1];
        }
    } 

    void InitializeBarbedBinders(const intvec &nBarbedBinders,const vec &BarbedBinderRates, 
        const vec &AlphaDimer, const vec &AlphaTrimer){
        
        int _nBarbedProts = nBarbedBinders.size();
        _nBarbedBinders = intvec(_nBarbedProts);
        std::memcpy(_nBarbedBinders.data(),nBarbedBinders.data(),nBarbedBinders.size()*sizeof(int));
        _BarbedBindersRates = vec(2*_nBarbedProts);
        std::memcpy(_BarbedBindersRates.data(),BarbedBinderRates.data(),BarbedBinderRates.size()*sizeof(int));
        
        _DimerMinus.resize(_nBarbedProts+1);
        _TrimerMinus.resize(_nBarbedProts+1);
        if (AlphaDimer.size() != _nBarbedProts){
            throw std::runtime_error("Size mismatch in barbed alpha array");
        }
        for (int iProt=1; iProt <= _nBarbedProts; iProt++){
            _DimerMinus[iProt] = _DimerMinus[0]*AlphaDimer[iProt-1];
            _TrimerMinus[iProt]= _TrimerMinus[0]*AlphaTrimer[iProt-1];
        }

        _nDimers.resize(_nBarbedProts+1); = intvec(_nBarbedProts+1,0);
        _nTrimers.resize(_nBarbedProts+1);
         // The rates with proteins attached to monomers or at barbed ends
    }
    
    void InitializeRateMatrices(const vec &DimerAlphas, const vec &TrimerAlphas, const vec &BarbedEndAlphas){
        // Matrices of rates (i=monomer bound protein, j = barbed bound protein)
        int nBarbedDim = _nBarbedBinders.size()+1;
        int nMonDim = _nMonBinders.size()+1;
        int ExpectedDim = nBarbedDim*nMonDim;
        if (DimerAlphas.size() != ExpectedDim || TrimerAlphas.size() != ExpectedDim || BarbedEndAlphas.size() != ExpectedDim){
            throw std::runtime_error("Size mismatch in matrix of alphas");
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
    }
        
   
    void InitializeBranchers(int nBranchers, const vec &BranchRates){
        _FreeBranchers = nBranchers;
        _BranchRates = vec(BranchRates.size());
        std::memcpy(_BranchRates.data(),BranchRates.data(),BranchRates.size()*sizeof(double)); 
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
                    _BarbedBindersRates, _PointedPlus, _BarbedPlus, _PointedMinus, _BarbedMinus, _BranchRates,
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
        intvec NToEach(_nMonBinders.size()+1,nFreeMon)
        if (_nMonBinders.size()==0){
            return NToEach;
        }
        // Now we are in the case when there are multiple proteins
        NMonProt = _nMonBinders.size();
        double denom=1;
        for (int j=0; j < NMonProt; j++){
            denom+=_MonBindKEqs[j]*_nMonBinders[j];
        }
        int nAlreadyBd=0;
        for (int j=NMonProt; j > 0; j--){
            NToEach[j]=1.0/denom*nFreeMon*_nMonBinders[j-1];
            nAlreadyBd+=NToEach[j];
        }
        NToEach[0] = nFreeMon-nAlreadyBd;
        std::cout << "Number free monomers (no protein) = " << NToEach[0] << std::endl;
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
                    FiberEventTimes[iFib] =  _Fibers[iFib]->NextEventTime(_nMonomersBoundToEach, _nBarbedBinders, _FreeBranchers);
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
            //std::cout << "Event time " << deltaT << " index " << index << std::endl;
            if (t+deltaT > dt || deltaT < 0){
                //std::cout << "Free formins and arps (to check) " << _FreeFormins << " , " << _FreeArp << std::endl;
                return;
            }
            if (index < numNucRxns){
                ProcessNucleateReaction(index);
            } else if (index < numNumRxns+numBindRxns){
                ProcessNucleateBarbedBindReaction(index);
            } else {
                int iFib = index - (numSingleRxns+1);
                int StructChange = _Fibers[iFib]->NeedsStructureChange();
                if (StructChange==0){
                    _Fibers[iFib]->ReactNextEvent(_nMonomersBoundToEach, _nBarbedBinders,_FreeBranchers);
                } else if (StructChange == -1){
                    // Delete structure or turn branched fiber into linear fiber
                    uint nFibers = _Fibers[iFib]->nFibers();
                    if (nFibers == 1){
                        _Fibers.erase(_Fibers.begin() + iFib);   
                        _nTrimers++;
                        _nFreeMon++;
                    } else {
                        // Convert branched filament to linear
                        int nMonomers=0; 
                        vec3 X0, tau;
                        int BarbedProtein = _Fibers[iFib]->PrepareForLinearSwitch(X0,tau,nMonomers);
                        if (nMonomers < 4){
                            std::cout << "Fiber should have more than 4 monomers! " << std::endl;
                        }
                        _Fibers[iFib] =  std::make_shared<Fiber>(_Lens, nMonomers, _nMonProts, _nBarbedProts, 
                            _BarbedBindersRates, _PointedPlus, _BarbedPlus, _PointedMinus, _BarbedMinus, _BranchRates,
                            _spacing, _a, _mu, _kbT,-1);
                        _Fibers[iFib]-> SetX0AndTau(X0,tau);
                        _Fibers[iFib]-> SetBoundBarbed(0,BarbedProtein);
                        _nFreeMon++;
                        _FreeBranchers++;
                    }
                } else if (StructChange==1){ // Converting linear filament to branched  
                    _Fibers[iFib] =  std::make_shared<BranchedFiber>(_Fibers[iFib].get());
                    _nFreeMon--;
                    _FreeBranchers--;
                } 
            }
            // Check conservation of monomers
            int TotMon=_nFreeMon+2*_nDimers+3*_nTrimers;
            for(auto&& FibObj: _Fibers){
                TotMon+=FibObj->TotalMonomers();
            }
            //std::cout << "Free arps " << _FreeBranchers << std::endl;
            //std::cout << "Free formins " << _nBarbedBinders[0] << std::endl;
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
    
    npInt NumMonOnEachFiber(){
        uint nFib = nTotalFibers();
        intvec Info(3+nFib);
        Info[0]=_nFreeMon
        Info[1]=_nDimers;
        Info[2]=_nTrimers;
        int iFib=0;
        for(auto&& FibObj: _Fibers){
            int nFibThis = FibObj->nFibers();
            for (int jFib=0; jFib < nFibThis; jFib++){
                Info[3+iFib]=FibObj->NumMonomers(jFib);
                iFib++;
            }
        }    
        return makePyArray(Info);
    }
    
    npInt NumMonOnEachStructure(){
        uint nStruct = _Fibers.size();
        intvec Info(3+nStruct);
        Info[0]=_nFreeMon;
        Info[1]=_nDimers;
        Info[2]=_nTrimers;
        int iFib=0;
        for(auto&& FibObj: _Fibers){
            Info[3+iFib]=FibObj->TotalMonomers();
            iFib++;
        }    
        return makePyArray(Info);
    }
    
    npInt BranchedOrLinear(bool MothersAreBranched){
        uint nFib = nTotalFibers();
        intvec Info(nFib);
        int iFib=0;
        for(auto&& FibObj: _Fibers){
            int nFibThis = FibObj->nFibers();
            for (int jFib=0; jFib < nFibThis; jFib++){
                Info[iFib]=(nFibThis > 1);
                if (jFib==0 && nFibThis > 1){
                    Info[iFib]=2; // 2 for mothers, 1 for branches
                }    
                iFib++;
            }
        }    
        return makePyArray(Info);
    }
    
    npInt BranchedOrLinearStruct(){
        uint nStruct = _Fibers.size();
        intvec Info(nStruct);
        int iFib=0;
        for(auto&& FibObj: _Fibers){
            int nFibThis = FibObj->nFibers();
            Info[iFib]=(nFibThis > 1);
            iFib++;
        }    
        return makePyArray(Info);
    }
    
    npInt BoundBarbedStates(){
        uint nFib = nTotalFibers();
        intvec Info(nFib);
        int iFib=0;
        for(auto&& FibObj: _Fibers){
            int nFibThis = FibObj->nFibers();
            for (int jFib=0; jFib < nFibThis; jFib++){
                Info[iFib]=FibObj->getBoundBarbed(jFib);
                iFib++;
            }
        }     
        return makePyArray(Info); 
    }
    
    uint nTotalFibers(){
        uint nFib = 0;
        for(auto&& FibObj: _Fibers){
            nFib+=FibObj->nFibers();
        }
        return nFib;     
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
        vec _DimerPlus _TrimerPlus, _BarbedPlus; // matrices of rates (i=monomer bound protein, j = barbed bound protein)
        vec _DimerMinus, _TrimerMinus; // vector of rates only depending on barbed protein
        vec _PointedPlus; // vector of rates only depending on monomer protein
        double _BarbedMinus, _PointedMinus;
        
        // Branchers
        int _FreeBranchers;
        vec _BranchRates;
               
        double TimeNucleationReactions(int &index, const intvec &NumBoundToEach){
            int NumPossBarbed = _nBarbedProts+1;
            int NumPossMon = _nMonProts+1;
            double deltaT = 1.0/0.0;
            for (int iBp = 0; iBp < NumPossBarbed; iBp++){
                // Reaction 0: formation of dimers
                for (int iMp = 0; iMp < NumPossMon; iMp++){
                    int matindex = MatrixIndex(iMp, iBp);
                    double RateDimerForm = _DimerPlus[matindex]*NumBoundToEach[iMp]*NumBoundToEach[iMp];
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
                double RateDimerBreakup = _DimerOffRate[iBp]*_nDimers[iBp];
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
                double RateTrimerBreakup = _TrimerOffRate[iBp]*_nTrimers[iBp];
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
                TryDeltaT = logrand()/RateBindToDimer; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = FirstIndex+iBp;
                }
                double RateUnbindFromDimer = _BarbedBindersRates[2*iBp+1]*_nDimers[iBp];
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
                double RateUnbindTrimer = _BarbedBindersRates[2*iBp+1]*_nTrimers[iBp];
                TryDeltaT = logrand()/RateUnbindTrimer; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index =  FirstIndex+3*_nBarbedProts+iBp;
                }             
            }
            return deltaT;
        }
        
        void ProcessNucleateReaction(const uint &index){
            int NumPossBarbed = _nBarbedProts+1;
            int ProtOnBarbed = index % NumPossBarbed;
            std::cout << "The index " << index << " so I'm putting " << ProtOnBarbed << " on the barbed end " << std::endl;
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
            } else if (index == 4*NumPossBarbed){
                // Forming tetramer - establish object for it
                    _nFreeMon--;
                    _nTrimers[ProtOnBarbed]--;
                    _Fibers.push_back(std::make_shared<Fiber>(_Lens, nMonomers, _nMonProts, _nBarbedProts, 
                            _BarbedBindersRates, _PointedPlus, _BarbedPlus, _PointedMinus, _BarbedMinus, _BranchRates,
                            _spacing, _a, _mu, _kbT,-1));
                    int iB = index-4;
                    _Fibers[_Fibers.size()-1]->SetBoundBarbed(0,iB);
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
                _nDimers[iBp]++;
                _nBarbedBinders[iBp]--;
            } else if (RelIndex < 2*_nBarbedProts){
                // Unbinding from dimer
                _Dimers[0]++;
                _nDimers[iBp]--;
                _nBarbedBinders[iBp]++;
            } else if (RelIndex < 3*_nBarbedProts){
                // Binding to trimer
                _nTrimers[0]--;
                _nTrimers[iBp]++;
                _nBarbedBinders[iBp]--;
            } else { 
                // Unbinding from trimer
                _nTrimers[0]++;
                _nTrimers[iBp]--;
                _nBarbedBinders[iBp]++;
            }
        }
                        
        int MatrixIndex(int iMB, int iBB){
            return (_nBarbedBinders.size()+1)*iMB+iBB;     
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
        .def("InitializeBranchers",&ActinMixedNucleates::InitializeBranchers)
        .def("InitBranchedStructure",&ActinMixedNucleates::InitBranchedStructure)
        .def("Diffuse",&ActinMixedNucleates::Diffuse)
        .def("React",&ActinMixedNucleates::React)
        .def("getX", &ActinMixedNucleates::getX)
        .def("AllX0",&ActinMixedNucleates::AllX0)
        .def("AllTaus",&ActinMixedNucleates::AllTaus)
        .def("NumMonOnEachFiber",&ActinMixedNucleates::NumMonOnEachFiber)
        .def("NumMonOnEachStructure",&ActinMixedNucleates::NumMonOnEachStructure)
        .def("nTotalFibers", &ActinMixedNucleates::nTotalFibers)
        .def("BranchedOrLinear", &ActinMixedNucleates::BranchedOrLinear)
        .def("BranchedOrLinearStruct",&ActinMixedNucleates::BranchedOrLinearStruct)
        .def("BoundBarbedStates",&ActinMixedNucleates::BoundBarbedStates);
}    
