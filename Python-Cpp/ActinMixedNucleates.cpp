#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ActinStructure.cpp"

namespace py=pybind11;
typedef py::array_t<int, py::array::c_style | py::array::forcecast> npInt; 
typedef py::array_t<double, py::array::c_style | py::array::forcecast> npDoub; 

class ActinMixedNucleates {
    
    public:
          
    ActinMixedNucleates(const uint NMon, const vec3 &Lengths, const vec &SpontaneousRxnRates,
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
        _TotalMonomers = NMon;
        _a = a;
        _kbT = kbT;
        _mu = mu;
        _spacing = spacingInTermsOfA*a;
        
        // Reaction variables
        // We assume all reaction rates are given in units of s^(-1)
        _DimerOnRate = SpontaneousRxnRates[0];
        _DimerOffRate = SpontaneousRxnRates[1];
        _TrimerOnRate = SpontaneousRxnRates[2];
        _TrimerOffRate = SpontaneousRxnRates[3];
        _TetramerOnRate = SpontaneousRxnRates[4]+SpontaneousRxnRates[6];
        
        // Fiber binding rates
        _FiberEndBindingRates = vec(4);
        _FiberEndBindingRates[0] = SpontaneousRxnRates[4]; // Barbed binding
        _FiberEndBindingRates[1] = SpontaneousRxnRates[5]; // Barbed unbinding
        _FiberEndBindingRates[2] = SpontaneousRxnRates[6]; // Pointed binding
        _FiberEndBindingRates[3] = SpontaneousRxnRates[7]; // Pointed unbinding
        
        // Initialize variables relating to barbed-end-bound and monomer bound proteins
        _AlphasPointed = vec(1,1.0);
        _AlphasBarbed = vec(1,1.0);
        _nMonomersBoundToEach = intvec(1,NMon);
        _nBarbedBinders = intvec(0);
        _BarbedBindersRates = vec(0);
        _BarbedNucRates = vec(0);
        _FreeBranchers = 0;
        _BranchRates = vec(0);
        
        rngu.seed(seed);
        unifdist = std::uniform_real_distribution<double>(0.0,1.0);
        rng.seed(seed);
        normaldist = std::normal_distribution<double>(0.0,1.0);
    }
    
    void InitializeMonomerBinders(const intvec &nMonBinders, const vec &AlphasPointed){
        _nMonBinders = intvec(nMonBinders.size());
        std::memcpy(_nMonBinders.data(),nMonBinders.data(),nMonBinders.size()*sizeof(int));
        _AlphasPointed = vec(AlphasPointed.size());
        std::memcpy(_AlphasPointed.data(),AlphasPointed.data(),AlphasPointed.size()*sizeof(double)); 
        _nMonomersBoundToEach.resize(nMonBinders.size(),0);
    } 

    void InitializeBarbedBinders(const intvec &nBarbedBinders, const vec &BarbedNucleationRates, 
        const vec &BarbedBindingRates, const vec &AlphasBarbed){
        _nBarbedBinders = intvec(nBarbedBinders.size());
        std::memcpy(_nBarbedBinders.data(),nBarbedBinders.data(),nBarbedBinders.size()*sizeof(int));
        _BarbedNucRates = vec(BarbedNucleationRates.size());
        std::memcpy(_BarbedNucRates.data(),BarbedNucleationRates.data(),BarbedNucleationRates.size()*sizeof(double));
        _BarbedBindersRates = vec(BarbedBindingRates.size());
        std::memcpy(_BarbedBindersRates.data(),BarbedBindingRates.data(),BarbedBindingRates.size()*sizeof(double)); 
        _AlphasBarbed = vec(AlphasBarbed.size());
        if (_AlphasBarbed.size() <= _nBarbedBinders.size()){
            throw std::runtime_error("Need a leading 1 in your alpha array");
        }
        std::memcpy(_AlphasBarbed.data(),AlphasBarbed.data(),AlphasBarbed.size()*sizeof(double)); 
    }
        
   
    void InitializeBranchers(int nBranchers, const vec &BranchRates){
        _FreeBranchers = nBranchers;
        _BranchRates = vec(BranchRates.size());
        std::memcpy(_BranchRates.data(),BranchRates.data(),BranchRates.size()*sizeof(double)); 
    }
    
    //void InitBranchedStructure(intvec nMonsPerBranch, intvec ForminsOn){
        /* 
        Initialize a branched filament (for debugging mostly). The intvector 
        nMonsPerBranch gives the number of monomers on each branch, and the 
        vector ForminsOn tells us if there is a formin on each branch 
        */
        /*
        int nMonsIn = std::accumulate(nMonsPerBranch.begin(), nMonsPerBranch.end(),0);
        _FreeMonomers-=nMonsIn;
        uint nBranches = nMonsPerBranch.size()-1;
        _FreeArp-=nBranches;
        // Initialize fiber
        vec3 tau = PointOnUnitSphere();
        vec3 X0 = UniformPointInBox();
        _Fibers.push_back(std::make_shared<Fiber>(X0, tau, nMonsPerBranch[0], _spacing, _a, _mu, _kbT,_FiberEndBindingRates));
        // Branches
        uint ind = _Fibers.size()-1;
        vec3 RandGauss;
        for (uint jBranch=0; jBranch < nBranches; jBranch++){
            for (int d=0; d<3;d++){
                RandGauss[d]=normaldist(rng);
            }
            double u1=unifdist(rngu);
            double u2=unifdist(rngu);
            if (jBranch==0){
                _Fibers[ind] =  std::make_shared<BranchedFiber>(_Fibers[ind].get(),0.99,RandGauss);
            } else {
                _Fibers[ind]->addBranch(u1,u2,RandGauss);
            }
            for (int iMon=1; iMon < nMonsPerBranch[jBranch+1]; iMon++){
                _Fibers[ind]->BarbedBindReaction(0.99,false); // (always add to the last branch)
            }
        }
        for (uint iB=0; iB <=nBranches; iB++){
            if (ForminsOn[iB] > 0){
                double DeltaBranch = 1.0/(nBranches+1);
                _Fibers[ind]->BindFormin((0.5+iB)*DeltaBranch);
                _FreeFormins--;
            }
        }
    }*/
        
    
    void Diffuse(double dt){
        // Diffuse the fibers
        for(auto&& FibObj: _Fibers){
            FibObj->Diffuse(dt);
        }
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
            double deltaT = TimeNucleationReactions(index);
            int numSingleRxns = 4+_nBarbedBinders.size();
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
                //std::cout << "Min event time " << FiberEventTimes[MinFibIndex] << " vs. nuc time " << deltaT << std::endl;
                if (FiberEventTimes[MinFibIndex] < deltaT){
                    deltaT = FiberEventTimes[MinFibIndex];
                    index = numSingleRxns+1+MinFibIndex;
                }
            }
            if (t+deltaT > dt || deltaT < 0){
                //std::cout << "Free formins and arps (to check) " << _FreeFormins << " , " << _FreeArp << std::endl;
                return;
            }
            if (index <= numSingleRxns){
                ProcessNucleationReaction(index);
            } else {
                int iFib = index - (numSingleRxns+1);
                int StructChange = _Fibers[iFib]->NeedsStructureChange();
                if (StructChange==0){
                    _Fibers[iFib]->ReactNextEvent(_nMonomersBoundToEach, _nMonBinders, _nBarbedBinders,_FreeBranchers);
                } else if (StructChange == -1){
                    // Delete structure or turn branched fiber into linear fiber
                    uint nFibers = _Fibers[iFib]->nFibers();
                    if (nFibers == 1){
                        _Fibers.erase(_Fibers.begin() + iFib);   
                        _nTrimers++;
                        _nMonomersBoundToEach[0]++;
                    } else {
                        // Convert branched filament to linear
                        int nMonomers=0; 
                        vec3 X0, tau;
                        int BarbedProtein = _Fibers[iFib]->PrepareForLinearSwitch(X0,tau,nMonomers);
                        if (nMonomers < 4){
                            std::cout << "Fiber should have more than 4 monomers! " << std::endl;
                        }
                        _Fibers[iFib] =  std::make_shared<Fiber>(_Lens, nMonomers, _FiberEndBindingRates, 
                            _AlphasPointed, _BarbedBindersRates, _AlphasBarbed, _BranchRates,
                            _spacing, _a, _mu, _kbT,-1);
                        _Fibers[iFib]-> SetX0AndTau(X0,tau);
                        _Fibers[iFib]-> SetBoundBarbed(BarbedProtein);
                        _nMonomersBoundToEach[0]++;
                        _FreeBranchers++;
                    }
                } else if (StructChange==1){ // Converting linear filament to branched  
                    _Fibers[iFib] =  std::make_shared<BranchedFiber>(_Fibers[iFib].get());
                    _nMonomersBoundToEach[0]--;
                    _FreeBranchers--;
                } 
            }
            // Check conservation of monomers
            if (_nMonomersBoundToEach[0] > _TotalMonomers){
                std::cout << "Error - more monomers than total!" << std::endl;
                return;
            }
            int TotMon=_nMonomersBoundToEach[0]+2*_nDimers+3*_nTrimers;
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
        Info[0]=_nMonomersBoundToEach[0];
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
        Info[0]=_nMonomersBoundToEach[0];
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
        std::vector <std::shared_ptr<Fiber>> _Fibers;
        int _TotalMonomers, _nDimers, _nTrimers;
        intvec _nMonomersBoundToEach, _nMonBinders, _nBarbedBinders;
        double _DimerOnRate, _DimerOffRate, _TrimerOnRate, _TrimerOffRate,_TetramerOnRate;
        
        vec _FiberEndBindingRates;
        vec _AlphasPointed, _BarbedNucRates, _BarbedBindersRates, _AlphasBarbed;
        
        int _FreeBranchers;
        vec _BranchRates;
        
        double _a, _kbT, _mu, _spacing;
        vec3 _Lens;
        vec _X; 
             
        std::normal_distribution<double> normaldist;
        std::mt19937_64 rng;
        std::uniform_real_distribution<double> unifdist;
        std::mt19937_64 rngu;
        
        double TimeNucleationReactions(int &index){
            double FreeMonomers = 1.0*_nMonomersBoundToEach[0];
            // Reaction 0: formation of dimers
            double RateDimerForm = _DimerOnRate*FreeMonomers*FreeMonomers;
            if (FreeMonomers < 2){
                RateDimerForm = 0;
            }    
            double deltaT = logrand()/RateDimerForm;
            // Reaction 1: break-up of dimers
            double RateDimerBreakup = _DimerOffRate*_nDimers;
            double TryDeltaT = logrand()/RateDimerBreakup; 
            if (TryDeltaT < deltaT){
                deltaT = TryDeltaT;
                index = 1;
            }
            // Reaction 2: trimer formation
            double RateTrimerForm = _TrimerOnRate*FreeMonomers*_nDimers;
            TryDeltaT = logrand()/RateTrimerForm; 
            if (TryDeltaT < deltaT){
                deltaT = TryDeltaT;
                index = 2;
            }
            // Reaction 3: trimer breakup
            double RateTrimerBreakup = _TrimerOffRate*_nTrimers;
            TryDeltaT = logrand()/RateTrimerBreakup; 
            if (TryDeltaT < deltaT){
                deltaT = TryDeltaT;
                index = 3;
            }
            // Reaction 4: tetramer formation
            double RateTetramerForm = _TetramerOnRate*_nTrimers*FreeMonomers;
            TryDeltaT = logrand()/RateTetramerForm; 
            if (TryDeltaT < deltaT){
                deltaT = TryDeltaT;
                index = 4;
            }
            // Reactions 5+: barbed end protein (i.e., formin) mediated nucleation
            for (uint iB = 0; iB < _nBarbedBinders.size(); iB++){
                double RateNucleate = _BarbedNucRates[iB]*FreeMonomers*FreeMonomers*_nBarbedBinders[iB];
                double TryDeltaT = logrand()/RateNucleate;
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = 5+iB;
                }
            }
            return deltaT;
        }
        
        void ProcessNucleationReaction(uint index){
            if (index == 0){
                // Dimer formation
                _nMonomersBoundToEach[0]-=2;
                _nDimers++;
            } else if (index == 1){
                // Dimer break
                _nMonomersBoundToEach[0]+=2;
                _nDimers--;
            } else if (index == 2){
                // Trimer formation
                _nMonomersBoundToEach[0]--;
                _nDimers--;
                _nTrimers++;
            } else if (index == 3){
                // Trimer breakup
                _nMonomersBoundToEach[0]++;
                _nDimers++;
                _nTrimers--;
            } else {
                // Forming tetramer - establish object for it
                if (index == 4){
                    _nMonomersBoundToEach[0]--;
                    _nTrimers--;
                    _Fibers.push_back(std::make_shared<Fiber>(_Lens, 4, _FiberEndBindingRates, 
                        _AlphasPointed, _BarbedBindersRates, _AlphasBarbed, _BranchRates, 
                        _spacing, _a, _mu, _kbT,-1));
                } else {
                    int iB = index-5;
                    _nBarbedBinders[iB]--;
                    _nMonomersBoundToEach[0]-=2;
                    _Fibers.push_back(std::make_shared<Fiber>(_Lens, 2, _FiberEndBindingRates, 
                        _AlphasPointed, _BarbedBindersRates, _AlphasBarbed, _BranchRates, 
                        _spacing, _a, _mu, _kbT,-1));
                    _Fibers[_Fibers.size()-1]->SetBoundBarbed(iB+1);
                }
            }
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
        //.def("InitBranchedStructure",&ActinMixedNucleates::InitBranchedStructure)
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
