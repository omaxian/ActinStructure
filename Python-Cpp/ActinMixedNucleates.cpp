#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ActinStructure.cpp"

namespace py=pybind11;
typedef py::array_t<int, py::array::c_style | py::array::forcecast> npInt; 
typedef py::array_t<double, py::array::c_style | py::array::forcecast> npDoub; 

class ActinMixedNucleates {
    
    public:
    
    ActinMixedNucleates(uint NMon, intvec NperFiber, vec3 Lengths, double a, double kbT, double mu, vec RxnRates, int seed, int nThr){
    
        /*
        Initialization
        Nmon = total number of monomers, including those attached to fibers
        NperFiber = integer vector of number of monomers on each of the nascent fibers
        Lengths = 3-vector of (x,y,z) periodic domain lengths
        a = hydrodynamic radius
        kbT = thermal energy
        mu = fluid viscosity
        RxnRates = vector with the rates (see below for ordering)
        seed = the seed for reproducibility
        */
        
        // Initialize random number generators
        rngu.seed(seed);
        unifdist = std::uniform_real_distribution<double>(0.0,1.0);
        rng.seed(seed);
        normaldist = std::normal_distribution<double>(0.0,1.0);
        for (int d=0; d<3;d++){
            _Lens[d]=Lengths[d];
        }
                
        // Diffusion variables
        _TotalMonomers = NMon;
        _a = a;
        _kbT = kbT;
        _mu = mu;
        _spacing = _a/2;
        
        // Reaction variables
        // We assume all reaction rates are given in units of s^(-1)
        _DimerOnRate = RxnRates[0];
        _DimerOffRate = RxnRates[1];
        _TrimerOnRate = RxnRates[2];
        _TrimerOffRate = RxnRates[3];
        _TetramerOnRate = RxnRates[4]+RxnRates[6];
        _FiberEndBindingRates = vec(5,1.0);
        _FiberEndBindingRates[0] = RxnRates[4]; // Barbed binding
        _FiberEndBindingRates[1] = RxnRates[5]; // Barbed unbinding
        _FiberEndBindingRates[2] = RxnRates[6]; // Pointed binding
        _FiberEndBindingRates[3] = RxnRates[7]; // Pointed unbinding
        _FiberEndBindingRates[4] = 1;
        _TotalFormins = 0;
        _FreeFormins = 0;
        _nForminNucleates = 0;
        _TotalArp = 0;
        _FreeArp = 0;
        
        // Initialize fibers at end of monomer list (just overwrite what is there already)
        uint nInFibs = std::accumulate(NperFiber.begin(), NperFiber.end(),0);
        _FreeMonomers = _TotalMonomers - nInFibs;
        _nDimers = 0;
        _nTrimers = 0;
        uint nFibs = NperFiber.size();
        for (uint iFib=0; iFib < nFibs; iFib++){
            // Initialize a fiber
            vec3 tau = PointOnUnitSphere();
            vec3 X0 = UniformPointInBox();
            uint ThisnMon = NperFiber[iFib];
            _Fibers.push_back(std::make_shared<Fiber>(X0, tau, ThisnMon, _spacing, _a, _mu, _kbT,_FiberEndBindingRates));
        }
       
    }
    
    void InitializeFormins(uint nFormin, vec ForminRates){
        _TotalFormins = nFormin;
        _FreeFormins = nFormin;
        _ForminNucleationRate = ForminRates[0];
        _ForminBindRate = ForminRates[1];
        _ForminUnbindRate = ForminRates[2];
        _FiberEndBindingRates[4] = ForminRates[3]; // formin enhancement
    }
    
    void InitializeArp(uint nArp, vec ArpRates){
        _TotalArp = nArp;
        _FreeArp = nArp;
        _MinMonForBranching = 4;
        _ArpMonomerFiberBindRate = ArpRates[0];
        _ArpMonomerFiberUnbindRate = ArpRates[1];
    }
        
    
    void Diffuse(double dt){
        // Diffuse the fibers
        for(auto&& FibObj: _Fibers){
            int nRand = FibObj->NumberRand();
            vec RandomNumbers(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=normaldist(rng);
            }
            FibObj->Diffuse(dt,RandomNumbers);
        }
    }
    
    void React(double dt){
        /*
        Event-driven simulation with well-mixed monomers 
        For now only doing FIBERS - no branching yet
        */       
        
        double t = 0;
        while (t < dt){
            uint index = 0;
            double deltaT = TimeNucleationReactions(index);
            deltaT = TimeForminNucleationReactions(index,deltaT);
            int nRxnsPerFiber = 9;
            int numSingleRxns = 6;
            deltaT = TimeFiberBindUnbindReactions(index,deltaT,numSingleRxns,nRxnsPerFiber);
            // All reactions listed -- process the next one
            bool EventHappened = false;
            //std::cout << "Time " << deltaT << " for total " << t << " and index " << index << std::endl;
            if (t+deltaT > dt){
                std::cout << "Free formins and arps (to check) " << _FreeFormins << " , " << _FreeArp << std::endl;
                return;
            }
            if (index < 5){
                // Nucleation reaction
                EventHappened = ProcessNucleationReaction(index);
            } else if (index == 5){
                EventHappened = ProcessForminNucleationReaction();
            } else { // Working with index > 5
                EventHappened = ProcessFiberBindUnbindReaction(index,numSingleRxns,nRxnsPerFiber);
            } // end fiber
            // Check conservation of monomers
            if (_FreeMonomers > _TotalMonomers){
                std::cout << "Error - more monomers than total!" << std::endl;
                return;
            }
            if (EventHappened){
                t+=deltaT;
            } 
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
    
    npInt NumMonOnEachFiber(){
        uint nFib = nTotalFibers();
        intvec Info(3+nFib);
        Info[0]=_FreeMonomers;
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
    
    npInt BranchedOrLinear(bool MothersAreBranched){
        uint nFib = nTotalFibers();
        intvec Info(nFib);
        int iFib=0;
        for(auto&& FibObj: _Fibers){
            int nFibThis = FibObj->nFibers();
            for (int jFib=0; jFib < nFibThis; jFib++){
                Info[iFib]=(nFibThis > 1);
                if (!MothersAreBranched && jFib==0){
                    Info[iFib]=0;
                }    
                iFib++;
            }
        }    
        return makePyArray(Info);
    }
    
    npInt BoundFormins(){
        uint nFib = nTotalFibers();
        intvec Info(nFib);
        int iFib=0;
        for(auto&& FibObj: _Fibers){
            int nFibThis = FibObj->nFibers();
            for (int jFib=0; jFib < nFibThis; jFib++){
                Info[iFib]=FibObj->ForminBound(jFib);
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
        uint _TotalMonomers, _FreeMonomers, _nDimers, _nTrimers;
        uint _TotalFormins, _FreeFormins, _nForminNucleates;
        uint _TotalArp, _FreeArp, _MinMonForBranching;
        double _a, _kbT, _mu, _spacing;
        vec3 _Lens;
        double _DimerOnRate, _DimerOffRate, _TrimerOnRate, _TrimerOffRate,_TetramerOnRate;
        vec _FiberEndBindingRates;
        double _ForminNucleationRate, _ForminBindRate, _ForminUnbindRate, _ForminEnhancement;
        double _ArpMonomerFiberBindRate, _ArpMonomerFiberUnbindRate;
        vec _X; 
        std::normal_distribution<double> normaldist;
        std::mt19937_64 rng;
        std::uniform_real_distribution<double> unifdist;
        std::mt19937_64 rngu;
        
        double TimeNucleationReactions(uint &index){
            // Reaction 0: formation of dimers
            //std::cout << "Time dimer formation " << deltaT << std::endl;
            double RateDimerForm = _DimerOnRate*_FreeMonomers*_FreeMonomers;
            if (_FreeMonomers < 2){
                RateDimerForm = 0;
            }    
            double deltaT = logrand()/RateDimerForm;
            // Reaction 1: break-up of dimers
            double RateDimerBreakup = _DimerOffRate*_nDimers;
            double TryDeltaT = logrand()/RateDimerBreakup; 
            //std::cout << "Time dimer breakup " << TryDeltaT << std::endl;
            if (TryDeltaT < deltaT){
                deltaT = TryDeltaT;
                index = 1;
            }
            // Reaction 2: trimer formation
            double RateTrimerForm = _TrimerOnRate*_FreeMonomers*_nDimers;
            TryDeltaT = logrand()/RateTrimerForm; 
            //std::cout << "Time trimer form " << TryDeltaT << std::endl;
            if (TryDeltaT < deltaT){
                deltaT = TryDeltaT;
                index = 2;
            }
            // Reaction 3: trimer breakup
            double RateTrimerBreakup = _TrimerOffRate*_nTrimers;
            TryDeltaT = logrand()/RateTrimerBreakup; 
            //std::cout << "Time trimer breakup " << TryDeltaT << std::endl;
            if (TryDeltaT < deltaT){
                deltaT = TryDeltaT;
                index = 3;
            }
            // Reaction 4: tetramer formation
            double RateTetramerForm = _TetramerOnRate*_nTrimers*_FreeMonomers;
            TryDeltaT = logrand()/RateTetramerForm; 
            //std::cout << "Time tetramer form " << TryDeltaT << std::endl;
            if (TryDeltaT < deltaT){
                deltaT = TryDeltaT;
                index = 4;
            }
            return deltaT;
        }
        
        bool ProcessNucleationReaction(uint index){
            if (index == 0){
                // Dimer formation
                //std::cout << "Forming dimer " << std::endl;
                _FreeMonomers-=2;
                _nDimers++;
            } else if (index == 1){
                //std::cout << "Removing dimer " << std::endl;
                _FreeMonomers+=2;
                _nDimers--;
            } else if (index == 2){
                //std::cout << "Forming trimer " << std::endl;
                _FreeMonomers--;
                _nDimers--;
                _nTrimers++;
            } else if (index == 3){
                //std::cout << "Removing trimer " << std::endl;
                _FreeMonomers++;
                _nDimers++;
                _nTrimers--;
            } else if (index == 4){
                //std::cout << "Forming tetramer and fiber object " << std::endl;
                // Forming tetramer - establish object for it
                _FreeMonomers--;
                _nTrimers--;
                vec3 tau = PointOnUnitSphere();
                vec3 X0 = UniformPointInBox();
                _Fibers.push_back(std::make_shared<Fiber>(X0, tau, 4, _spacing, _a, _mu, _kbT,_FiberEndBindingRates));
            } else {
                std::cout << "Doing a nucleation reaction but can't find the option?!" << std::endl;
            }
            return true;
        }
        
        double TimeForminNucleationReactions(uint &index, double deltaT){
            // Reaction 5: formin nucleate formation (irreversible reaction; instantiates a fiber object with 2 monomers)
            if (_TotalFormins > 0){
                double RateForminNucleate = _ForminNucleationRate*_FreeMonomers*_FreeMonomers*_FreeFormins;
                double TryDeltaT = logrand()/RateForminNucleate;
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = 5;
                }
            }
            return deltaT;
        }
        
        bool ProcessForminNucleationReaction(){
            _FreeMonomers-=2;
            _FreeFormins--;
            // Initialize a fiber and bind a formin to it
            vec3 tau = PointOnUnitSphere();
            vec3 X0 = UniformPointInBox();
            _Fibers.push_back(std::make_shared<Fiber>(X0, tau, 2, _spacing, _a, _mu, _kbT,_FiberEndBindingRates));
            _Fibers[_Fibers.size()-1]->BindFormin(0);
            return true;
        }
        
        double TimeFiberBindUnbindReactions(uint &index, double deltaT, int numBefore, int nRxnsPerFiber){
            // Add/subtract from existing fibers
            uint nFib = _Fibers.size();
            //std::cout << "Number of fibers " << nFib << std::endl;
            for (uint iFib = 0; iFib < nFib; iFib++){
                // Rxn 0: add to pointed end
                double RateAddPointed =  _Fibers[iFib]->PointedBindingRate()*_FreeMonomers;
                double TryDeltaT = logrand()/RateAddPointed; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = numBefore+nRxnsPerFiber*iFib;
                }
                // Rxn 1: add to barbed end
                double RateAddBarbedNF =  _Fibers[iFib]->BarbedBindingRateNoFormin()*_FreeMonomers;
                TryDeltaT = logrand()/RateAddBarbedNF; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = numBefore+nRxnsPerFiber*iFib+1;
                }
                // Rxn 2: add to formin bound barbed end
                double RateAddBarbedF = _Fibers[iFib]->BarbedBindingRateWithFormin()*_FreeMonomers;
                TryDeltaT = logrand()/RateAddBarbedF; 
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = numBefore+nRxnsPerFiber*iFib+2;
                }
                // Rxn 3: subtract from pointed end
                double RateSubtractPt = _Fibers[iFib]->PointedUnbindingRate();
                TryDeltaT = logrand()/RateSubtractPt; 
                //std::cout << "Time remove from fiber " << TryDeltaT << std::endl;
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = numBefore+nRxnsPerFiber*iFib+3;
                }
                // Rxn 4: subtract from barbed end
                double RateSubtractBrb = _Fibers[iFib]->BarbedUnbindingRate();
                TryDeltaT = logrand()/RateSubtractBrb; 
                //std::cout << "Time remove from fiber " << TryDeltaT << std::endl;
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = numBefore+nRxnsPerFiber*iFib+4;
                }
                // Reactions with formin
                if (_TotalFormins > 0){
                    // Reaction 5: bind formin
                    double ForminBindRate = _Fibers[iFib]->ForminBindRate(_ForminBindRate)*_FreeFormins;
                    TryDeltaT = logrand()/ForminBindRate;
                    if (TryDeltaT < deltaT){
                        deltaT = TryDeltaT;
                        index = numBefore+nRxnsPerFiber*iFib+5;
                    }
                    // Reaction 6: unbind formin
                    double ForminUnbindRate = _Fibers[iFib]->ForminUnbindRate(_ForminUnbindRate);
                    TryDeltaT = logrand()/ForminUnbindRate;
                    if (TryDeltaT < deltaT){
                        deltaT = TryDeltaT;
                        index = numBefore+nRxnsPerFiber*iFib+6;
                    } 
                }
                // Branching reactions
                if (_TotalArp > 0){ 
                    // Reaction 7: bind arp 2/3
                    double BranchRate =  _ArpMonomerFiberBindRate*_FreeArp*_FreeMonomers*
                        _Fibers[iFib]->nFibersEligibleForBranching(_MinMonForBranching);
                    TryDeltaT = logrand()/BranchRate;
                    if (TryDeltaT < deltaT){
                        deltaT = TryDeltaT;
                        index = numBefore+nRxnsPerFiber*iFib+7;
                    } 
                    // Arp23 bind and unbind. There is one reaction. If chosen, randomly select binding and unbinding. 
                    // If binding, randomly select a place to bind and a new tangent vector, and convert from fiber
                    // to branched fiber object if necessary
                    // If unbinding, randomly select one of the arps/branches and remove it. Return it as a fiber objet
                    // which gets copied to the list of fibers 
                }   
            }
            return deltaT;
        }
        
        bool ProcessFiberBindUnbindReaction(uint index, int numBefore, int nRxnsPerFiber){
            // Identify fiber number and if it's addition or subtraction
            int FibNum = (index-numBefore)/nRxnsPerFiber;
            int RxnType = (index-numBefore) % nRxnsPerFiber; // 0 for addition, 1 for subtraction, 2 for formin
            double u = unifdist(rngu);
            if (RxnType==0){ // Add to pointed end
                _FreeMonomers--;
                _Fibers[FibNum]->PointedBindReaction();
            } else if (RxnType < 3){ // Add to barbed end
                _FreeMonomers--;
                _Fibers[FibNum]->BarbedBindReaction(u, RxnType==2); // type 2 means choose one of the formin ends, 1 the free ends
            } else if (RxnType < 5){ // remove from pointed or barbed end
                // =======================================================
                // EXCEPTIONS
                // Handle this in the fiber class (DEAL WITH IT LATER)
                int nMon = _Fibers[FibNum]->NumMonomers(0); // number of monomers on the mother
                int TotalMons = _Fibers[FibNum]->TotalMonomers();
                bool HasFormin = _Fibers[FibNum]->ForminBound(0);
                if (nMon == 4 && !HasFormin){ 
                    //std::cout << "Breaking up a tetramer " << std::endl;
                    _Fibers.erase(_Fibers.begin() + FibNum);   
                    _FreeMonomers+=TotalMons-3;
                    _nTrimers++;
                } else if (nMon == 2 && HasFormin){
                    //std::cout << "Can't break formin bound dimer - event rejected" << std::endl;
                    return false;
                } 
                // ===========================================================
                else {
                    _FreeMonomers++;
                    _Fibers[FibNum]->UnbindReaction(unifdist(rngu),RxnType==4); // 4 is pointed end, 3 is barbed end
                }
            } else if (RxnType==5) {// Formin binding
                _Fibers[FibNum]->BindFormin(unifdist(rngu));
                _FreeFormins--;
            } else if (RxnType==6) {// Formin unbinding
                _Fibers[FibNum]->UnbindFormin(unifdist(rngu));
                _FreeFormins++;
            } else if (RxnType==7) { // Arp 2/3 binding
                // Arp23 binding
                _FreeArp--;
                _FreeMonomers--;
                uint nFibs = _Fibers[FibNum]->nFibers();
                vec3 RandGauss;
                for (int d=0; d<3;d++){
                    RandGauss[d]=normaldist(rng);
                }
                double u1=unifdist(rngu);
                if (nFibs==1){
                    // Convert to a branched fiber
                    _Fibers[FibNum] =  std::make_shared<BranchedFiber>(_Fibers[FibNum].get(),u1,RandGauss);
                } else {
                   // Already a BranchedFiber object. Add a new branch.
                   double u2=unifdist(rngu);
                   _Fibers[FibNum]->addBranch(u1,u2, RandGauss, _MinMonForBranching);
                }
            }
            return true;
        }
               
        vec3 PointOnUnitSphere(){
            vec3 tau;
            for (int d=0; d < 3; d++){
                tau[d] = normaldist(rng);
            }
            normalize(tau);    
            return tau;
        }
        
        vec3 PointInSphere(double r){
            vec3 tau = PointOnUnitSphere();
            double u = unifdist(rngu);
            double cbrtu = std::cbrt(u);
            for (int d=0; d < 3; d++){
                tau[d] = r*cbrtu*tau[d];
            }
            return tau;
        }
        
        vec3 UniformPointInBox(){
            vec3 X0;
            for (int d=0; d < 3; d++){
                X0[d]=unifdist(rngu)*_Lens[d];
            }
            return X0;
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
        .def(py::init<uint, intvec, vec3, double, double, double, vec,int,int>())
        .def("InitializeFormins",&ActinMixedNucleates::InitializeFormins)
        .def("InitializeArp",&ActinMixedNucleates::InitializeArp)
        .def("Diffuse",&ActinMixedNucleates::Diffuse)
        .def("React",&ActinMixedNucleates::React)
        .def("getX", &ActinMixedNucleates::getX)
        .def("NumMonOnEachFiber", &ActinMixedNucleates::NumMonOnEachFiber)
        .def("nTotalFibers", &ActinMixedNucleates::nTotalFibers)
        .def("BranchedOrLinear", &ActinMixedNucleates::BranchedOrLinear)
        .def("BoundFormins",&ActinMixedNucleates::BoundFormins);
}    
