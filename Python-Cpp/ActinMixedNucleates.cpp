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
        _spacing = 2*_a;
        
       
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
            _Fibers.push_back(Fiber(X0, tau, ThisnMon, _spacing, _a, _mu, _kbT));
        }

        // Reaction variables
        // We assume all reaction rates are given in units of s^(-1)
        _DimerOnRate = RxnRates[0];
        _DimerOffRate = RxnRates[1];
        _TrimerOnRate = RxnRates[2];
        _TrimerOffRate = RxnRates[3];
        _BarbedBindingRate = RxnRates[4];
        _BarbedUnbindingRate = RxnRates[5];
        _PointedBindingRate = RxnRates[6];
        _PointedUnbindingRate = RxnRates[7];
        _TotalFormins = 0;
        _FreeFormins = 0;
        _nForminNucleates = 0;
        
    }
    
    void InitializeFormins(uint nFormin, vec ForminRates){
        _TotalFormins = nFormin;
        _FreeFormins = nFormin;
        _ForminNucleationRate = ForminRates[0];
        _ForminBindRate = ForminRates[1];
        _ForminUnbindRate = ForminRates[2];
        _ForminEnhancement = ForminRates[3];
    }
        
    
    void Diffuse(double dt){
        // Diffuse the fibers
        for (uint i=0; i < _Fibers.size(); i++){
            int nRand = _Fibers[i].NumberRand();
            vec RandomNumbers(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=normaldist(rng);
            }
            _Fibers[i].Diffuse(dt,RandomNumbers);
        }
        // Diffuse the branched fibers
        for (uint i=0; i < _BranchedFibers.size(); i++){
            int nRand = _BranchedFibers[i].NumberRand();
            vec RandomNumbers(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=normaldist(rng);
            }
            _BranchedFibers[i].Diffuse(dt,RandomNumbers);
        }
    }
    
    int React(double dt){
        /*
        Event-driven simulation with well-mixed monomers 
        For now only doing FIBERS - no branching yet
        */       
        
        double t = 0;
        while (t < dt){
            uint index = 0;
            double deltaT = TimeNucleationReactions(index);
            deltaT = TimeForminNucleationReactions(index,deltaT);
            int nRxnsPerFiber = 3;
            int numSingleRxns = 6;
            deltaT = TimeFiberBindUnbindReactions(index,deltaT,numSingleRxns,nRxnsPerFiber);
            // All reactions listed -- process the next one
            bool EventHappened = false;
            //std::cout << "Time " << deltaT << " for total " << t << " and index " << index << std::endl;
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
                return 0;
            }
            if (EventHappened){
                t+=deltaT;
            } 
        } // end while loop   
        return _Fibers.size();           
    }
    
    npDoub getX(){
        uint nFib = _Fibers.size();
        vec AllX(0);
        for (uint iFib=0; iFib < nFib; iFib++){
            vec XFib = _Fibers[iFib].getX();
            AllX.insert(AllX.end(), XFib.begin(), XFib.end());
        }
        return makePyDoubleArray(AllX);
    }
    
    npInt getStructureInfo(){
        uint nFib = _Fibers.size();
        intvec Info(3+nFib);
        Info[0]=_FreeMonomers;
        Info[1]=_nDimers;
        Info[2]=_nTrimers;
        for (uint iFib=0; iFib < nFib; iFib++){
            Info[3+iFib]=_Fibers[iFib].NumMonomers();
        }    
        return makePyArray(Info);
    }
    
    npInt getBoundFormins(){
        uint nFib = _Fibers.size();
        intvec Info(nFib);
        for (uint iFib=0; iFib < nFib; iFib++){
            Info[iFib]=_Fibers[iFib].ForminBound();
        }    
        return makePyArray(Info); 
    }
    
        
    private:
        std::vector <Fiber> _Fibers;
        std::vector <BranchedFiber> _BranchedFibers;
        uint _TotalMonomers, _FreeMonomers, _nDimers, _nTrimers;
        uint _TotalFormins, _FreeFormins, _nForminNucleates;
        double _a, _kbT, _mu, _spacing;
        vec3 _Lens;
        double _DimerOnRate, _DimerOffRate, _TrimerOnRate, _TrimerOffRate;
        double _BarbedBindingRate, _BarbedUnbindingRate, _PointedBindingRate, _PointedUnbindingRate;
        double _ForminNucleationRate, _ForminBindRate, _ForminUnbindRate, _ForminEnhancement;
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
            double RateTetramerForm = (_BarbedBindingRate + _PointedBindingRate)*_nTrimers*_FreeMonomers;
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
                _Fibers.push_back(Fiber(X0, tau, 4, _spacing, _a, _mu, _kbT));
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
            _Fibers.push_back(Fiber(X0, tau, 2, _spacing, _a, _mu, _kbT));
            _Fibers[_Fibers.size()-1].BindFormin();
            return true;
        }
        
        double TimeFiberBindUnbindReactions(uint &index, double deltaT, int numBefore, int nRxnsPerFiber){
            // Add/subtract from existing fibers
            uint nFib = _Fibers.size();
            //std::cout << "Number of fibers " << nFib << std::endl;
            for (uint iFib = 0; iFib < nFib; iFib++){
                bool BoundFormin = _Fibers[iFib].ForminBound();
                double BarbedBindRate = _BarbedBindingRate;
                if (BoundFormin){
                    BarbedBindRate*=_ForminEnhancement;
                }    
                double RateAdd = (BarbedBindRate + _PointedBindingRate)*_FreeMonomers;
                int nMon = _Fibers[iFib].NumMonomers();
                /*if (nMon == 5){
                    std::cout << "MAX 5 MONOMERS " << std::endl;
                    RateAdd = 0;
                }*/
                double TryDeltaT = logrand()/RateAdd; 
                //std::cout << "Time add to fiber " << TryDeltaT << std::endl;
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = numBefore+nRxnsPerFiber*iFib;
                }
                double RateSubtract = _BarbedUnbindingRate + _PointedUnbindingRate;
                TryDeltaT = logrand()/RateSubtract; 
                //std::cout << "Time remove from fiber " << TryDeltaT << std::endl;
                if (TryDeltaT < deltaT){
                    deltaT = TryDeltaT;
                    index = numBefore+nRxnsPerFiber*iFib+1;
                }
                // Reaction with formin
                if (_TotalFormins > 0){
                    // There is one rate. Unbinding if already bound, binding otherwise.
                    double ForminRate = _ForminUnbindRate;
                    if (!BoundFormin){
                        ForminRate = _ForminBindRate*_FreeFormins;
                    } else if (nMon < 4){
                        ForminRate = 0;
                    }
                    TryDeltaT = logrand()/ForminRate;
                    if (TryDeltaT < deltaT){
                        deltaT = TryDeltaT;
                        index = numBefore+nRxnsPerFiber*iFib+2;
                    } 
                }
            }
            return deltaT;
        }
        
        bool ProcessFiberBindUnbindReaction(uint index, int numBefore, int nRxnsPerFiber){
            // Identify fiber number and if it's addition or subtraction
            int FibNum = (index-numBefore)/nRxnsPerFiber;
            int RxnType = (index-numBefore) % nRxnsPerFiber; // 0 for addition, 1 for subtraction, 2 for formin
            int nMon = _Fibers[FibNum].NumMonomers();
            bool HasFormin = _Fibers[FibNum].ForminBound();
            //std::cout << "Fiber number " << FibNum << std::endl;
            if (RxnType==2){
                if (HasFormin){
                    _FreeFormins++;
                    _Fibers[FibNum].UnbindFormin();  
                } else {
                    _FreeFormins--;
                    _Fibers[FibNum].BindFormin();
                }  
            } else if (RxnType==1){
                //std::cout << "Subtraction!" << std::endl;
                // If it's a tetramer that is breaking up, then remove it from the list and be done
                if (nMon == 4 && !HasFormin){ 
                    //std::cout << "Breaking up a tetramer " << std::endl;
                    _Fibers.erase(_Fibers.begin() + FibNum);   
                    _FreeMonomers++;
                    _nTrimers++;
                } else if (nMon == 2 && HasFormin){
                    //std::cout << "Can't break formin bound dimer - event rejected" << std::endl;
                    return false;
                } else { // Structure stays the same
                    _FreeMonomers++;
                    //std::cout << "Removing monomer from longer fiber " << std::endl;
                    double pPointed = _PointedUnbindingRate/(_PointedUnbindingRate+_BarbedUnbindingRate);
                    //std::cout << "Prob pointed end " << pPointed << std::endl;
                    bool FromPointedEnd = unifdist(rngu) < pPointed;
                    //std::cout << "FromPointedEnd " << FromPointedEnd << std::endl;
                    _Fibers[FibNum].removeMonomer(FromPointedEnd);
                }
            } else if (RxnType==0){ // Addition
                _FreeMonomers--;
                //std::cout << "Addition!" << std::endl;
                bool BoundFormin = _Fibers[FibNum].ForminBound();
                double BarbedBindRate = _BarbedBindingRate;
                if (BoundFormin){
                    BarbedBindRate*=_ForminEnhancement;
                }    
                double pPointed = _PointedBindingRate/(_PointedBindingRate+BarbedBindRate);
                //std::cout << "Prob pointed end " << pPointed << std::endl;
                bool ToPointedEnd = unifdist(rngu) < pPointed;
                //std::cout << "ToPointedEnd " << ToPointedEnd  << std::endl;
                _Fibers[FibNum].addMonomer(ToPointedEnd);
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
        .def("Diffuse",&ActinMixedNucleates::Diffuse)
        .def("React",&ActinMixedNucleates::React)
        .def("getX", &ActinMixedNucleates::getX)
        .def("getStructureInfo", &ActinMixedNucleates::getStructureInfo)
        .def("getBoundFormins",&ActinMixedNucleates::getBoundFormins);
}    
