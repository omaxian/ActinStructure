#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ActinStructure.cpp"
#include "NeighborSearcher.cpp"

namespace py=pybind11;
typedef py::array_t<int, py::array::c_style | py::array::forcecast> npInt; 
typedef py::array_t<double, py::array::c_style | py::array::forcecast> npDoub; 

class AllTheActin{
    
    public:
    
    AllTheActin(uint NMon, intvec NperFiber, vec3 Lengths, double a, double kbT, double mu, vec RxnRates, int seed, int nThr)
    :_ReactionNeighbors(){
    
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
        
                
        // Diffusion variables
        _TotalMonomers = NMon;
        _a = a;
        _kbT = kbT;
        _mu = mu;
        _X = vec(3*_TotalMonomers);
        _spacing = 2*_a;
        
        // Initialize monomers
        for (uint i=0; i < _TotalMonomers; i++){
            vec3 MonPt;
            for (int d=0; d<3; d++){
               MonPt[d]=unifdist(rngu)*Lengths[d];
            }
            _Monomers.push_back(Monomer(MonPt,a,mu,kbT));
        }
        _StructureIndex = intvec(_TotalMonomers,-1);
        _BarbedEnd = std::vector<bool>(_TotalMonomers,false);
        _PointedEnd = std::vector<bool>(_TotalMonomers,false);
        
        // Initialize fibers at end of monomer list (just overwrite what is there already)
        uint nInFibs = std::accumulate(NperFiber.begin(), NperFiber.end(),0);
        uint nFibs = NperFiber.size();
        uint MonStart = NMon-nInFibs;
        for (uint iFib=0; iFib < nFibs; iFib++){
            // Initialize a fiber
            vec3 tau = PointOnUnitSphere();
            vec MonP = _Monomers[MonStart].getX();
            vec3 X0; 
            for (int d=0; d < 3; d++){
                X0[d]=MonP[d];
            }
            uint ThisnMon = NperFiber[iFib];
            intvec MonIndices(ThisnMon);
            for (uint iM=0; iM < ThisnMon; iM++){
                MonIndices[iM]=MonStart+iM;
                _StructureIndex[MonStart+iM]=iFib;
            }
            _PointedEnd[MonStart]=true;
            _BarbedEnd[MonStart+ThisnMon-1]=true;
            _Fibers.push_back(Fiber(X0, tau, ThisnMon, MonIndices, _spacing, _a, _mu, _kbT));
            MonStart+=ThisnMon;
        }
        UpdateMonomerLocations();

        // Reaction variables
        _BindingRadius = RxnRates[0];
        _TwoMonRate = RxnRates[1];
        _TwoMonOffRate = RxnRates[2];
        _BarbedBindingRate = RxnRates[3];
        _BarbedUnbindingRate = RxnRates[4];
        _PointedBindingRate = RxnRates[5];
        _PointedUnbindingRate = RxnRates[6];
        _ReactionNeighbors = NeighborSearcher(_TotalMonomers, Lengths, _BindingRadius, nThr);
        
    }
    
    void Diffuse(double dt){
        int start=0; 
        // Diffuse monomers
        for (uint i=0; i < _TotalMonomers; i++){
            int nRand = _Monomers[i].NumberRand();
            vec RandomNumbers(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=normaldist(rng);
            }
            start+=nRand;    
            if (_StructureIndex[i]==-1){
                _Monomers[i].Diffuse(dt,RandomNumbers);
            }
        }
        // Diffuse the fibers
        for (uint i=0; i < _Fibers.size(); i++){
            int nRand = _Fibers[i].NumberRand();
            vec RandomNumbers(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=normaldist(rng);
            }
            start+=nRand;
            _Fibers[i].Diffuse(dt,RandomNumbers);
        }
        // Diffuse the branched fibers
        for (uint i=0; i < _BranchedFibers.size(); i++){
            int nRand = _BranchedFibers[i].NumberRand();
            vec RandomNumbers(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=normaldist(rng);
            }
            start+=nRand;    
            _BranchedFibers[i].Diffuse(dt,RandomNumbers);
        }
        // Keep monomers in sync
        UpdateMonomerLocations();
    }
    
    int React(double dt){
        /*
        FOR NOW WE ARE ONLY CONSIDERING LINEAR FILAMENTS AND MONOMERS
        This uses first order kinetics -- it processes binding reactions
        with rate r, then the reaction occurs with probability r*dt during 
        the time step. Reactions can be blocked 
        */
        AssembleX();
        intvec NeighborList = _ReactionNeighbors.NeighborList(_X);
        uint nPair = NeighborList.size()/2;
         
        int newStructIndex =  _Fibers.size()+_BranchedFibers.size();
        
        // Binding reactions
        for (uint iPair = 0; iPair < nPair; iPair++){
            int pt1 = NeighborList[2*iPair];
            int pt2 = NeighborList[2*iPair+1];
            int struct1 = _StructureIndex[pt1];
            int struct2 = _StructureIndex[pt2];
            if (struct1 >=0 && struct2 >=0){
                // Do nothing (two fibers can't come together)
                ////std::cout << "Two fibers; doing nothing" << std::endl;
            } else if (struct1 == -1 && struct2 == -1){
                // 2 monomers forming a fiber
                bool Reacted = FormFiberFromMonomers(pt1, pt2, dt,newStructIndex);
                if (Reacted){
                    //std::cout << "Add mon+mon, New structure index " << newStructIndex << std::endl;
                }
            } else { // There is one monomer and a fiber
                bool Reacted = AddMonomerToFiber(struct1, pt1, struct2, pt2, dt);
                if (Reacted){
                    //std::cout << "Add mon+fib, New structure index " << newStructIndex << std::endl;
                }
            } 
        }
        
        uint nFib = _Fibers.size();
        // Process unbinding reactions
        for (int iFib = nFib-1; iFib >= 0 ; iFib--){
            NMonomerFibBreakup(iFib, dt);
        }
        
        return _Fibers.size();
    }
    
    npDoub getX(){
        AssembleX();
        return makePyDoubleArray(_X);
    }
    
    npInt getStructureIDs(){
        return makePyArray(_StructureIndex);
    }
           
        
    private:
        std::vector <Monomer> _Monomers;
        std::vector <Fiber> _Fibers;
        std::vector <BranchedFiber> _BranchedFibers;
        uint _TotalMonomers;
        double _a, _kbT, _mu, _spacing;
        double _BindingRadius, _TwoMonRate, _TwoMonOffRate, _BarbedBindingRate, _BarbedUnbindingRate, _PointedBindingRate, _PointedUnbindingRate;
        intvec _StructureIndex;
        std::vector<bool>  _BarbedEnd, _PointedEnd;
        vec _X; 
        std::normal_distribution<double> normaldist;
        std::mt19937_64 rng;
        std::uniform_real_distribution<double> unifdist;
        std::mt19937_64 rngu;
        NeighborSearcher _ReactionNeighbors;
        
        
        void AssembleX(){
            for (uint i=0; i < _TotalMonomers; i++){
                vec XStruct = _Monomers[i].getX();
                for (uint j=0; j < 3; j++){
                    _X[3*i+j]=XStruct[j];
                }
            }
        } 
        
        void UpdateMonomerLocations(){
            // Copy from fibers
            for (uint i=0; i < _Fibers.size(); i++){
                vec FibX = _Fibers[i].getX();
                intvec fibMons = _Fibers[i].getMonomerIndices();
                for (uint i =0; i < fibMons.size(); i++){
                    uint MonIndex =  fibMons[i];
                    vec3 XForMon;
                    for (int d=0; d< 3; d++){
                        XForMon[d] = FibX[3*i+d];
                    }
                    _Monomers[MonIndex].setX0(XForMon);
                }
            }
            // Copy from branched fibers
            for (uint i=0; i < _BranchedFibers.size(); i++){
                vec FibX = _BranchedFibers[i].getX();
                intvec fibMons = _BranchedFibers[i].getMonomerIndices();
                for (uint i =0; i < fibMons.size(); i++){
                    uint MonIndex =  fibMons[i];
                    vec3 XForMon;
                    for (int d=0; d< 3; d++){
                        XForMon[d] = FibX[3*i+d];
                    }
                    _Monomers[MonIndex].setX0(XForMon);
                }
            }
        }
            
        
        bool FormFiberFromMonomers(int pt1, int pt2, double dt, int &newStructIndex){
            double pForm = dt*_TwoMonRate;
            //std::cout << "Two monomers; forming fiber number " << newStructIndex << " with probability " << pForm << std::endl;
            double r = unifdist(rngu);
            if (r < pForm) { 
                // Randomly switch the pointed end
                double re = unifdist(rngu);
                if (re < 0.5){
                   int temp = pt1;
                   pt1 = pt2; 
                   pt2 = temp;
                }
                // The pointed end is pt1
                ////std::cout << "Pointed end at pt " << pt1 << std::endl;
                vec X0 = _Monomers[pt1].getX();
                vec3 tau, X0Dim;
                for (int d = 0; d < 3; d++){
                    tau[d] = normaldist(rng);
                    X0Dim[d] = X0[d];
                }
                normalize(tau);
                intvec MonomerInds = {pt1,pt2};
                _Fibers.push_back(Fiber(X0Dim,tau,2,MonomerInds,_spacing,_a,_mu,_kbT));
                // Update monomer locations for barbed end
                _Monomers[pt2].setX0(_Fibers[newStructIndex].getBarbedEnd());
                _StructureIndex[pt1] = newStructIndex;
                _StructureIndex[pt2] = newStructIndex;
                _BarbedEnd[pt2] = true;
                ////std::cout << "Barbed end at point " << pt2 << std::endl;
                _PointedEnd[pt1] = true;
                ////std::cout << "Pointed end at point " << pt1 << std::endl;
                newStructIndex++;
                return true;
            }
            return false;
        }
        
        bool AddMonomerToFiber(const int struct1, const int pt1, const int struct2, const int pt2, double dt){
            // Find the monomer
            int MonIndex = pt2;
            int FibIndex = pt1;
            if (struct1 == -1){
                MonIndex = pt1;
                FibIndex = pt2;
            }
            ////std::cout << "The fiber is at point " << FibIndex << std::endl;
            double RxnRate = 0;
            if (_BarbedEnd[FibIndex]){
                RxnRate = _BarbedBindingRate;  
                ////std::cout << "Barbed end at point " << FibIndex << std::endl;
            } else if (_PointedEnd[FibIndex]){
                RxnRate = _PointedBindingRate;
                ////std::cout << "Pointed end at point " << FibIndex << std::endl;
            }  
            // Only do it if the number of monomers > 4 (TEMPORARY)
            double pForm = dt*RxnRate;
            int cExistingFibIndex = _StructureIndex[FibIndex]; // The structure to join
            int nMon = _Fibers[cExistingFibIndex].NumMonomers();
            if (nMon ==5){
                pForm=0;
                //std::cout << "5 monomer max!" << std::endl;
            }
            ////std::cout << "Adding monomer to fiber with probability " << pForm << std::endl;
            double r = unifdist(rngu);
            if (r < pForm) { 
                // Process the reaction
                ////std::cout << "Adding to existing fiber " << std::endl;
                int ExistingFibIndex = _StructureIndex[FibIndex]; // The structure to join
                ////std::cout << "Point " << FibIndex << " adding monomer # " << MonIndex << " to fiber " << ExistingFibIndex << std::endl;
                _Fibers[ExistingFibIndex].addMonomer(_PointedEnd[FibIndex], MonIndex);
                _StructureIndex[MonIndex] = ExistingFibIndex;
                if (_BarbedEnd[FibIndex]){
                    ////std::cout << "New barbed end" << std::endl;
                    _BarbedEnd[FibIndex]=false;
                    _BarbedEnd[MonIndex]=true;
                } else if (_PointedEnd[FibIndex]){
                    ////std::cout << "New pointed end" << std::endl;
                    _PointedEnd[FibIndex]=false;
                    _PointedEnd[MonIndex]=true;
                }
                return true;
            }
            return false;
        }
        
        void TwoMonomerFibBreakup(const int iFib){
            ////std::cout << "Breaking up fiber " << std::endl;
            vec Xfib = _Fibers[iFib].getX();
            intvec MonIndices = _Fibers[iFib].getMonomerIndices();
            ////std::cout << "Breaking up mons " << MonIndices[0] << " , " << MonIndices[1] << std::endl;
            //for (int i=0; i < 6; i++){
            //    //std::cout << Xfib[i] << " , ";
            //}
            ////std::cout << std::endl;
            vec3 FixedMon = _Fibers[iFib].getPointedEnd();
            int FixedIndex = MonIndices[0];
            int MovedIndex = MonIndices[1];
            double re = unifdist(rngu);
            if (re < 0.5){
                FixedMon = _Fibers[iFib].getBarbedEnd();
                FixedIndex = MonIndices[1];
                MovedIndex = MonIndices[0];
            }
            // Birth another monomer within the reactive sphere
            vec3 deltaX = PointInSphere(_BindingRadius); 
            ////std::cout << "Delta X = " << deltaX[0] << " , " << deltaX[1] << " , " << deltaX[2] << std::endl;
            vec3 NewMon;
            for (int d =0; d < 3; d++){
                NewMon[d] = FixedMon[d]+deltaX[d]; 
            }
            ////std::cout << "New monomer " << FixedIndex << " = " << FixedMon[0] << " , " << FixedMon[1] << " , " << FixedMon[2] << std::endl;
            ////std::cout << "New monomer " << MovedIndex << " = " << NewMon[0] << " , " << NewMon[1] << " , " << NewMon[2] << std::endl;
            // Reset structure index and barbed/pointed ends
            _StructureIndex[FixedIndex]=-1;
            _BarbedEnd[FixedIndex] = false;
            _PointedEnd[FixedIndex] = false;
            _StructureIndex[MovedIndex]=-1;
            _BarbedEnd[MovedIndex] = false;
            _PointedEnd[MovedIndex] = false;
            _Monomers[MovedIndex].setX0(NewMon);
            _Monomers[FixedIndex].setX0(FixedMon);
            //std::cout << "Removing fiber " << iFib << std::endl;
            _Fibers.erase(_Fibers.begin() + iFib);   
            // Update structure factors ahead of it
            for (uint iMon=0; iMon < _TotalMonomers; iMon++){
                if (_StructureIndex[iMon] > iFib){
                    _StructureIndex[iMon]--;
                }       
            }
        }
        
        void NMonomerFibBreakup(const int iFib, double dt){
            // Start with barbed end
            double BarbedOffProb = _BarbedUnbindingRate*dt;
            int nMonomers = _Fibers[iFib].NumMonomers();
            if (nMonomers == 2){
                BarbedOffProb = 0.5*_TwoMonOffRate*dt;
            }
            //std::cout << "Barbed end off with probability " << BarbedOffProb << std::endl;
            double r = unifdist(rngu);
            if (r < BarbedOffProb){
                // Break up the fiber
                int nMonomers = _Fibers[iFib].NumMonomers();
                ////std::cout << "Barbed end coming off of " << nMonomers << " monomers " << std::endl;
                if (nMonomers == 2){
                    TwoMonomerFibBreakup(iFib);
                    return;
                }
                // More than 2 monomers  
                vec3 BarbedEnd = _Fibers[iFib].getBarbedEnd();
                int MonIndex = _Fibers[iFib].BarbedIndex();
                ////std::cout << "Old barbed end #" << MonIndex << " = " << BarbedEnd[0] << " , " << BarbedEnd[1] << " , " << BarbedEnd[2] << std::endl;
                _Fibers[iFib].removeMonomer(false); // remove from barbed end
                _BarbedEnd[MonIndex] = false;
                _StructureIndex[MonIndex] = -1;
                int NewBarbedIndex = _Fibers[iFib].BarbedIndex();
                _BarbedEnd[NewBarbedIndex] = true;
                // Add monomer to random pt
                vec3 DeltaX = PointInSphere(_BindingRadius);
                BarbedEnd = _Fibers[iFib].getBarbedEnd();
                ////std::cout << "New barbed end " << NewBarbedIndex << " = " << BarbedEnd[0] << " , " << BarbedEnd[1] << " , " << BarbedEnd[2] << std::endl;
                vec3 NewPt;
                for (int d =0; d < 3; d++){
                    NewPt[d] = BarbedEnd[d]+DeltaX[d];
                }
                ////std::cout << "New pt " << MonIndex << " = " << NewPt[0] << " , " << NewPt[1] << " , " << NewPt[2] << std::endl;
                _Monomers[MonIndex].setX0(NewPt);
            }
            
            double PointedOffProb = _PointedUnbindingRate*dt;
            nMonomers = _Fibers[iFib].NumMonomers();
            if (nMonomers == 2){
                PointedOffProb = 0.5*_TwoMonOffRate*dt;
            }
            //std::cout << "Pointed end off with probability " << PointedOffProb << std::endl;
            r = unifdist(rngu);
            if (r < PointedOffProb){
                // Break up the fiber
                int nMonomers = _Fibers[iFib].NumMonomers();
                ////std::cout << "Pointed end comining off of " << nMonomers << " monomers " << std::endl;
                if (nMonomers == 2){
                    TwoMonomerFibBreakup(iFib);
                    return;
                }
                // More than 2 monomers  
                vec3 PointedEnd = _Fibers[iFib].getPointedEnd();
                int MonIndex = _Fibers[iFib].PointedIndex();
                ////std::cout << "Old pointed end #" << MonIndex << " = " << PointedEnd[0] << " , " << PointedEnd[1] << " , " << PointedEnd[2] << std::endl;
                _Fibers[iFib].removeMonomer(true); // remove from barbed end
                _PointedEnd[MonIndex] = false;
                _StructureIndex[MonIndex] = -1;
                int NewPointedIndex = _Fibers[iFib].PointedIndex();
                _PointedEnd[NewPointedIndex] = true;
                // Add monomer to random pt
                vec3 DeltaX = PointInSphere(_BindingRadius);
                PointedEnd = _Fibers[iFib].getPointedEnd();
                ////std::cout << "New pointed end #" << NewPointedIndex << " = " << PointedEnd[0] << " , " << PointedEnd[1] << " , " << PointedEnd[2] << std::endl;
                vec3 NewPt;
                for (int d =0; d < 3; d++){
                    NewPt[d] = PointedEnd[d]+DeltaX[d];
                }
                ////std::cout << "New pt #" << MonIndex << " = " << NewPt[0] << " , " << NewPt[1] << " , " << NewPt[2] << std::endl;
                _Monomers[MonIndex].setX0(NewPt);
            }
            ////std::cout << "End method " << std::endl;
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

PYBIND11_MODULE(AllTheActin, m) {
    py::class_<AllTheActin>(m, "AllTheActin")
        .def(py::init<uint, intvec, vec3, double, double, double, vec,int,int>())
        .def("Diffuse",&AllTheActin::Diffuse)
        .def("React",&AllTheActin::React)
        .def("getX", &AllTheActin::getX)
        .def("getStructureIDs",&AllTheActin::getStructureIDs);
}    
