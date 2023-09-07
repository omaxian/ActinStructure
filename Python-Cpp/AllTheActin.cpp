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
    
    AllTheActin(npDoub pyPositions, vec3 Lengths, double a, double kbT, double mu, vec RxnRates, int seed)
    :_ReactionNeighbors(){
        // Initialize with only monomers
        vec Positions(pyPositions.size());
        std::memcpy(Positions.data(),pyPositions.data(),pyPositions.size()*sizeof(double));  
        _TotalMonomers = Positions.size()/3;
        for (uint i=0; i < _TotalMonomers; i++){
            vec3 MonPt;
            for (int d=0; d<3; d++){
               MonPt[d]=Positions[3*i+d];
            }
            _Monomers.push_back(Monomer(MonPt,a,mu,kbT));
        }
        _StructureIndex = intvec(_TotalMonomers,-1);
        _BarbedEnd = std::vector<bool>(_TotalMonomers,false);
        _PointedEnd = std::vector<bool>(_TotalMonomers,false);
        
        // Diffusion variables
        _a = a;
        _kbT = kbT;
        _mu = mu;
        _X = vec(3*_TotalMonomers);
        rng.seed(seed);
        normaldist = std::normal_distribution<double>(0.0,1.0);
        
        // Reaction variables
        rngu.seed(seed+1);
        unifdist = std::uniform_real_distribution<double>(0.0,1.0);
        _BindingRadius = RxnRates[0];
        _TwoMonRate = RxnRates[1];
        _BarbedBindingRate = RxnRates[2];
        _BarbedUnbindingRate = RxnRates[3];
        _PointedBindingRate = RxnRates[4];
        _PointedUnbindingRate = RxnRates[5];
        int nThr = 1;
        _ReactionNeighbors = NeighborSearcher(_TotalMonomers, Lengths, _BindingRadius, nThr);
        
    }
    
    void Diffuse(double dt){
        int start=0; 
        // Diffuse monomers
        for (uint i=0; i < _Monomers.size(); i++){
            int nRand = _Monomers[i].NumberRand();
            vec RandomNumbers(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=normaldist(rng);
            }
            start+=nRand;    
            _Monomers[i].Diffuse(dt,RandomNumbers);
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
         
        intvec MonomersToDelete;
        int newStructIndex =  _StructureIndex[_TotalMonomers-1]+1;
        
        // Binding reactions
        for (uint iPair = 0; iPair < nPair; iPair++){
            int pt1 = NeighborList[2*iPair];
            int pt2 = NeighborList[2*iPair+1];
            int struct1 = _StructureIndex[pt1];
            int struct2 = _StructureIndex[pt2];
            if (struct1 >=0 && struct2 >=0){
                // Do nothing (two fibers can't come together)
                //std::cout << "Two fibers; doing nothing" << std::endl;
            } else if (struct1 == -1 && struct2 == -1){
                // 2 monomers forming a fiber
                bool Reacted = FormFiberFromMonomers(pt1, pt2, dt, MonomersToDelete,newStructIndex);
                if (Reacted){
                    //std::cout << "Add mon+mon, New structure index " << newStructIndex << std::endl;
                }
            } else { // There is one monomer and a fiber
                bool Reacted = AddMonomerToFiber(struct1, pt1, struct2, pt2, dt, MonomersToDelete);
                if (Reacted){
                    //std::cout << "Add mon+fib, New structure index " << newStructIndex << std::endl;
                }
            } 
        }
        std::sort(MonomersToDelete.begin(), MonomersToDelete.end(), std::greater<int>());
        for (uint i=0; i < MonomersToDelete.size(); i++){
            int indexToDelete = MonomersToDelete[i];
            //std::cout << "Deleting index " << indexToDelete << " from monomer list." << std::endl;
            _Monomers.erase(_Monomers.begin() + indexToDelete);
        }
        //std::cout << "Number of free monomers " << _Monomers.size() << std::endl;
        /*vec X = _Fibers[0].getX();
        for (int i=0; i < _Fibers[0].NumMonomers(); i++){
            std:: cout << X[3*i] << " " << X[3*i+1] << " " << X[3*i+2] << std::endl;
        }*/
        
        uint nFib = _Fibers.size();
        // Process unbinding reactions
        for (int iFib = nFib-1; iFib >= 0 ; iFib--){
            NMonomerFibBreakup(iFib, dt);
        }
        
        //std::cout << "Number monomers " << _Monomers.size() << std::endl;
        //std::cout << "Number fibers " << _Fibers.size() << std::endl;
        return _Monomers.size();
    }
    
    npDoub getX(){
        AssembleX();
        return makePyDoubleArray(_X);
    }
    
    npInt getStructureIDs(){
        AssembleX();
        return makePyArray(_StructureIndex);
    }
           
        
    private:
        std::vector <Monomer> _Monomers;
        std::vector <Fiber> _Fibers;
        std::vector <BranchedFiber> _BranchedFibers;
        uint _TotalMonomers;
        double _a, _kbT, _mu;
        double _BindingRadius, _TwoMonRate, _BarbedBindingRate, _BarbedUnbindingRate, _PointedBindingRate, _PointedUnbindingRate;
        intvec _StructureIndex;
        std::vector<bool>  _BarbedEnd, _PointedEnd;
        vec _X; 
        std::normal_distribution<double> normaldist;
        std::mt19937_64 rng;
        std::uniform_real_distribution<double> unifdist;
        std::mt19937_64 rngu;
        NeighborSearcher _ReactionNeighbors;
        
        
        void AssembleX(){
            uint start=0;
            for (uint i=0; i < _Monomers.size(); i++){
                vec XStruct = _Monomers[i].getX();
                for (uint j=0; j < XStruct.size(); j++){
                    _X[start+j]=XStruct[j];
                }
                start+=XStruct.size();
                _StructureIndex[i] = -1;
                _BarbedEnd[i] = false;
                _PointedEnd[i] = false;
            }
            int StructIndex=0;
            int StartMon = _Monomers.size();
            for (uint i=0; i < _Fibers.size(); i++){
                vec XStruct = _Fibers[i].getX();
                for (uint j=0; j < XStruct.size(); j++){
                    _X[start+j]=XStruct[j];
                }
                start+=XStruct.size();
                uint nMon = XStruct.size()/3;
                for (uint iMon=0; iMon < nMon; iMon++){
                    _StructureIndex[StartMon+iMon] = StructIndex;
                    _BarbedEnd[StartMon+iMon] = _Fibers[i].isBarbedEnd(iMon);
                    _PointedEnd[StartMon+iMon] = _Fibers[i].isPointedEnd(iMon);
                }  
                StructIndex++;
                StartMon+=nMon;
            }  
            for (uint i=0; i < _BranchedFibers.size(); i++){
                vec XStruct = _BranchedFibers[i].getX();
                for (uint j=0; j < XStruct.size(); j++){
                    _X[start+j]=XStruct[j];
                }
                start+=XStruct.size();
                uint nMon = XStruct.size()/3;
                for (uint iMon=0; iMon < nMon; iMon++){
                    _StructureIndex[StartMon+iMon] = StructIndex;
                    _BarbedEnd[StartMon+iMon] = _BranchedFibers[i].isBarbedEnd(iMon);
                    _PointedEnd[StartMon+iMon] = _BranchedFibers[i].isPointedEnd(iMon);
                }  
                StructIndex++;
                StartMon+=nMon;
            }
        } 
        
        bool FormFiberFromMonomers(uint pt1, uint pt2, double dt, intvec &MonomersToDelete, int &newStructIndex){
            double pForm = dt*_TwoMonRate;
            //std::cout << "Two monomers; forming fiber number " << newStructIndex << " with probability " << pForm << std::endl;
            double r = unifdist(rngu);
            if (r < pForm) { 
                //std::cout << "Random val " << r << " -> reaction proceeding" << std::endl;
                // Randomly switch the pointed end
                double re = unifdist(rngu);
                if (re < 0.5){
                   int temp = pt1;
                   pt1 = pt2; 
                   pt2 = temp;
                }
                // The pointed end is pt1
                //std::cout << "Pointed end at pt " << pt1 << std::endl;
                MonomersToDelete.push_back(pt1);
                MonomersToDelete.push_back(pt2);
                vec X0 = _Monomers[pt1].getX();
                vec3 tau, X0Dim;
                for (int d = 0; d < 3; d++){
                    tau[d] = normaldist(rng);
                    X0Dim[d] = X0[d];
                }
                normalize(tau);
                _Fibers.push_back(Fiber(X0Dim,tau,2,2*_a,_a,_mu,_kbT));
                _StructureIndex[pt1] = newStructIndex;
                _StructureIndex[pt2] = newStructIndex;
                _BarbedEnd[pt2] = true;
                //std::cout << "Barbed end at point " << pt2 << std::endl;
                _PointedEnd[pt1] = true;
                //std::cout << "Pointed end at point " << pt1 << std::endl;
                newStructIndex++;
                return true;
            }
            return false;
        }
        
        bool AddMonomerToFiber(const int struct1, const int pt1, const int struct2, const int pt2, double dt, intvec &MonomersToDelete){
            // Find the monomer
            //std::cout << struct1 << " , " << pt1 << " , " << struct2 << "," << pt2 << std::endl;
            int MonIndex = pt2;
            int FibIndex = pt1;
            if (struct1 == -1){
                MonIndex = pt1;
                FibIndex = pt2;
            }
            //std::cout << "The fiber is at point " << FibIndex << std::endl;
            double RxnRate = 0;
            if (_BarbedEnd[FibIndex]){
                RxnRate = _BarbedBindingRate;  
                //std::cout << "Barbed end at point " << FibIndex << std::endl;
            } else if (_PointedEnd[FibIndex]){
                RxnRate = _PointedBindingRate;
                //std::cout << "Pointed end at point " << FibIndex << std::endl;
            }  
            double pForm = dt*RxnRate;
            //std::cout << "Adding monomer to fiber with probability " << pForm << std::endl;
            double r = unifdist(rngu);
            if (r < pForm) { 
                // Process the reaction
                //std::cout << "Adding to existing fiber " << std::endl;
                MonomersToDelete.push_back(MonIndex); 
                int ExistingFibIndex = _StructureIndex[FibIndex]; // The structure to join
                //std::cout << "Point " << FibIndex << " adding monomer # " << MonIndex << " to fiber " << ExistingFibIndex << std::endl;
                _Fibers[ExistingFibIndex].addMonomer(_PointedEnd[FibIndex]);
                _StructureIndex[MonIndex] = ExistingFibIndex;
                if (_BarbedEnd[FibIndex]){
                    //std::cout << "New barbed end" << std::endl;
                    _BarbedEnd[FibIndex]=false;
                    _BarbedEnd[MonIndex]=true;
                } else if (_PointedEnd[FibIndex]){
                    //std::cout << "New pointed end" << std::endl;
                    _PointedEnd[FibIndex]=false;
                    _PointedEnd[MonIndex]=true;
                }
                return true;
            }
            return false;
        }
        
        void TwoMonomerFibBreakup(const int iFib){
            //std::cout << "Breaking up fiber " << std::endl;
            vec Xfib = _Fibers[iFib].getX();
            //for (int i=0; i < 6; i++){
                //std::cout << Xfib[i] << " , ";
            //}
            //std::cout << std::endl;
            vec3 FixedMon = _Fibers[iFib].getPointedEnd();
            double re = unifdist(rngu);
            if (re < 0.5){
                FixedMon = _Fibers[iFib].getBarbedEnd();
            }
            // Birth another monomer within the reactive sphere
            vec3 deltaX = PointInSphere(_BindingRadius); 
            //std::cout << "Delta X = " << deltaX[0] << " , " << deltaX[1] << " , " << deltaX[2] << std::endl;
            vec3 NewMon;
            for (int d =0; d < 3; d++){
                NewMon[d] = FixedMon[d]+deltaX[d]; 
            }
            //std::cout << "New monomer = " << FixedMon[0] << " , " << FixedMon[1] << " , " << FixedMon[2] << std::endl;
            //std::cout << "New monomer = " << NewMon[0] << " , " << NewMon[1] << " , " << NewMon[2] << std::endl;
            _Fibers.erase(_Fibers.begin() + iFib);    
            _Monomers.push_back(Monomer(FixedMon,_a,_mu,_kbT));
            _Monomers.push_back(Monomer(NewMon,_a,_mu,_kbT));
        }
        
        void NMonomerFibBreakup(const int iFib, double dt){
            // Start with barbed end
            double BarbedOffProb = _BarbedUnbindingRate*dt;
            //std::cout << "Barbed end off with probability " << BarbedOffProb << std::endl;
            double r = unifdist(rngu);
            if (r < BarbedOffProb){
                // Break up the fiber
                int nMonomers = _Fibers[iFib].NumMonomers();
                //std::cout << "Barbed end coming off of " << nMonomers << " monomers " << std::endl;
                if (nMonomers == 2){
                    TwoMonomerFibBreakup(iFib);
                    return;
                }
                // More than 2 monomers  
                vec3 BarbedEnd = _Fibers[iFib].getBarbedEnd();
                //std::cout << "Old barbed end " << BarbedEnd[0] << " , " << BarbedEnd[1] << " , " << BarbedEnd[2] << std::endl;
                _Fibers[iFib].removeMonomer(false); // remove from barbed end
                // Add monomer to random pt
                vec3 DeltaX = PointInSphere(_BindingRadius);
                BarbedEnd = _Fibers[iFib].getBarbedEnd();
                //std::cout << "New barbed end " << BarbedEnd[0] << " , " << BarbedEnd[1] << " , " << BarbedEnd[2] << std::endl;
                vec3 NewPt;
                for (int d =0; d < 3; d++){
                    NewPt[d] = BarbedEnd[d]+DeltaX[d];
                }
                //std::cout << "New pt " << NewPt[0] << " , " << NewPt[1] << " , " << NewPt[2] << std::endl;
                _Monomers.push_back(Monomer(NewPt,_a,_mu,_kbT));
            }
            
            double PointedOffProb = _PointedUnbindingRate*dt;
            //std::cout << "Pointed end off with probability " << PointedOffProb << std::endl;
            r = unifdist(rngu);
            if (r < PointedOffProb){
                // Break up the fiber
                int nMonomers = _Fibers[iFib].NumMonomers();
                //std::cout << "Pointed end comining off of " << nMonomers << " monomers " << std::endl;
                if (nMonomers == 2){
                    TwoMonomerFibBreakup(iFib);
                    return;
                }
                // More than 2 monomers  
                vec3 PointedEnd = _Fibers[iFib].getPointedEnd();
                //std::cout << "Removed monomer from fiber " << nMonomers << " monomers " << std::endl;
                //std::cout << "Old pointed end " << PointedEnd[0] << " , " << PointedEnd[1] << " , " << PointedEnd[2] << std::endl;
                _Fibers[iFib].removeMonomer(true); // remove from barbed end
                // Add monomer to random pt
                vec3 DeltaX = PointInSphere(_BindingRadius);
                PointedEnd = _Fibers[iFib].getPointedEnd();
                //std::cout << "New pointed end " << PointedEnd[0] << " , " << PointedEnd[1] << " , " << PointedEnd[2] << std::endl;
                vec3 NewPt;
                for (int d =0; d < 3; d++){
                    NewPt[d] = PointedEnd[d]+DeltaX[d];
                }
                //std::cout << "New pt " << NewPt[0] << " , " << NewPt[1] << " , " << NewPt[2] << std::endl;
                _Monomers.push_back(Monomer(NewPt,_a,_mu,_kbT));
            }
            //std::cout << "End method " << std::endl;
        }  
        
        vec3 PointInSphere(double r){
            double u = unifdist(rngu);
            vec3 tau;
            for (int d=0; d < 3; d++){
                tau[d] = normaldist(rng);
            }
            normalize(tau);
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
        .def(py::init<npDoub, vec3, double, double, double, vec,int>())
        .def("Diffuse",&AllTheActin::Diffuse)
        .def("React",&AllTheActin::React)
        .def("getX", &AllTheActin::getX)
        .def("getStructureIDs",&AllTheActin::getStructureIDs);
}    
