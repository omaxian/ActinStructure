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
    
    AllTheActin(npDoub pyPositions, vec3 Lengths, double a, double kbT, double mu, int seed)
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
        double rxnradius = 4*_a;
        int nThr = 1;
        _ReactionNeighbors = NeighborSearcher(_TotalMonomers, Lengths, rxnradius, nThr);
        
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
    
    npInt React(double dt){
        AssembleX();
        intvec NeighborList = _ReactionNeighbors.NeighborList(_X);
        uint nPairs = NeighborList.size()/2;
        // Here is where you would make a set of rates and assemble an array of rates, etc.
        // For now we are going to just assume that any possible reaction happens (diffusion limited)
        int newStructIndex =  _StructureIndex[_TotalMonomers-1]+1;
        intvec MonomersToDelete;
        for (uint iPair=0; iPair < nPairs; iPair++){
            int pt1 = NeighborList[2*iPair];
            int pt2 = NeighborList[2*iPair+1];
            int struct1 = _StructureIndex[pt1];
            int struct2 = _StructureIndex[pt2];
            // Here you would check if a reaction could actually occur
            if (struct1 >=0 && struct2 >=0){
                // Do nothing (two fibers can't come together)
                //std::cout << "Two fibers; doing nothing" << std::endl;
            } else if (struct1 == -1 && struct2 == -1){
                // 2 monomers forming a fiber
                std::cout << "Two monomers; forming fiber number " << newStructIndex << std::endl;
                // Randomly switch the barbed end
                double r = unifdist(rngu);
                if (r < 0.5){
                   int temp = pt1;
                   pt1 = pt2; 
                   pt2 = temp;
                }
                std::cout << "Barbed end at pt " << pt2 << std::endl;
                MonomersToDelete.push_back(pt1);
                MonomersToDelete.push_back(pt2);
                vec X0 = _Monomers[pt1].getX();
                vec X1 = _Monomers[pt2].getX();
                vec3 tau, X0Dim;
                for (int d = 0; d < 3; d++){
                    tau[d] = X1[d]-X0[d];
                    X0Dim[d] = X0[d];
                }
                normalize(tau);
                _Fibers.push_back(Fiber(X0Dim,tau,2,2*_a,_a,_mu,_kbT));
                _StructureIndex[pt1] = newStructIndex;
                _StructureIndex[pt2] = newStructIndex;
                _BarbedEnd[pt2] = true;
                newStructIndex++;
            } else if (_BarbedEnd[pt1] || _BarbedEnd[pt2]){ // one fiber and one monomer
                // One is a monomer, one is a fiber. Add the monomer to the 
                // barbed end. Need to check if the monomer is actually at the barbed end
                // Identify the barbed end
                std::cout << "Adding to existing fiber " << std::endl;
                int BEndIndex = pt1;
                int FreeIndex = pt2;
                if (_BarbedEnd[pt2]){
                    BEndIndex = pt2;
                    FreeIndex = pt1;
                }
                MonomersToDelete.push_back(FreeIndex); 
                int ExistingFibIndex = _StructureIndex[BEndIndex]; // The structure to join
                std::cout << "Barbed end " << BEndIndex << " adding monomer # " << FreeIndex << " to fiber " << ExistingFibIndex << std::endl;
                _Fibers[ExistingFibIndex].addMonomer();
                _StructureIndex[FreeIndex] = ExistingFibIndex;
                _BarbedEnd[BEndIndex]=false;
                _BarbedEnd[FreeIndex]=true;
            } else {
                std::cout << "Rejected because monomer is not close to BARBED end" << std::endl;
                std::cout << "Points " << pt1 << " , " << pt2 << std::endl;
                std::cout << "Structures " << struct1 << " , " << struct2 << std::endl;
                std::cout << "Barbed ends " << _BarbedEnd[pt1]  << " , " << _BarbedEnd[pt2]  << std::endl;
            }
        }
        std::sort(MonomersToDelete.begin(), MonomersToDelete.end(), std::greater<int>());
        for (uint i=0; i < MonomersToDelete.size(); i++){
            int indexToDelete = MonomersToDelete[i];
            std::cout << "Deleting index " << indexToDelete << " from monomer list." << std::endl;
            _Monomers.erase(_Monomers.begin() + indexToDelete);
        }
        std::cout << "Number of free monomers " << _Monomers.size() << std::endl;
        return makePyArray(NeighborList);
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
        intvec _StructureIndex;
        std::vector<bool>  _BarbedEnd;
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
                }  
                StructIndex++;
                StartMon+=nMon;
            }
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
        .def(py::init<npDoub, vec3, double, double, double, int>())
        .def("Diffuse",&AllTheActin::Diffuse)
        .def("React",&AllTheActin::React)
        .def("getX", &AllTheActin::getX)
        .def("getStructureIDs",&AllTheActin::getStructureIDs);
}    
