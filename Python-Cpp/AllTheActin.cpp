#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ActinStructure.cpp"

namespace py=pybind11;
typedef py::array_t<int, py::array::c_style | py::array::forcecast> npInt; 
typedef py::array_t<double, py::array::c_style | py::array::forcecast> npDoub; 

class AllTheActin{
    
    public:
    
    AllTheActin(npDoub pyPositions, double a, double kbT, double mu){
        // Initialize with only monomers
        vec Positions(pyPositions.size());
        std::memcpy(Positions.data(),pyPositions.data(),pyPositions.size()*sizeof(double));  
        _TotalMonomers = Positions.size()/3;
        _a = a;
        _kbT = kbT;
        _mu = mu;
        _X = vec(3*_TotalMonomers);
        // Initialize the monomers 
        for (uint i=0; i < _TotalMonomers; i++){
            vec3 MonPt;
            for (int d=0; d<3; d++){
               MonPt[d]=Positions[3*i+d];
            }
            _Monomers.push_back(Monomer(MonPt,a,mu,kbT));
        }
    }
    
    void Diffuse(double dt, npDoub pyRandVec){
        vec RandVec(pyRandVec.size());
        std::memcpy(RandVec.data(),pyRandVec.data(),pyRandVec.size()*sizeof(double));  
        int start=0; 
        // Diffuse monomers
        for (uint i=0; i < _Monomers.size(); i++){
            int nRand = _Monomers[i].NumberRand();
            vec RandomNumbers(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=RandVec[start+j];
            }
            start+=nRand;    
            _Monomers[i].Diffuse(dt,RandomNumbers);
        }
        // Diffuse the fibers
        for (uint i=0; i < _Fibers.size(); i++){
            int nRand = _Fibers[i].NumberRand();
            vec RandomNumbers(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=RandVec[start+j];
            }
            start+=nRand;    
            _Fibers[i].Diffuse(dt,RandomNumbers);
        }
        // Diffuse the branched fibers
        for (uint i=0; i < _BranchedFibers.size(); i++){
            int nRand = _BranchedFibers[i].NumberRand();
            vec RandomNumbers(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=RandVec[start+j];
            }
            start+=nRand;    
            _BranchedFibers[i].Diffuse(dt,RandomNumbers);
        }
    }
    
    npDoub getX(){
        AssembleX();
        return makePyDoubleArray(_X);
    }
           
        
    private:
        std::vector <Monomer> _Monomers;
        std::vector <Fiber> _Fibers;
        std::vector <BranchedFiber> _BranchedFibers;
        uint _TotalMonomers;
        double _a, _kbT, _mu;
        vec _X; 
        
        void AssembleX(){
            uint start=0;
            for (uint i=0; i < _Monomers.size(); i++){
                vec XStruct = _Monomers[i].getX();
                for (uint j=0; j < XStruct.size(); j++){
                    _X[start+j]=XStruct[j];
                }
                start+=XStruct.size();
            }
            for (uint i=0; i < _Fibers.size(); i++){
                vec XStruct = _Fibers[i].getX();
                for (uint j=0; j < XStruct.size(); j++){
                    _X[start+j]=XStruct[j];
                }
                start+=XStruct.size();
            }  
            for (uint i=0; i < _BranchedFibers.size(); i++){
                vec XStruct = _BranchedFibers[i].getX();
                for (uint j=0; j < XStruct.size(); j++){
                    _X[start+j]=XStruct[j];
                }
                start+=XStruct.size();
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

};

PYBIND11_MODULE(AllTheActin, m) {
    py::class_<AllTheActin>(m, "AllTheActin")
        .def(py::init<npDoub, double, double, double>())
        .def("Diffuse",&AllTheActin::Diffuse)
        .def("getX", &AllTheActin::getX);
}    
