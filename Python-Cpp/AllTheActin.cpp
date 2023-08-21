#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "ActinStructure.cpp"

namespace py=pybind11;
typedef py::array_t<int, py::array::c_style | py::array::forcecast> npInt; 
typedef py::array_t<double, py::array::c_style | py::array::forcecast> npDoub; 

class AllTheActin{
    
    public:
    
    AllTheActin(){
        // A fiber
        _TotalMonomers = 28;
        vec3 X0 = {0,0,0};
        int nForFiber = 10;
        double a = 4e-3;
        double kbT = 4.1e-3;
        double mu = 1.0;
        _Fibers.push_back(Fiber(X0,nForFiber,2*a,a,mu,kbT));
        // A branched fiber
        int nForBranch = 18;
        int nLinFib = 3;
        intvec BranchStartIndex = {0,10,15,nForBranch};
        intvec AttachPoints = {0, 4, 12};
        _BranchedFibers.push_back(BranchedFiber(X0,nForBranch,2*a,a,mu,kbT,nLinFib,BranchStartIndex,AttachPoints));
        _X = vec(3*_TotalMonomers);
    }
    
    void Diffuse(double dt, npDoub pyRandVec, npDoub pyRandDrift, bool DoDrift){
        vec RandVec(pyRandVec.size());
        std::memcpy(RandVec.data(),pyRandVec.data(),pyRandVec.size()*sizeof(double));  
        vec RandDriftVec(pyRandDrift.size());
        std::memcpy(RandDriftVec.data(),pyRandDrift.data(),pyRandDrift.size()*sizeof(double));  
        int start=0; 
        // Diffuse the fibers
        for (int i=0; i < _Fibers.size(); i++){
            int nRand = _Fibers[i].NumberRand();
            vec RandomNumbers(nRand);
            vec RandomDrift(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=RandVec[start+j];
                RandomDrift[j]=RandDriftVec[start+j];
            }
            start+=nRand;    
            _Fibers[i].Diffuse(dt,RandomNumbers,DoDrift,RandomDrift);
        }
        // Diffuse the branched fibers
        for (int i=0; i < _BranchedFibers.size(); i++){
            int nRand = _BranchedFibers[i].NumberRand();
            vec RandomNumbers(nRand);
            vec RandomDrift(nRand);
            for (int j=0; j < nRand; j++){
                RandomNumbers[j]=RandVec[start+j];
                RandomDrift[j]=RandDriftVec[start+j];
            }
            start+=nRand;    
            _BranchedFibers[i].Diffuse(dt,RandomNumbers,DoDrift,RandomDrift);
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
        double _TotalMonomers;
        vec _X; 
        
        void AssembleX(){
            int start=0;
            for (int i=0; i < _Fibers.size(); i++){
                vec XStruct = _Fibers[i].getX();
                for (int j=0; j < XStruct.size(); j++){
                    _X[start+j]=XStruct[j];
                }
                start+=XStruct.size();
            }  
            for (int i=0; i < _BranchedFibers.size(); i++){
                vec XStruct = _BranchedFibers[i].getX();
                for (int j=0; j < XStruct.size(); j++){
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
        .def(py::init<>())
        .def("Diffuse",&AllTheActin::Diffuse)
        .def("getX", &AllTheActin::getX);
}    
