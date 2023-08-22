#include <stdio.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "VectorMethods.cpp"
#include "types.h"
#include <random>

/**
NeighborSearch.cpp
**/

namespace py=pybind11;
typedef py::array_t<int, py::array::c_style | py::array::forcecast> npInt; 
typedef py::array_t<double, py::array::c_style | py::array::forcecast> npDoub; 

class NeighborSearcher {
    
    public:
    
    NeighborSearcher(){
    }
    
    NeighborSearcher(uint nPoints, vec3 Lengths, double cutoff, int nThr){
        _cutoff = cutoff;
        _nPts = nPoints;
        std::memcpy(_Lens.data(),Lengths.data(),Lengths.size()*sizeof(double));
        _nBins = intvec(3);
        int maxBins = 20;
        for (int d=0; d < 3; d++){  
            _nBins[d] = int(Lengths[d]/_cutoff);
            if (_nBins[d] > maxBins){
                _nBins[d] = maxBins;
            }
            if (_nBins[d] < 3){
                std::cout << "Neighbor search will give duplicates because # bins < 3" << std::endl;
            }
            _DPerBin[d] = Lengths[d]/_nBins[d];
        }
        std::cout << "Number of bins " << _nBins[0] << " , " << _nBins[1] << " , " << _nBins[2] << std::endl;
        _FirstByBin = intvec(_nBins[0]*_nBins[1]*_nBins[2]);
        _NextByPt = intvec(_nPts);
        _BinEachPt = intvec(_nPts);
        _X = vec(3*_nPts);
        _nThr = nThr;
    }
    
    intvec NeighborList(const vec &X){
        BinPoints(X);
        intvec AllNeighbors;
        // Loop over bins
        #pragma omp parallel num_threads(_nThr)
        {
        intvec privateNeighbors;
        #pragma omp for nowait
        for (uint iB=0; iB < _FirstByBin.size(); iB++){
            // Find Neighbor bins
            intvec Neighbors = NeighborBins(iB);
            // Loop over points in the box and then over neighbor bins
            int iPt = _FirstByBin[iB];
            while (iPt!=-1){ // loop over points in box
                // Look at neighbor bins
                for (uint jB = 0; jB < Neighbors.size(); jB++){
                    int NeighborNumber =  Neighbors[jB];
                    int jPt = _FirstByBin[NeighborNumber];
                    while (jPt!=-1){
                        if (iPt < jPt){
                            vec3 disp;
                            for (int d=0; d < 3; d++){
                                disp[d] = X[3*iPt+d]-X[3*jPt+d];
                                double shift = std::round(disp[d]/_Lens[d]);
                                disp[d] -= _Lens[d]*shift;
                            }
                            double r = normalize(disp);
                            if (r < _cutoff){
                                privateNeighbors.push_back(iPt);
                                privateNeighbors.push_back(jPt);
                            }
                        }
                        jPt = _NextByPt[jPt];
                    }
                }
                iPt = _NextByPt[iPt];
            }
        }
        #pragma omp critical
        AllNeighbors.insert(AllNeighbors.end(), privateNeighbors.begin(), privateNeighbors.end());
        }
        return AllNeighbors;
    }
    
    // Python wrapper for debugging     
    npInt pyNeighborList(npDoub pyX){
        vec X(pyX.size());
        std::memcpy(X.data(), pyX.data(),pyX.size()*sizeof(double));
        intvec Neighbors = NeighborList(X);
        return makePyArray(Neighbors);
    }
                    
    private:
    
    int _nPts, _nThr;
    vec _X;
    vec3 _Lens, _DPerBin;
    intvec _BinEachPt,_FirstByBin, _NextByPt;
    double _cutoff;
    intvec _nBins;
    
    void BinPoints(const vec &X){
        // Reset vectors
        std::fill(_FirstByBin.begin(), _FirstByBin.end(), -1);
        std::fill(_NextByPt.begin(), _NextByPt.end(), -1);
        // Bin the points 
        for (int i=0; i < _nPts; i++){
            vec3 xog;
            intvec BinNums(3);
            for (int d=0; d < 3; d++){
                xog[d] = X[3*i+d];
            }
            vec3 XInBox = ShiftX(xog);
            for (int d=0; d < 3; d++){    
                BinNums[d] = int(XInBox[d]/_DPerBin[d]);
                BinNums[d] = positive_modulo(BinNums[d],_nBins[d]); // takes care of rounding issues
            }
            _BinEachPt[i]=BinIndex(BinNums);
        }
        // Create lists first and next
        for (int i=0; i < _nPts; i++){
            int bin = _BinEachPt[i];
            if (_FirstByBin[bin] == -1){
                _FirstByBin[bin] = i;
            } else {
                int jPt = _FirstByBin[bin];
                int jPtprev = jPt;
                while (jPt!=-1){
                    jPtprev = jPt;
                    jPt = _NextByPt[jPt];  
                }
                _NextByPt[jPtprev]=i;
            }
        }
    }  
    
    vec3 ShiftX(const vec3 xog){
        vec3 XInBox;
        for (int d=0; d< 3; d++){
            XInBox[d] = positive_modulo(xog[d],_Lens[d]);
            if (XInBox[d] < 0 || XInBox[d] >=_Lens[d]){
                std::cout << "Bad x " << XInBox[d] << std::endl;
            }
        } 
        return XInBox;
    }   
    
    int BinIndex(const intvec &ThreeBin){
        return ThreeBin[0]+_nBins[0]*ThreeBin[1]+_nBins[0]*_nBins[1]*ThreeBin[2];
    }
    
    intvec ThreeBinIndex(int BinIndexIn){
        intvec ThreeBin(3);
        ThreeBin[0] = positive_modulo(BinIndexIn, _nBins[0]);
        int leftover = (BinIndexIn-ThreeBin[0])/_nBins[0];
        ThreeBin[1] = positive_modulo(leftover, _nBins[1]);
        int secondleftover = (leftover-ThreeBin[1])/_nBins[1];
        ThreeBin[2]=secondleftover;
        // Check 
        return ThreeBin;
    }
    
    intvec NeighborBins(int iBin){
        intvec Neighbors;
        intvec ThreeIndex = ThreeBinIndex(iBin);
        intvec ThisNeighbor(3);
        for (int xi=-1; xi < 2; xi++){
            ThisNeighbor[0] = positive_modulo(ThreeIndex[0]+xi,_nBins[0]);
            for (int yi=-1; yi < 2; yi++){
                ThisNeighbor[1] = positive_modulo(ThreeIndex[1]+yi, _nBins[1]);
                for (int zi=-1; zi < 2; zi++){
                    ThisNeighbor[2] = positive_modulo(ThreeIndex[2]+zi,_nBins[2]); 
                    int NeighborIndex = BinIndex(ThisNeighbor);
                    Neighbors.push_back(NeighborIndex);
                }
            }
        }
        return Neighbors;
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
    
    int positive_modulo(int i, int n) {
        return (i % n + n) % n;
    } 
    
    double positive_modulo(double x, double y){
        return std::fmod(std::fmod(x,y)+y,y);
    }
    
};

PYBIND11_MODULE(NeighborSearcher, m) {
    py::class_<NeighborSearcher>(m, "NeighborSearcher")
        .def(py::init<int, vec3, double,int>())
        .def("pyNeighborList",&NeighborSearcher::pyNeighborList);
}    


   
  


