#include <stdio.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "VectorMethods.cpp"
#include "types.h"
#include <random>

/**
ActinStructure.cpp

This needs to be rewritten so that we never create an array with nMonomers size.
It should just be the tangent vectors we are tracking and the length (in monomers)
**/

// Global variables for periodic ActinStructure
class ActinStructure {
    
    public:
    
    ActinStructure(){

    }
    
    int NumberRand(){
        return 0;
    }
    
    void Diffuse(double dt, const vec &W){
    }
    
    vec getX(){
        return _X;
    }
    
    protected:
        vec _X; // location
        vec3 _X0;
        double _DCoeff, _kbT, _Mobility;
    
};

class Monomer: public ActinStructure{

    public:
    
    Monomer(vec3 X, double a, double mu, double kbT){
        std::memcpy(_X0.data(),X.data(),X.size()*sizeof(double));
        _kbT = kbT;
        _Mobility = 1.0/(6*M_PI*mu*a);
        _DCoeff = 2.0*_kbT*_Mobility;
        _X = vec(3);
    }
    
    int NumberRand(){
        return 3;
    }
    
    void Diffuse(double dt,const vec &W){
        for (int d=0; d< 3; d++){
            double disp = sqrt(dt*_DCoeff)*W[d];
            _X0[d]+=disp;
        }
    }
    
    vec getX(){
        std::memcpy(_X.data(),_X0.data(),_X0.size()*sizeof(double));  
        return _X;   
    }
    
    void setX0(vec3 XIn){
        std::memcpy(_X0.data(),XIn.data(),XIn.size()*sizeof(double));      
    }
        
};

class Fiber: public ActinStructure{
    
    public:
    
    // Initialize
    Fiber(vec3 X0, vec3 tau, uint nMonomers, intvec MonomerInds, double spacing, double a, double mu, double kbT){
        _tau = vec(3);
        for (int d=0; d< 3; d++){
            _tau[d] = tau[d];
        }
        std::memcpy(_X0.data(),X0.data(),X0.size()*sizeof(double)); 
        _MonomerIndices = intvec(MonomerInds.size());
        std::memcpy(_MonomerIndices.data(),MonomerInds.data(),MonomerInds.size()*sizeof(int)); 
        _spacing = spacing;
        _nMonomers = nMonomers;
        _Mobility = 1.0/(6*M_PI*mu*a);
        _kbT = kbT;
        _X = vec(3*_nMonomers);
    }
    
    Fiber(vec3 X0, vec3 tau, uint nMonomers, double spacing, double a, double mu, double kbT, const vec &BindingRates){
        _tau = vec(3);
        for (int d=0; d< 3; d++){
            _tau[d] = tau[d];
        }
        std::memcpy(_X0.data(),X0.data(),X0.size()*sizeof(double)); 
        _spacing = spacing;
        _nMonomers = nMonomers;
        _Mobility = 1.0/(6*M_PI*mu*a);
        _kbT = kbT;
        _X = vec(3*_nMonomers);
        _ForminOn = false;
        _BaseBarbedBindRate = BindingRates[0];
        _BarbedBindRate = _BaseBarbedBindRate;
        _BarbedUnbindRate = BindingRates[1];
        _PointedBindRate = BindingRates[2];
        _PointedUnbindRate = BindingRates[3];
        _ForminEnhancement = BindingRates[4];
    }
       
    int NumberRand(){
        return 6;
    }
    
    void Diffuse(double dt,const vec &W){
        vec K(3*_nMonomers*6);
        computeX(_X0,_tau,_X);
        vec3 XCOM = calcKRigid(_X,K);
        vec3 X0FromCOM; 
        for (int d=0; d < 3; d++){
            X0FromCOM[d]=_X0[d]-XCOM[d];
        }
        vec KTK(6*6);
        BlasMatrixProduct(6, 3*_nMonomers,6, 1.0/_Mobility,0.0,K,true,K,KTK); // K'*M^(-1)*K
        // Check off diagonal blocks
        double maxOff=0;
        for (int iR=0; iR < 3; iR++){
            for (int d=0; d< 3; d++){
                maxOff=std::max(fabs(KTK[18+6*iR+d]),maxOff);
            }
        }
        if (maxOff > 1e-10){
            throw std::runtime_error("Trans-rot mobility coupling is not zero");
        }   
        bool normalize = true;
        bool half = true;
        double pinvtol= 1e-10;
        vec OmegaU(6);
        SolveWithPseudoInverse(6, 6,KTK, W, OmegaU, pinvtol,normalize,half,6); // N^(1/2)
        for (int d=0; d< 6; d++){
            OmegaU[d]*=sqrt(2*_kbT/dt); // sqrt(2kbT/dt)*N^(1/2)*W
        }             
        // Update center of mass
        vec3 Omega;
        for (int d=0; d< 3; d++){
            Omega[d]=dt*OmegaU[d];
        }
        rotateTangents(Omega,_tau);
        rotateVector(Omega,X0FromCOM);
        for (int d=0; d< 3; d++){
            _X0[d]=XCOM[d]+dt*OmegaU[3+d]+X0FromCOM[d];
        }
    }
        
       
    vec getX(){
        computeX(_X0,_tau,_X);
        return _X;
    }
    
    virtual void BindFormin(){
        _BarbedBindRate=_BaseBarbedBindRate*_ForminEnhancement;
        _ForminOn = true;
    }
    
    virtual void UnbindFormin(){
        _BarbedBindRate=_BaseBarbedBindRate;
        _ForminOn = false;
    }
    
    virtual double TotalBindingRate(){    
        double RateAdd = (_BarbedBindRate + _PointedBindRate);
        if (_nMonomers == 5){
           // std::cout << "MAX 5 MONOMERS " << std::endl;
            RateAdd = 0;
        }
        return RateAdd;
    }
    
    virtual void BindReaction(double r){
        double pPointed = _PointedBindRate/(_PointedBindRate+_BarbedBindRate);
        //std::cout << "Prob pointed end " << pPointed << std::endl;
        bool ToPointedEnd = r < pPointed;
        //std::cout << "ToPointedEnd " << ToPointedEnd  << std::endl;
        addMonomer(ToPointedEnd);
    }
    
    virtual double TotalUnbindingRate(){
        return (_BarbedUnbindRate+_PointedUnbindRate);
    }
    
    virtual void UnbindReaction(double r){
        double pPointed = _PointedUnbindRate/(_PointedUnbindRate+_BarbedUnbindRate);
        bool FromPointedEnd = r < pPointed;
        removeMonomer(FromPointedEnd);
    }
    
    virtual double TotalForminRate(int nFreeFormins, double ForminBindRate, double ForminUnbindRate){
        double ForminRate = ForminUnbindRate;
        if (!_ForminOn){
            ForminRate = ForminBindRate*nFreeFormins;
        } else if (_nMonomers < 4){ // Formin cannot unbind from filaments with less than 4 monomers
            ForminRate = 0;
        }
        return ForminRate;
    }
    
    virtual void ForminReaction(uint &nFreeFormins){
        if (_ForminOn){
            nFreeFormins++;
            UnbindFormin();
        } else {
            nFreeFormins--;
            BindFormin();
        }
    } 
        
    virtual bool isBarbedEnd(int MonIndex){
        return (MonIndex==_nMonomers-1);
    }
    
    virtual bool isPointedEnd(int MonIndex){
        return (MonIndex==0);
    }
    
    int NumMonomers(){
        return _nMonomers;
    }
    
    vec3 getPointedEnd(){
        return _X0;
    } 
       
    vec3 getBarbedEnd(){
        vec3 Be;
        for (int d =0; d < 3; d++){
            Be[d] = _X0[d]+_tau[d]*(_nMonomers-1)*_spacing;
        }  
        return Be;  
    }
    
    bool ForminBound(){
        return _ForminOn;
    }
    
    intvec getMonomerIndices(){
        return _MonomerIndices;
    }
    
    int BarbedIndex(){
        return _MonomerIndices[_nMonomers-1];
    }
    
    int PointedIndex(){
        return _MonomerIndices[0];
    }
    
    void addMonomer(bool PointedEnd, int MonomerInd){
        _nMonomers++;
        _X.resize(3*_nMonomers);
        if (PointedEnd) {// adjust the start point
            for (int d=0; d < 3; d++){
                _X0[d]-=_spacing*_tau[d];
            }
            _MonomerIndices.insert(_MonomerIndices.begin(), MonomerInd);
        } else {
            _MonomerIndices.push_back(MonomerInd);
        }
    }
    
    void addMonomer(bool PointedEnd){
        _nMonomers++;
        _X.resize(3*_nMonomers);
        if (PointedEnd) {// adjust the start point
            for (int d=0; d < 3; d++){
                _X0[d]-=_spacing*_tau[d];
            }
        }
    }
    
    void removeMonomer(bool PointedEnd, int &MonIndex){
        _nMonomers--;
        _X.resize(3*_nMonomers);
        if (PointedEnd) {// adjust the start point
            for (int d=0; d < 3; d++){
                _X0[d]+=_spacing*_tau[d];
            }
            MonIndex=_MonomerIndices[0];
            _MonomerIndices.erase(_MonomerIndices.begin());
        } else {
            MonIndex=_MonomerIndices[_nMonomers-1];
            _MonomerIndices.erase(_MonomerIndices.begin()+_nMonomers);
        }
    }
    
    void removeMonomer(bool PointedEnd){
        _nMonomers--;
        _X.resize(3*_nMonomers);
        if (PointedEnd) {// adjust the start point
            for (int d=0; d < 3; d++){
                _X0[d]+=_spacing*_tau[d];
            }
        }
    }
    
    
    protected: 
        int _nMonomers;
        double _spacing;
        vec _tau;
        intvec _MonomerIndices;
        double _BaseBarbedBindRate, _BarbedBindRate, _BarbedUnbindRate;
        double _PointedBindRate, _PointedUnbindRate,_ForminEnhancement;
        
        vec3 calcKRigid(const vec &X, vec &K){
            /*
            Kinematic matrix for rigid fibers
            */
            // First calculate cross product matrix
            double OneOverN = 1.0/_nMonomers;
            vec3 XCOM={0,0,0};
            for (int iPt=0; iPt < _nMonomers; iPt++){
                for (int d=0; d < 3; d++){
                    XCOM[d]+=X[3*iPt+d];
                }
            }
            for (int d=0; d < 3; d++){
                XCOM[d]*=OneOverN;
            }
            
            vec CPMatrix(9*_nMonomers);
            for (int iPt = 0; iPt < _nMonomers; iPt++){
                vec3 disp;
                for (int d=0; d < 3; d++){
                    disp[d]=X[3*iPt+d]-XCOM[d];
                }
                CPMatrix[9*iPt+1] = -disp[2];
                CPMatrix[9*iPt+2] =  disp[1];
                CPMatrix[9*iPt+3]= disp[2];
                CPMatrix[9*iPt+5]= -disp[0];
                CPMatrix[9*iPt+6]= -disp[1];
                CPMatrix[9*iPt+7]= disp[0];
            }
            // Add the identity part to K
            for (int iPt=0; iPt < _nMonomers; iPt++){
                for (int iD=0; iD < 3; iD++){
                    int rowstart = (3*iPt+iD)*6;
                    int Xrowstart = (3*iPt+iD)*3;
                    // Copy the row entry from XMatTimesCPMat
                    for (int iEntry=0; iEntry < 3; iEntry++){
                        K[rowstart+iEntry]=-CPMatrix[Xrowstart+iEntry];
                    }
                }
            } 
            for (int iPt=0; iPt < _nMonomers; iPt++){
                for (int iD=0; iD < 3; iD++){
                    int rowstart = (3*iPt+iD)*6;
                    K[rowstart+3+iD]=1;
                }
            }
            return XCOM;
        }
        
    private:
        bool _ForminOn;
        
        virtual void computeX(const vec3 &X0, const vec &tau, vec &X){
            for (int i=0; i < _nMonomers; i++){
                for (int d =0; d< 3; d++){
                    X[3*i+d] = X0[d]+tau[d]*i*_spacing;
                }
           }
        }
        
        virtual void rotateTangents(vec3 &Omega, vec &Tau){
            vec3 ThisTau;
            for (int d=0; d < 3; d++){
                ThisTau[d] = Tau[d];
            }
            rotateVector(Omega, ThisTau);
            for (int d=0; d < 3; d++){
                Tau[d] = ThisTau[d];
            }
       }   
          
        
};  

class BranchedFiber: public Fiber{
    
    public:
    
    BranchedFiber(vec3 X0, vec3 tau0, int nMonomers, double spacing, double a, double mu, double kbT, 
        int nLinFib, intvec BranchStartIndex,intvec AttachPoints):Fiber(X0,tau0,nMonomers,{0,0},spacing,a,mu,kbT){
        
        _nLinearFib = nLinFib;
        _BranchStartIndex = intvec(BranchStartIndex.size()); // has nLinFib+1 entries
        std::memcpy(_BranchStartIndex.data(),BranchStartIndex.data(),BranchStartIndex.size()*sizeof(int));
        _AttachPoints = intvec(AttachPoints.size());
        std::memcpy(_AttachPoints.data(),AttachPoints.data(),AttachPoints.size()*sizeof(int));
        // Tangent vectors
        _tau = vec(3*nLinFib);
        std::random_device rd; 
        std::mt19937 gen(rd());
        std::normal_distribution<float> NormalDist(0,1.0); 
        vec3 FirstTau;
        for (int d=0; d< 3; d++){
            FirstTau[d] = NormalDist(gen);
        }
        normalize(FirstTau);
        for (int d=0; d< 3; d++){
            _tau[d] = FirstTau[d];
        }
        for (int i=1; i < nLinFib; i++){
            vec3 AttachedTau;
            // Find the fiber it is attached to
            int BaseFib = 0;
            for (int j=0; j < i; j++){
                if (_BranchStartIndex[j+1] > _AttachPoints[i]){
                    BaseFib=j;
                    break;
                }
            }
            for (int d =0; d < 3; d++){
                AttachedTau[d] = _tau[3*BaseFib+d];
            }
            vec3 NewTau = RandomSeventyDegreeRotation(AttachedTau);
            for (int d=0; d< 3; d++){
                _tau[3*i+d] =NewTau[d];
            }
        }
        _tau[0]=-0.195726;
        _tau[1]= -0.585448;
        _tau[2]=0.786728;
        _tau[3]=0.804808;
        _tau[4]= -0.0597155;
        _tau[5]=0.590524;
        _tau[6]=0.825568 ;
        _tau[7]= 0.206785;
        _tau[8]=-0.525051 ;
        std::cout << "Tangent vectors (ALWAYS SAME!): " << std::endl;
        for (int i=0; i < nLinFib; i++){
            for (int d=0; d< 3; d++){    
               std::cout << _tau[3*i+d] << " , ";
            }
            std::cout << std::endl;
        }
    }
    
    bool isBarbedEnd(int MonIndex) override{
        std::cout << "Check this method " << std::endl;
        for (int iFib =0; iFib < _nLinearFib; iFib++){
            if (MonIndex == _BranchStartIndex[iFib+1]-1){
                return true;
            }
        }
        return false;
    } 
    
    private:
       intvec _BranchStartIndex, _AttachPoints; 
       int _nLinearFib;
       
       void rotateTangents(vec3 &Omega, vec &Tau) override{
            for (int iFib=0; iFib < _nLinearFib; iFib++){
                vec3 ThisTau;
                for (int d=0; d < 3; d++){
                    ThisTau[d] = Tau[3*iFib+d];
                }
                rotateVector(Omega, ThisTau);
                for (int d=0; d < 3; d++){
                    Tau[3*iFib+d] = ThisTau[d];
                }
            }
       }
       
       void computeX(const vec3 &X0, const vec &Taus, vec &X) override{
            // Mother filament
            for (int d=0; d<3; d++){
                X[d]=X0[d];
            }
            for (int iFib=0; iFib < _nLinearFib; iFib++){
                vec3 ThisTau, ThisX0;
                for (int d=0; d < 3; d++){
                    ThisTau[d] = Taus[3*iFib+d];
                    ThisX0[d] = X[3*_AttachPoints[iFib]+d];
                }
                int nMonThisBranch = _BranchStartIndex[iFib+1]-_BranchStartIndex[iFib];
                int startIndex=1; 
                if (iFib==0){
                    startIndex=0;
                }                
                for (int i=startIndex; i < startIndex+nMonThisBranch; i++){
                    int xIndex = _BranchStartIndex[iFib]+(i-1);
                    if (iFib==0){
                        xIndex = _BranchStartIndex[iFib]+i;
                    }
                    for (int d =0; d< 3; d++){
                        X[3*xIndex+d] = ThisX0[d]+ThisTau[d]*i*_spacing;
                    }
                }
            }
        }
        
        vec3 RandomSeventyDegreeRotation(const vec3 &tau){
            // Find a random vector an angle of 70 degrees 
            double seventy = 70.0*M_PI/180.0;
            std::random_device rd; 
            std::mt19937 gen(rd());
            std::normal_distribution<float> NormalDist(0,1.0); 
            vec3 v; 
            for (int d=0; d < 3; d++){
                v[d] = NormalDist(gen);
            }
            double vDotTau = dot(v,tau);
            for (int d=0; d < 3; d++){
                v[d]-=vDotTau*tau[d];
            }
            normalize(v);
            vec3 tau2;
            for (int d=0; d < 3; d++){
                tau2[d] = sin(seventy)*v[d] + cos(seventy)*tau[d];
            }
            return tau2;
        }
};
