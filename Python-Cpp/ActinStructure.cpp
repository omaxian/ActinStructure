#include <stdio.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "VectorMethods.cpp"
#include "types.h"
#include <random>

/**
ActinStructure.cpp
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
        _DCoeff = 2.0*_kbT/_Mobility;
        
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
    
    vec3 getX(){
        return _X0;
    }
        
};

class Fiber: public ActinStructure{
    
    public:
    
    // Initialize
    Fiber(vec3 X0, int nMonomers, double spacing, double a, double mu, double kbT){
        std::random_device rd; 
        std::mt19937 gen(0); // rd()
        std::normal_distribution<float> NormalDist(0,1.0); 
        _tau = vec(3);
        vec3 ThisTau;
        for (int d=0; d< 3; d++){
            ThisTau[d] = NormalDist(gen);
        }
        normalize(ThisTau);
        for (int d=0; d< 3; d++){
            _tau[d] = ThisTau[d];
        }
        std::memcpy(_X0.data(),X0.data(),X0.size()*sizeof(double)); 
        _spacing = spacing;
        _nMonomers = nMonomers;
        _Mobility = 1.0/(6*M_PI*mu*a);
        _kbT = kbT;
        _X = vec(3*_nMonomers);
    }
    
    int NumberRand(){
        return 6;
    }
    
    void Diffuse(double dt,const vec &W,bool drift,const vec &WDrift){
        vec K(3*_nMonomers*6);
        computeX(_X0,_tau,_X);
        calcKRigid(_X,K);
        vec KTK(6*6), KTKCopy(6*6);
        BlasMatrixProduct(6, 3*_nMonomers,6, 1.0/_Mobility,0.0,K,true,K,KTK); // K'*M^(-1)*K
        std::memcpy(KTKCopy.data(),KTK.data(),KTK.size()*sizeof(double));
        bool normalize = true;
        bool half = true;
        double pinvtol= 1e-10;
        vec OmegaU(6);
        SolveWithPseudoInverse(6, 6,KTK, W, OmegaU, pinvtol,normalize,half,6); // N^(1/2)
        for (int d=0; d< 6; d++){
            OmegaU[d]*=sqrt(2*_kbT/dt); // sqrt(2kbT/dt)*N^(1/2)*W
        }
        if (drift){ // (N(Xtilde)-N(x))*Wtilde
            double deltaRFD = 1e-5;
            vec XTilde(3*_nMonomers), KTilde(3*_nMonomers*6);
            GetXTilde(_X0,_tau,WDrift,deltaRFD,XTilde);
            calcKRigid(XTilde,KTilde);
            vec KTKTilde(6*6);
            BlasMatrixProduct(6, 3*_nMonomers,6, 1.0/_Mobility,0.0,KTilde,true,KTilde,KTKTilde); // K'*M^(-1)*K
            bool half = false;
            vec DeltaOmegaU1(6);
            SolveWithPseudoInverse(6, 6,KTKTilde, WDrift, DeltaOmegaU1, pinvtol,normalize,half,6); // NTilde*WTilde
            vec DeltaOmegaU2(6);
            SolveWithPseudoInverse(6, 6,KTKCopy, WDrift, DeltaOmegaU2, pinvtol,normalize,half,6); // N*Wtilde
            for (int d=0; d< 6; d++){
                OmegaU[d]+=_kbT/deltaRFD*(DeltaOmegaU1[d]-DeltaOmegaU2[d]);
            }
        }
             
        // Update center of mass
        vec3 Omega;
        for (int d=0; d< 3; d++){
            _X0[d]+=dt*OmegaU[3+d];
            Omega[d]=dt*OmegaU[d];
        }
        rotateTangents(Omega,_tau);
    }
        
       
    vec getX(){
        computeX(_X0,_tau,_X);
        return _X;
    }
    
    protected: 
        int _nMonomers;
        double _spacing;
        vec _tau;
        
        void calcKRigid(const vec &X, vec &K){
            /*
            Kinematic matrix for rigid fibers
            */
            // First calculate cross product matrix
            vec CPMatrix(9*_nMonomers);
            for (int iPt = 0; iPt < _nMonomers; iPt++){
                vec3 disp;
                for (int d=0; d < 3; d++){
                    disp[d]=X[3*iPt+d]-X[d];
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
                        K[rowstart+iEntry]=CPMatrix[Xrowstart+iEntry];
                    }
                }
            } 
            for (int iPt=0; iPt < _nMonomers; iPt++){
                for (int iD=0; iD < 3; iD++){
                    int rowstart = (3*iPt+iD)*6;
                    K[rowstart+3+iD]=1;
                }
            }
        }
        
        void GetXTilde(const vec3 &X0,const vec &tau,const vec &WDrift,double deltaRFD, vec &XTilde){
            vec3 X0Tilde, OmegaTilde;
            vec TauTilde(tau.size());
            std::memcpy(TauTilde.data(),tau.data(),tau.size()*sizeof(double));
            for (int d=0; d < 3; d++){
                X0Tilde[d] = X0[d] + deltaRFD*WDrift[3+d];
                OmegaTilde[d] = deltaRFD*WDrift[d]; 
            }
            rotateTangents(OmegaTilde, TauTilde);
            computeX(X0Tilde,TauTilde,XTilde);
        }   
        
    private:
        
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
    
    BranchedFiber(vec3 X0, int nMonomers, double spacing, double a, double mu, double kbT, 
        int nLinFib, intvec BranchStartIndex,intvec AttachPoints):Fiber(X0,nMonomers,spacing,a,mu,kbT){
        
        _nLinearFib = nLinFib;
        _BranchStartIndex = intvec(BranchStartIndex.size());
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
