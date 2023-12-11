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
static const int NMonomersForFormin=4;

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
        _BarbedBindRate = BindingRates[0];
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
        _X.resize(3*_nMonomers);
        computeX(_X0,_tau,_X);
        return _X;
    }
        
    virtual double PointedBindingRate(){    
        if (_nMonomers == 5){
           // std::cout << "MAX 5 MONOMERS " << std::endl;
            return 0;
        }
        return _PointedBindRate;
    }
    
    virtual void PointedBindReaction(){
        _nMonomers++;
        for (int d=0; d < 3; d++){
            _X0[d]-=_spacing*_tau[d];
        }
    }
    
    virtual double BarbedBindingRateNoFormin(){    
        if (_nMonomers == 5 || _ForminOn){
           // std::cout << "MAX 5 MONOMERS " << std::endl;
            return 0;
        }
        return _BarbedBindRate;
    }
    
    virtual double BarbedBindingRateWithFormin(){
        if (_nMonomers == 5 || !_ForminOn){
            return 0;
        }
        return _BarbedBindRate*_ForminEnhancement;
    }
    
    virtual void BarbedBindReaction(double r, bool ForminEnd){
        _nMonomers++;
    }
    
    double PointedUnbindingRate(){
        return _PointedUnbindRate;
    }
    
    virtual double BarbedUnbindingRate(){
        return _BarbedUnbindRate;
    }
    
    virtual void UnbindReaction(double r, bool PointedEnd){
        _nMonomers--;
        if (PointedEnd) {// adjust the start point
            for (int d=0; d < 3; d++){
                _X0[d]+=_spacing*_tau[d];
            }
        }
    }
    
    virtual double ForminBindRate(double ForminBRate){
        if (_ForminOn || _nMonomers < NMonomersForFormin){
            return 0;
        } 
        return ForminBRate;
    }
    
    virtual double ForminUnbindRate(double ForminUBRate){
        // Formin cannot unbind from filaments with less than 4 monomers
        if (!_ForminOn || _nMonomers < NMonomersForFormin){
            return 0;
        } 
        return ForminUBRate;
    }
    
    virtual void BindFormin(double rU){
        _ForminOn = true;
    }
    
    virtual void UnbindFormin(double rU){
        _ForminOn = false;
    }
            
    virtual int NumMonomers(int FibIndex){
        return _nMonomers;
    }
    
    virtual int TotalMonomers(){
        return _nMonomers;
    }
    
    vec3 getPointedEnd(){
        return _X0;
    } 
    
    virtual void addBranch(double rUF, double rUM, vec3 RandomGaussian, int MinMonForAttach){
        throw std::runtime_error("Cannot add branch to non-branched fiber");
   }
       
    vec3 getBarbedEnd(){
        vec3 Be;
        for (int d =0; d < 3; d++){
            Be[d] = _X0[d]+_tau[d]*(_nMonomers-1)*_spacing;
        }  
        return Be;  
    }
    
    virtual bool ForminBound(int FibIndex){
        return _ForminOn;
    }
    
    virtual uint nFibersEligibleForBranching(int MinForBranching){
        if (_nMonomers >= MinForBranching){
            return 1;
        } else {
            return 0;
        }
    }
    
    virtual uint nFibers(){
        return 1;
    }
        
    protected: 
        int _nMonomers;
        double _spacing;
        vec _tau;
        intvec _MonomerIndices;
        double _BarbedBindRate, _BarbedUnbindRate, _PointedBindRate, _PointedUnbindRate,_ForminEnhancement;
        bool _ForminOn;
        
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
        
        virtual void computeX(const vec3 &X0, const vec &tau, vec &X){
            X.resize(3*_nMonomers);
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
    
    BranchedFiber(vec3 X0, vec3 tau0, intvec nMonomersPerFib, intvec Mothers,intvec AttachPoints, 
        double spacing, double a, double mu, double kbT, const vec &BindingRates, 
        int seed):Fiber(X0,tau0,nMonomersPerFib[0],spacing,a,mu,kbT,BindingRates){
        
        _nLinearFib = nMonomersPerFib.size();
        _nMonomersPerFib = intvec(nMonomersPerFib.size());
        std::memcpy(_nMonomersPerFib.data(),nMonomersPerFib.data(),nMonomersPerFib.size()*sizeof(int));
        _Mothers = intvec(Mothers.size());
        std::memcpy(_Mothers.data(),Mothers.data(),Mothers.size()*sizeof(int));
        _AttachPoints = intvec(AttachPoints.size());
        std::memcpy(_AttachPoints.data(),AttachPoints.data(),AttachPoints.size()*sizeof(int));
        _ForminsOn = boolvec(_nLinearFib,false);
        // Tangent vectors
        _tau = vec(3*_nLinearFib);
        rng.seed(seed);
        normaldist = std::normal_distribution<double>(0.0,1.0);
        for (int d=0; d< 3; d++){
            _tau[d] = tau0[d];
        }
        for (int i=1; i < _nLinearFib; i++){
            vec3 AttachedTau;
            // Find the fiber it is attached to
            int BaseFib = Mothers[i];
            vec3 Randn;
            for (int d =0; d < 3; d++){
                AttachedTau[d] = _tau[3*BaseFib+d];
                Randn[d] = normaldist(rng);
            }
            
            vec3 NewTau = RandomSeventyDegreeRotation(AttachedTau,Randn);
            for (int d=0; d< 3; d++){
                _tau[3*i+d] =NewTau[d];
            }
        }
    }
    
    BranchedFiber(Fiber *Mother, double rU, vec3 RandomGaussian):Fiber(*Mother){
        //std::cout << "Initializing from mother with " << _nMonomers << std::endl;
        _nLinearFib = 2;
        _nMonomersPerFib = intvec(2);
        _nMonomersPerFib[0] = _nMonomers;
        _nMonomersPerFib[1] = 1;
        _AttachPoints = intvec(2);
        _AttachPoints[0] = 0;
        int AttachPt = int(rU*_nMonomers);
        _AttachPoints[1] = AttachPt;
        //std::cout << "Attaching at point " << AttachPt << std::endl;
        _Mothers = intvec(2);
        _Mothers[0] = 0;
        _Mothers[1] = 0;
        _ForminsOn = boolvec(2);
        _ForminsOn[0] = _ForminOn;
        _ForminsOn[1] = false;
        //std::cout << "The old tau " << _tau[0] << " , " << _tau[1] << " , " << _tau[2] << std::endl;
        vec3 MotherTau;
        for (int d=0; d < 3; d++){
            MotherTau[d] = _tau[d];
        }
        vec3 NewTau = RandomSeventyDegreeRotation(MotherTau,RandomGaussian);
        _tau.insert(_tau.end(), NewTau.begin(), NewTau.end());
        //std::cout << "The new tau1 " << _tau[0] << " , " << _tau[1] << " , " << _tau[2] << std::endl;
        //std::cout << "The new tau2 " << _tau[3] << " , " << _tau[4] << " , " << _tau[5] << std::endl;
    }
    
    void addBranch(double rUF, double rUM, vec3 RandomGaussian, int MinMonForAttach) override{
        // Choose the fiber and the point
        intvec EligibleFibs;
        for (int iFib=0; iFib < _nLinearFib; iFib++){ 
            if (_nMonomersPerFib[iFib] >= MinMonForAttach){
                EligibleFibs.push_back(iFib);
            }
        }
        int AttachFib = EligibleFibs[floor(rUF*EligibleFibs.size())];
        if (_nMonomersPerFib[AttachFib] < MinMonForAttach){
            std::cout << "Mistake in choosing fiber for branching" << std::endl;
        }
        _Mothers.push_back(AttachFib);
        _AttachPoints.push_back(int(rUM*_nMonomersPerFib[AttachFib]));
        _nMonomersPerFib.push_back(1);
        _ForminsOn.push_back(false);
        _nLinearFib++;
        vec3 MotherTau;
        for (int d=0; d < 3; d++){
            MotherTau[d] = _tau[3*AttachFib+d];
        }
        vec3 NewTau = RandomSeventyDegreeRotation(MotherTau,RandomGaussian);
        _tau.insert(_tau.end(), NewTau.begin(), NewTau.end());
        //std::cout << "The mother tau " << MotherTau[0] << " , " << MotherTau[1] << " , " << MotherTau[2] << std::endl;
        //std::cout << "The new tau2 " << _tau[3*_nLinearFib-3] << " , " << _tau[3*_nLinearFib-2] << " , " << _tau[3*_nLinearFib-1] << std::endl;
    }
    
    uint nFibersEligibleForBranching(int MinForBranching) override{
        int nElig=0;
        for (int iFib=0; iFib < _nLinearFib; iFib++){ 
            if (_nMonomersPerFib[iFib] >= MinForBranching){
                nElig++;
            }
        }
        return nElig;
    }
    
    uint nFibers() override{
        return _nLinearFib;
    }
    
    bool ForminBound(int FibIndex) override{
        return _ForminsOn[FibIndex];
    }
    
    int NumMonomers(int FibIndex) override{
        return _nMonomersPerFib[FibIndex];
    }
    
    int TotalMonomers() override{
        return  std::accumulate(_nMonomersPerFib.begin(), _nMonomersPerFib.end(),0);
    }
    
    double PointedBindingRate() override{    
        if (_nMonomersPerFib[0] == 5){
           // std::cout << "MAX 5 MONOMERS " << std::endl;
            return 0;
        }
        return _PointedBindRate; // Only 1 pointed end (the mother filament)
    }
    
    void PointedBindReaction() override{
        _nMonomersPerFib[0]++; 
        for (int d=0; d < 3; d++){
            _X0[d]-=_spacing*_tau[d];
        }
    }
    
    double BarbedBindingRateNoFormin() override{    
        double TotalRate=0;
        for (int iBranch=0; iBranch < _nLinearFib; iBranch++){
            if (!_ForminsOn[iBranch] && _nMonomersPerFib[iBranch] < 5){
                TotalRate+=_BarbedBindRate;
            } 
        } 
        return TotalRate;
    }
    
    double BarbedBindingRateWithFormin() override{
        double TotalRate=0;
        for (int iBranch=0; iBranch < _nLinearFib; iBranch++){
            if (_ForminsOn[iBranch] && _nMonomersPerFib[iBranch] < 5){
                TotalRate+=_BarbedBindRate*_ForminEnhancement;
            } 
        } 
        return TotalRate;
    }
    
    void BarbedBindReaction(double r, bool ForminEnd){
        intvec EligibleBranches;
        for (int iBranch=0; iBranch < _nLinearFib; iBranch++){
            if (_ForminsOn[iBranch]==ForminEnd && _nMonomersPerFib[iBranch] < 5){
                EligibleBranches.push_back(iBranch);
            } 
        }
        int BranchNum = EligibleBranches[int(r*EligibleBranches.size())];
        _nMonomersPerFib[BranchNum]++;
        //std::cout << "The chosen branched for barbed binding " << BranchNum << " formin " 
        //    << _ForminsOn[BranchNum] << " compare with " << ForminEnd << std::endl;
    }
    
    double BarbedUnbindingRate() override{
        return _BarbedUnbindRate*_nLinearFib;
    }
    
    void UnbindReaction(double r, bool PointedEnd) override{
        // Need to deal with exceptions and all that here
        if (PointedEnd){
            _nMonomersPerFib[0]--; 
            for (int d=0; d < 3; d++){
                _X0[d]+=_spacing*_tau[d];
            }
        } else {
            // Choose the barbed end at random
            int FibNumber = int(r*_nLinearFib);
            _nMonomersPerFib[FibNumber]--;
            // Remove the branch if it's down to 0 monomers (later)
        }
            
    }
  
    double ForminBindRate(double ForminBRate) override{
        double TotalRate=0;
        for (int iBranch=0; iBranch < _nLinearFib; iBranch++){
            if (!_ForminsOn[iBranch] && _nMonomersPerFib[iBranch]>=NMonomersForFormin){
                TotalRate+=ForminBRate;
            } 
        } 
        return TotalRate;
    }
    
    double ForminUnbindRate(double ForminUBRate) override{
        double TotalRate=0;
        for (int iBranch=0; iBranch < _nLinearFib; iBranch++){
            if (_ForminsOn[iBranch] && _nMonomersPerFib[iBranch]>=NMonomersForFormin){
                TotalRate+=ForminUBRate;
            } 
        }
        return TotalRate;
    }
    
    void BindFormin(double rU) override{
        // Make list of eligible branched
        intvec EligibleBranches;
        for (int iBranch=0; iBranch < _nLinearFib; iBranch++){
            if (!_ForminsOn[iBranch]  && _nMonomersPerFib[iBranch]>=NMonomersForFormin){
                EligibleBranches.push_back(iBranch);
            }
        }
        int TheBranch = EligibleBranches[int(rU*EligibleBranches.size())];
        _ForminsOn[TheBranch] = true;
    }
    
    void UnbindFormin(double rU) override{
        // Make list of eligible branched
        intvec EligibleBranches;
        for (int iBranch=0; iBranch < _nLinearFib; iBranch++){
            if (_ForminsOn[iBranch] && _nMonomersPerFib[iBranch]>=NMonomersForFormin){
                EligibleBranches.push_back(iBranch);
            }
        }
        int TheBranch = EligibleBranches[int(rU*EligibleBranches.size())];
        _ForminsOn[TheBranch] = false;
    }
        
   
    
    private:
       intvec _nMonomersPerFib, _Mothers, _AttachPoints; 
       int _nLinearFib;
       boolvec _ForminsOn;
       std::normal_distribution<double> normaldist;
       std::mt19937_64 rng;
       double _SeventyDegrees = 70.0*M_PI/180.0;
       
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
            int TotalMons = TotalMonomers();
            X.resize(3*TotalMons);
            // Mother filament
            for (int d=0; d<3; d++){
                X[d]=X0[d];
            }
            int nSoFar=0;
            for (int iFib=0; iFib < _nLinearFib; iFib++){
                vec3 ThisTau, ThisX0;
                int Mother = _Mothers[iFib];
                int IndexOnMother = _AttachPoints[iFib];
                int StartIndex = IndexOnMother;
                for (int jFib=0; jFib < Mother; jFib++){
                    StartIndex+=_nMonomersPerFib[jFib];
                }
                for (int d=0; d < 3; d++){
                    ThisTau[d] = Taus[3*iFib+d];
                    ThisX0[d] = X[3*StartIndex+d];
                }
                int nMonThisBranch =_nMonomersPerFib[iFib];
                int FirstMon=1; 
                if (iFib==0){
                    FirstMon=0;
                }                
                for (int i=FirstMon; i < FirstMon+nMonThisBranch; i++){
                    for (int d =0; d< 3; d++){
                        X[3*nSoFar+d] = ThisX0[d]+ThisTau[d]*i*_spacing;
                    }
                    nSoFar++;
                }
            }
        }
         
        
        vec3 RandomSeventyDegreeRotation(const vec3 &tau, vec3 &Randn){
            // Find a random vector an angle of 70 degrees 
            double vDotTau = dot(Randn,tau);
            for (int d=0; d < 3; d++){
                Randn[d]-=vDotTau*tau[d];
            }
            normalize(Randn);
            vec3 tau2;
            for (int d=0; d < 3; d++){
                tau2[d] = sin(_SeventyDegrees)*Randn[d] + cos(_SeventyDegrees)*tau[d];
            }
            return tau2;
        }
};
