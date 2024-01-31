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
static const int NMonForFiber = 4;
static const int NMonForBranch = 4;
static const bool MaxFive = false;

// Global variables for periodic ActinStructure
class ActinStructure {
    
    public:
    
    ActinStructure(){

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
    std::normal_distribution<double> normaldist;
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> unifdist;
    std::mt19937_64 rngu;
    
    double logrand(){
        return -log(1.0-unifdist(rngu));
    }
    
    vec3 PointOnUnitSphere(){
        vec3 tau;
        for (int d=0; d < 3; d++){
            tau[d] = normaldist(rng);
        }
        normalize(tau);    
        return tau;
    }
    
    vec3 UniformPointInBox(vec3 Lens){
        vec3 X0;
        for (int d=0; d < 3; d++){
            X0[d]=unifdist(rngu)*Lens[d];
        }
        return X0;
    }
    
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
    
    Fiber(vec3 Lens, uint nMonomers, const int &nMonProts, const int &nBarbedProts, 
          const vec &BarbedBindersRates, const vec &PointedPlus, const vec &BarbedPlus, 
          double PointedMinus, const vec &BarbedMinus, const vec &BranchOnRates, double BranchOffRate,
          double spacing, double a, double mu, double kbT, int seed){

        // Initialize random distributions
        if (seed >=0){
            rngu.seed(seed);
            rng.seed(seed);
        } else {
            auto seed1 = std::random_device{}();
            auto seed2 = std::random_device{}();
            rngu.seed(seed1);
            rng.seed(seed2);
        }
            
        unifdist = std::uniform_real_distribution<double>(0.0,1.0);
        normaldist = std::normal_distribution<double>(0.0,1.0);
        
        // Initialize Tau and X0
        vec3 tau = PointOnUnitSphere();
        _tau = vec(3);
        for (int d=0; d< 3; d++){
            _tau[d] = tau[d];
        }
        _X0 = UniformPointInBox(Lens);
        
        _spacing = spacing;
        _nMonomers = nMonomers;
        _Mobility = 1.0/(6*M_PI*mu*a);
        _kbT = kbT;
        _X = vec(3*_nMonomers);
        
        // These are the rates for unbinding
        _PointedUnbindRate = PointedMinus;

        _ProteinAtBarbedEnd = 0;
        _nMonProts = nMonProts;
        _nBarbedProts = nBarbedProts;
        _PointedBindRate = vec(PointedPlus.size());
        _BarbedBindRate = vec(BarbedPlus.size());
        _BarbedUnbindRate = vec(BarbedMinus.size());
        std::memcpy(_PointedBindRate.data(),PointedPlus.data(),PointedPlus.size()*sizeof(double));
        std::memcpy(_BarbedBindRate.data(),BarbedPlus.data(),BarbedPlus.size()*sizeof(double));
        std::memcpy(_BarbedUnbindRate.data(),BarbedMinus.data(),BarbedMinus.size()*sizeof(double));
        
        _BarbedBindersOnOffRates = vec(BarbedBindersRates.size());
        std::memcpy(_BarbedBindersOnOffRates.data(),BarbedBindersRates.data(),BarbedBindersRates.size()*sizeof(double)); 
        
        _BranchOnRates = vec(BranchOnRates.size());
        std::memcpy(_BranchOnRates.data(),BranchOnRates.data(),BranchOnRates.size()*sizeof(double)); 
        _BranchOffRate = BranchOffRate;
        _NextReactionBranch = -1;
    }
    
    void SetX0AndTau(const vec3 &X0, const vec3 &tau){
        for (int d=0; d< 3; d++){
            _tau[d] = tau[d];
        }
        std::memcpy(_X0.data(),X0.data(),X0.size()*sizeof(double)); 
    }
       
    void Diffuse(double dt){
        vec W(6);
        for (int j=0; j < 6; j++){
            W[j]=normaldist(rng);
        }
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
    
    virtual double NextEventTime(const intvec &nActinMonomers, const intvec &nBarbedBinders, const int &nBranchers){
        // Polymerization on pointed end: Rxns [0,_nMonProts)
        _NextReactionIndex = 0;
        _NextReactionBranch = -1;
        int NumPossMon = _nMonProts+1;
        double deltaT = 1.0/0.0;
        double TryDt = deltaT;
        int BranchIndex=-1;
        for (int iMonProtein = 0; iMonProtein < NumPossMon; iMonProtein++){
            TryDt = PointedBindingTime(nActinMonomers[iMonProtein],iMonProtein);
            if (TryDt < deltaT){
                deltaT = TryDt;
                _NextReactionIndex = iMonProtein;
                _NextReactionBranch = 0;
            }
        }
        // Polymerization on barbed end: Rxns [NumPossMon,2*NumPossMon)
        for (int iMonProtein = 0; iMonProtein < NumPossMon; iMonProtein++){
            TryDt = BarbedBindingTime(nActinMonomers[iMonProtein],iMonProtein,BranchIndex);
            if (TryDt < deltaT){
                deltaT = TryDt;
                _NextReactionIndex = NumPossMon+iMonProtein;
                _NextReactionBranch = BranchIndex;
            }
        }
        // Depolymerization on pointed end: Rxn 2*NumPossMon
        TryDt = PointedUnbindingTime();
        if (TryDt < deltaT){
            deltaT = TryDt;
            _NextReactionIndex = 2*NumPossMon;
            _NextReactionBranch = 0;
        }
        // Depolymerization on barbed end: Rxn 2*NumPossMon+1
        TryDt = BarbedUnbindingTime(BranchIndex);
        if (TryDt < deltaT){
            deltaT = TryDt;
            _NextReactionIndex = 2*NumPossMon+1;
            _NextReactionBranch = BranchIndex;
        }
        // Reactions with things bound/unbound from barbed end
        if (_nBarbedProts > 0){
            TryDt = BarbedBoundProteinUnbindTime(BranchIndex);
            if (TryDt < deltaT){
                deltaT = TryDt;
                _NextReactionIndex = 2*NumPossMon+2;
                _NextReactionBranch = BranchIndex;
            }
            for (int iBarbedProt = 0; iBarbedProt < _nBarbedProts; iBarbedProt++){
                TryDt = BarbedProteinBindTime(iBarbedProt, nBarbedBinders[iBarbedProt], BranchIndex);
                if (TryDt < deltaT){
                    deltaT = TryDt;
                    _NextReactionIndex = 2*NumPossMon+3+iBarbedProt;
                    _NextReactionBranch = BranchIndex;
                }
            }
        }
        // Forming/removing branches
        if (_BranchOnRates.size() > 0){
            for (int iMonProtein = 0; iMonProtein < NumPossMon; iMonProtein++){
                TryDt = BranchFormingTime(nBranchers,nActinMonomers[iMonProtein],iMonProtein,BranchIndex);
                if (TryDt < deltaT){
                    deltaT = TryDt;
                    _NextReactionIndex = 2*NumPossMon+3+_nBarbedProts;
                    _NextReactionBranch = BranchIndex;
                }
            }
            TryDt = BranchRemovalTime(BranchIndex);
            if (TryDt < deltaT){
                deltaT = TryDt;
                _NextReactionIndex = 2*NumPossMon+3+_nBarbedProts+1;
                _NextReactionBranch = BranchIndex;
            }
        }
        return deltaT;
    }
       
    virtual int NeedsStructureChange(){
        int NumPossMon = _nMonProts+1;
        if (_NextReactionIndex == 2*NumPossMon+3+_nBarbedProts){ // branch formation
            return 1;
        }
        if (_NextReactionIndex ==2*NumPossMon || _NextReactionIndex ==2*NumPossMon+1){ // remove me
            if (_nMonomers == NMonForFiber){
                return -1;
            }
        }
        return 0;
    }
    
    virtual void ReactNextEvent(int &nFreeMonomers,intvec &nBarbedBinders, int &FreeBranchers){
        int NumPossMon = _nMonProts+1;
        if (_NextReactionIndex < NumPossMon){ // Polymerization on pointed end
            nFreeMonomers--;
            PointedBindReaction();
        } else if (_NextReactionIndex < 2*NumPossMon){
            nFreeMonomers--;
            BarbedBindReaction();
        } else if (_NextReactionIndex < 2*NumPossMon+2){
            nFreeMonomers++;
            bool Pointed = (_NextReactionIndex==2*NumPossMon);
            UnbindReaction(Pointed);
        } else if (_NextReactionIndex == 2*NumPossMon+2){
            BoundBarbedProteinUnbindReaction(nBarbedBinders);
        } else if (_NextReactionIndex < 2*NumPossMon+3+_nBarbedProts){ 
            int ProtToBind = _NextReactionIndex - (2*NumPossMon+2);
            BindBarbedProteinReaction(ProtToBind, nBarbedBinders);
        } else if (_NextReactionIndex == 2*NumPossMon+3+_nBarbedProts){
            FreeBranchers--;
            nFreeMonomers--;
            addBranch(_NextReactionBranch);
        } else {
            FreeBranchers++;
            nFreeMonomers++;
            removeBranch(_NextReactionBranch, nBarbedBinders);
        }
    }
           
    vec getX(){
        _X.resize(3*_nMonomers);
        computeX(_X0,_tau,_X);
        return _X;
    }
        
            
    virtual int NumMonomers(int FibIndex){
        return _nMonomers;
    }
    
    virtual void SetAttachPt(const int &iBranch, const int &AttachPt){
        //std::cout << "Set attach pts is purely abstract for linear fibers " << std::endl;
    }
    
    virtual void SetNumMons(const int &iBranch, const int &nMons){
        _nMonomers = nMons;
    }
    
    virtual int TotalMonomers(){
        return _nMonomers;
    }
    
    virtual vec getX0(){
        return vec(std::begin(_X0), std::end(_X0));;
    } 
    
    vec getTau(){
        return _tau;
    }
    
    virtual void SetBoundBarbed(int jFib, int iB){
        _ProteinAtBarbedEnd=iB;
    }
    
    virtual int getBoundBarbed(int jFib){
        return _ProteinAtBarbedEnd;
    }
       
    vec3 getBarbedEnd(){
        vec3 Be;
        for (int d =0; d < 3; d++){
            Be[d] = _X0[d]+_tau[d]*(_nMonomers-1)*_spacing;
        }  
        return Be;  
    }
    
    virtual int PrepareForLinearSwitch(vec3 &X0, vec3 &Tau){
        return 0;
    }
    
    virtual uint nFibers(){
        return 1;
    }
    
    virtual void addBranch(int Mother){
        std::cout << "This method (add branch) just virtual - should never be called " << std::endl;
    }
        
    protected: 
        int _nMonomers;
        int _ProteinAtBarbedEnd, _nBarbedProts, _nMonProts;
        int _NextReactionIndex, _NextReactionBranch;
        double _spacing;
        vec _tau, _BarbedBindersOnOffRates, _BranchOnRates;
        double _PointedUnbindRate, _BranchOffRate;
        vec _BarbedBindRate, _PointedBindRate, _BarbedUnbindRate;
        
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
        
        int MatrixIndex(int iMB, int iBB){
            return (_nBarbedProts+1)*iMB+iBB;     
        }
        
        virtual double PointedBindingTime(int nMons, int MonBoundProtein){    
            if (MaxFive && _nMonomers == 5){
                return 1.0/0.0;
            }
            double Rate = _PointedBindRate[MonBoundProtein]*nMons;
            return logrand()/Rate;
        }
    
        virtual void PointedBindReaction(){
            _nMonomers++;
            for (int d=0; d < 3; d++){
                _X0[d]-=_spacing*_tau[d];
            }
        }
        
        virtual double BarbedBindingTime(const int &nMons, const int &MonBoundProtein, int &BranchIndex){    
            if (MaxFive && _nMonomers == 5){
                return 1.0/0.0;
            }
            int MatIndex = MatrixIndex(MonBoundProtein, _ProteinAtBarbedEnd);
            double Rate = nMons*_BarbedBindRate[MatIndex];
            return logrand()/Rate;
        }
        
        
        virtual void BarbedBindReaction(){
            _nMonomers++;
        }
        
        virtual double PointedUnbindingTime(){
            return logrand()/ _PointedUnbindRate;
        }
        
        virtual double BarbedUnbindingTime(int &BranchIndex){
            return logrand()/_BarbedUnbindRate[_ProteinAtBarbedEnd];
        }
        
        virtual void UnbindReaction(bool PointedEnd){
            _nMonomers--;
            if (PointedEnd) {// adjust the start point
                for (int d=0; d < 3; d++){
                    _X0[d]+=_spacing*_tau[d];
                }
            }
        }
        
        virtual double BarbedBoundProteinUnbindTime(int &BranchIndex){
            if (_ProteinAtBarbedEnd == 0){
                return 1.0/0.0;
            }    
            // Unbinding reaction
            return logrand()/(_BarbedBindersOnOffRates[2*_ProteinAtBarbedEnd-1]);
        }
        
        virtual void BoundBarbedProteinUnbindReaction(intvec &nBarbedBinders){
            nBarbedBinders[_ProteinAtBarbedEnd-1]++;
            _ProteinAtBarbedEnd = 0;
        }
        
        virtual double BarbedProteinBindTime(const int &iBarbedProt, const int &numBarbedProts, int &BranchIndex){
            if (_ProteinAtBarbedEnd > 0){
                return 1.0/0.0;
            } 
            return logrand()/(_BarbedBindersOnOffRates[2*iBarbedProt]*numBarbedProts);  
        }
        
        virtual void BindBarbedProteinReaction(const int &ProtToBind, intvec &nBarbedBinders){
            _ProteinAtBarbedEnd = ProtToBind;
            nBarbedBinders[ProtToBind-1]--;
        }
        
        virtual double BranchFormingTime(const int &nBranchers, const int &nMons, const int &MonBoundProtein, int &BranchIndex){
            double BranchRate = 0;
            if (_nMonomers >= NMonForBranch){
                BranchRate =  _BranchOnRates[MonBoundProtein]*nBranchers*nMons*_nMonomers;
            }
            return logrand()/BranchRate;
        }
                
        virtual double BranchRemovalTime(int &BranchIndex){
            // Linear fibers don't have branches to remove. Time is infinity
            return 1.0/0.0;
        }
        
        virtual void removeBranch(int BranchNum, intvec &nBarbedBinders){
            std::cout << "This method (remove branch) just virtual - should never be called " << std::endl;
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
    
        BranchedFiber(Fiber *Mother):Fiber(*Mother){
            //std::cout << "Initializing from mother with " << _nMonomers << std::endl;
            _nLinearFib = 2;
            _nMonomersPerFib = intvec(2);
            _nMonomersPerFib[0] = _nMonomers;
            _nMonomersPerFib[1] = 1;
            _AttachPoints = intvec(2);
            _AttachPoints[0] = 0;
            double ru = unifdist(rngu);
            int AttachPt = NMonForBranch+int(ru*(_nMonomersPerFib[0]-NMonForBranch+1))-1;
            _AttachPoints[1] = AttachPt;
            //std::cout << "Random " << rU << "->Attaching at point " << _AttachPoints[1] << " out of " << _nMonomers << std::endl;
            _Mothers = intvec(2);
            _Mothers[0] = 0;
            _Mothers[1] = 0;
            _ProteinAtBarbedEnds = intvec(2);
            _ProteinAtBarbedEnds[0] = _ProteinAtBarbedEnd;
            _ProteinAtBarbedEnds[1]  = 0;
            vec3 MotherTau;
            for (int d=0; d < 3; d++){
                MotherTau[d] = _tau[d];
            }
            vec3 NewTau = RandomSeventyDegreeRotation(MotherTau);
            _tau.insert(_tau.end(), NewTau.begin(), NewTau.end());
        }
            
        int PrepareForLinearSwitch(vec3 &X0, vec3 &Tau) override{
            _tau.resize(3);
            for (int d=0; d< 3; d++){
                Tau[d]=_tau[d];
                X0[d]=_X0[d];
            }
            return _nMonomersPerFib[0];
            //std::cout << "Number monomers " << _nMonomers << std::endl;
            //std::cout << "First entry of X0 and Tau: " << _X0[0] << " , " << _tau[0] << std::endl;
        }
              
        int NeedsStructureChange() override{
            // Has to be branch removal and only 1 branch
            if (_nLinearFib==2 && _NextReactionIndex==2*(_nMonProts+1)+3+_nBarbedProts+1){
                return -1;
            }
            return 0;
        }
           
        void SetAttachPt(const int &iBranch, const int &AttachPt) override{
            _AttachPoints[iBranch]=AttachPt;
        }
        
        void SetNumMons(const int &iBranch, const int &nMons) override{
            _nMonomersPerFib[iBranch] = nMons; 
        }
        
        uint nFibers() override{
            return _nLinearFib;
        }
        
        vec getX0() override{
            // This will return the first point on each branch (NOT the mother point)
            vec AllX0(3*_nLinearFib);
            for (int d=0; d<3; d++){
                AllX0[d]=_X0[d];
            }
            for (int iFib=1; iFib < _nLinearFib; iFib++){
                int Mother = _Mothers[iFib];
                int IndexOnMother = _AttachPoints[iFib];               
                for (int d =0; d< 3; d++){
                    AllX0[3*iFib+d] = AllX0[3*Mother+d]+_tau[3*Mother+d]*_spacing*IndexOnMother;
                    AllX0[3*iFib+d]+=_tau[3*iFib+d]*_spacing; // add the current fiber spacing
                }
            }
            return AllX0;
        }
        
        void SetBoundBarbed(int jFib, int iB) override{
            _ProteinAtBarbedEnds[jFib] = iB;
        }
    
        int getBoundBarbed(int jFib) override{
            return _ProteinAtBarbedEnds[jFib];
        }
        
        int NumMonomers(int FibIndex) override{
            return _nMonomersPerFib[FibIndex];
        }
        
        int TotalMonomers() override{
            return  std::accumulate(_nMonomersPerFib.begin(), _nMonomersPerFib.end(),0);
        }
        
        void addBranch(int Mother) override{
            // Choose the fiber and the point
            _Mothers.push_back(Mother);
            double rUM = unifdist(rngu);
            int AttachPt = NMonForBranch+int(rUM*(_nMonomersPerFib[Mother]-NMonForBranch+1))-1;
            //std::cout << "Random " << rUM << "->Attaching at point " << AttachPt << " out of " << _nMonomersPerFib[AttachFib] << std::endl;
            _AttachPoints.push_back(AttachPt);
            _nMonomersPerFib.push_back(1);
            _ProteinAtBarbedEnds.push_back(0);
            _nLinearFib++;
            vec3 MotherTau;
            for (int d=0; d < 3; d++){
                MotherTau[d] = _tau[3*Mother+d];
            }
            vec3 NewTau = RandomSeventyDegreeRotation(MotherTau);
            _tau.insert(_tau.end(), NewTau.begin(), NewTau.end());
            //std::cout << "The mother tau " << MotherTau[0] << " , " << MotherTau[1] << " , " << MotherTau[2] << std::endl;
            //std::cout << "The new tau2 " << _tau[3*_nLinearFib-3] << " , " << _tau[3*_nLinearFib-2] << " , " << _tau[3*_nLinearFib-1] << std::endl;
        }
           
    private:
        intvec _nMonomersPerFib, _Mothers, _AttachPoints, _ProteinAtBarbedEnds;
        int _nLinearFib;
        double _SeventyDegrees = 70.0*M_PI/180.0;
        
        double PointedBindingTime(int nMons, int MonBoundProtein) override{    
            if (MaxFive && _nMonomersPerFib[0] == 5){
                return 1.0/0.0;
            }
            double Rate = _PointedBindRate[MonBoundProtein]*nMons;
            return logrand()/Rate;
        }      
    
        void PointedBindReaction() override{
            _nMonomersPerFib[0]++; 
            for (int d=0; d < 3; d++){
                _X0[d]-=_spacing*_tau[d];
            }
            // Change the indicies of the branches bound to the mother
            for (int jBranch=1; jBranch < _nLinearFib; jBranch++){
                if (_Mothers[jBranch]==0){
                    _AttachPoints[jBranch]++;
                }
            }
        }
        
        double BarbedBindingTime(const int &nMons, const int &MonBoundProtein, int &BranchIndex) override{ 
            double MinTime = 1.0/0.0;
            for (int iBranch=0; iBranch < _nLinearFib; iBranch++){   
                bool Eligible = !MaxFive || (_nMonomersPerFib[iBranch] < 5); 
                int MatIndex = MatrixIndex(MonBoundProtein, _ProteinAtBarbedEnds[iBranch]);
                double Rate = nMons*_BarbedBindRate[MatIndex];
                double ThisTime = logrand()/Rate;
                if (ThisTime < MinTime && Eligible){
                    MinTime = ThisTime;
                    BranchIndex = iBranch;
                }
            }
            return MinTime;
        }
                
        void BarbedBindReaction() override{
            _nMonomersPerFib[_NextReactionBranch]++;
            if (_nMonomersPerFib[_NextReactionBranch] > 5){
                std::cout << "In barbed bind too many " << std::endl;
                throw std::runtime_error("Bad");
            }
        }
        
        double PointedUnbindingTime() override{
            // Return 0 if there is a branch sitting on the fourth monomer. Otherwise allow unbind
            for (int jBranch=1; jBranch < _nLinearFib; jBranch++){
                if (_Mothers[jBranch]==0 &&_AttachPoints[jBranch]==NMonForBranch-1){
                    return 1.0/0.0;
                } 
            }
            return  logrand()/_PointedUnbindRate;
        }
                
        double BarbedUnbindingTime(int &BranchIndex) override{
            // We only allow branches of length 2 or more to depolymerize from the barbed end
            // Unbinding is only allowed if there is no branch sitting on that monomer
            double MinTime = 1.0/0.0;
            for (int iBranch=0; iBranch < _nLinearFib; iBranch++){
                bool BranchAtBarbedEnd = false;
                for (int jBranch=1; jBranch < _nLinearFib; jBranch++){
                    if (_Mothers[jBranch]==iBranch &&
                        _AttachPoints[jBranch]==_nMonomersPerFib[iBranch]-1){
                        //std::cout << "Disallowing unbinding on fiber " << iBranch << " because there is branch " << jBranch 
                        //<< " attached there  at point " << _AttachPoints[jBranch] << " and iBranch has " 
                        //<< _nMonomersPerFib[iBranch] << " monomers." << std::endl;
                        BranchAtBarbedEnd=true;
                    }
                }
                //std::cout << "Number of monomers on branch " << iBranch << " = " << _nMonomersPerFib[iBranch] << std::endl;
                if (_nMonomersPerFib[iBranch] > 1 && !BranchAtBarbedEnd){
                    double ThisTime = logrand()/_BarbedUnbindRate[_ProteinAtBarbedEnds[iBranch]];
                    if (ThisTime < MinTime){
                        MinTime = ThisTime;
                        BranchIndex = iBranch;
                    }         
                } 
            }
            return MinTime;
        }
        
        void UnbindReaction(bool PointedEnd) override{
            if (PointedEnd){
                //std::cout << "Pointed end unbind" << std::endl;
                _nMonomersPerFib[0]--; 
                for (int d=0; d < 3; d++){
                    _X0[d]+=_spacing*_tau[d];
                }
                // Change the indicies of the branches bound to the mother
                for (int jBranch=1; jBranch < _nLinearFib; jBranch++){
                    if (_Mothers[jBranch]==0){
                        _AttachPoints[jBranch]--;
                    }
                }
            } else {
                //std::cout << "Unbinding at branch " << BranchNum << std::endl;
                _nMonomersPerFib[_NextReactionBranch]--;
            }
            if (_nMonomersPerFib[0] <NMonForBranch){
                std::cout << "Number monomers on mother has dropped below 4" << std::endl;
                std::cout << "Number of fibers " << _nLinearFib << std::endl;
                std::cout << "Attachment point on mother " << _AttachPoints[1] << std::endl;
                std::cout << "Type of unbind " << PointedEnd << std::endl;
                throw std::runtime_error("STOP");
            }
        }
        
        double BarbedBoundProteinUnbindTime(int &BranchIndex) override{
            double MinTime = 1.0/0.0;
            for (int iFib = 0; iFib < _nLinearFib; iFib++){
                int ProteinBound = _ProteinAtBarbedEnds[iFib];
                if (ProteinBound > 0){
                    double ThisTime = logrand()/_BarbedBindersOnOffRates[2*ProteinBound-1];
                    if (ThisTime < MinTime){
                        MinTime = ThisTime;
                        BranchIndex = iFib;
                    }    
                }
            }    
            return MinTime;
        }
        
        void BoundBarbedProteinUnbindReaction(intvec &nBarbedBinders) override{
            nBarbedBinders[_ProteinAtBarbedEnds[_NextReactionBranch]-1]++;
            _ProteinAtBarbedEnds[_NextReactionBranch] = 0;
        }
        
        double BarbedProteinBindTime(const int &iBarbedProt, const int &numBarbedProts, int &BranchIndex) override{
            double MinTime = 1.0/0.0;
            for (int iFib = 0; iFib < _nLinearFib; iFib++){
                if (_ProteinAtBarbedEnds[iFib]== 0){
                    double ThisTime = logrand()/(_BarbedBindersOnOffRates[2*iBarbedProt]*numBarbedProts);  
                    if (ThisTime < MinTime){
                        MinTime = ThisTime;
                        BranchIndex = iFib;
                    }
                } 
            }
            return MinTime;
        }
        
        void BindBarbedProteinReaction(const int &ProtToBind, intvec &nBarbedBinders) override{
            _ProteinAtBarbedEnds[_NextReactionBranch] = ProtToBind;
            nBarbedBinders[ProtToBind-1]--;
        }
        
        double BranchFormingTime(const int &nBranchers, const int &nMons, const int &MonBoundProtein, int &BranchIndex) override{
            double MinTime = 1.0/0.0;
            double BranchRate =  _BranchOnRates[MonBoundProtein]*nBranchers*nMons;
            for (int iFib=0; iFib < _nLinearFib; iFib++){
                if (_nMonomersPerFib[iFib] >= NMonForBranch){
                    double ThisTime = logrand()/(_nMonomersPerFib[iFib]*BranchRate);
                    if (ThisTime < MinTime){
                        MinTime = ThisTime;
                        BranchIndex = iFib; 
                    }   
                }
            }
            return MinTime;
        }
                
        double BranchRemovalTime(int &BranchIndex) override{
            double MinTime = 1.0/0.0;
            for (int iFib=1; iFib < _nLinearFib; iFib++){
                if (_nMonomersPerFib[iFib]==1){
                    double ThisTime = logrand()/_BranchOffRate;
                    if (ThisTime < MinTime){
                        MinTime = ThisTime;
                        BranchIndex = iFib; 
                    }   
                }
            }
            return MinTime;
        }
        
        void removeBranch(int BranchNum, intvec &nBarbedBinders) override{
            _Mothers.erase(_Mothers.begin() + BranchNum); 
            _nMonomersPerFib.erase(_nMonomersPerFib.begin() + BranchNum); 
            if (_ProteinAtBarbedEnds[BranchNum] > 0){
                nBarbedBinders[_ProteinAtBarbedEnds[BranchNum]-1]++;
            }
            _ProteinAtBarbedEnds.erase(_ProteinAtBarbedEnds.begin() + BranchNum); 
            _AttachPoints.erase(_AttachPoints.begin() + BranchNum); 
            _nLinearFib--;
            //std::cout << "Removed branch " << BranchNum << std::endl;
            for (int d=2; d >=0; d--){
                //std::cout << "Erasing element " << 3*BranchNum+d << " from tau " << std::endl;
                _tau.erase(_tau.begin()+3*BranchNum+d);
            }
            //std::cout << "There are now " << _nLinearFib  << " fibers and tau has size " << _tau.size() << std::endl;
        }
        
        
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
         
        
        vec3 RandomSeventyDegreeRotation(const vec3 &tau){
            // Find a random vector an angle of 70 degrees 
            vec3 Randn;
            for (int d=0; d<3;d++){
                Randn[d]=normaldist(rng);
            }
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
