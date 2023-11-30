nTrial = 10;
nSeed = 20;
nSteps = 2000;
RotStat = zeros(nSeed,nSteps);
AllRotStats = zeros(nTrial,nSteps);
branched=1;
for iTrial=1:nTrial
MeanDrift=zeros(3,1);
for seed=1:nSeed
a = 4e-3;
X0 = [0 0 0];
Tau = [1 0 0];
NMon = 10;
if (branched)
nPerFiber = [10; 5; 3];
AttachedTo = [-1; 1; 2]; % 1 is the mother
branchpts = [-1; 5; 3];
NMon=sum(nPerFiber);
Tau = [-0.195726 , -0.585448 , 0.786728;
0.804808 , -0.0597155 , 0.590524;
0.825568 , 0.206785 , -0.525051];
end
% Diffusion
kbT = 4.1e-3;
dt = 1e-4;
delta = 1e-4;
MInverse = 6*pi*a*eye(3*NMon);
for iT=1:nSteps
    X = BranchedFiberPoints(X0,NMon,[],[],Tau,a);
    if (branched)
        X = BranchedFiberPoints(X0,nPerFiber,AttachedTo,branchpts,Tau,a);
    end
    Xb = X(1,:);
    K = RigidRTMatrix(X);
    Nnow = pinv(K'*MInverse*K);
    if (iT==1) 
        XAll0=X;
        C=Nnow(1:3,1:3);
        [U,V]=eig(C);
        alpha1 = kbT*(V(2,2)+V(3,3));
        u1 = U(:,1);
        Tau0 = Tau(1,:)';
        if (branched)
            Coords = Tau' \ u1;
        end
    end
    vtilde = randn(3,1);
    %vtilde = load('WDrift.txt');
    %vtilde = vtilde(7:end);
    % Update X0 and rotate taus
    TausTilde = rotate(Tau,delta*vtilde);
    XTilde = BranchedFiberPoints(X0,NMon,[],[],TausTilde,a);
    if (branched)
        XTilde = BranchedFiberPoints(X0,nPerFiber,AttachedTo,branchpts,TausTilde,a);
    end
    KTilde = RigidRTMatrix(XTilde);
    NTilde = pinv(KTilde'*MInverse*KTilde);
    Wn =randn(3,1);%load('W1.txt');
    %Wn=Wn(7:12);
    MeanDrift = MeanDrift+ kbT/delta*(NTilde-Nnow)*vtilde;
    v = sqrt(2*kbT/dt)*real(Nnow^(1/2))*Wn + kbT/delta*(NTilde-Nnow)*vtilde;
    Tau = rotate(Tau,dt*v(1:3)');
    X = BranchedFiberPoints([0 0 0],NMon,[],[],Tau,a);
    if (branched)
        X = BranchedFiberPoints(X0,nPerFiber,AttachedTo,branchpts,Tau,a);
    end
    % Subtract so barbed end stays in place
    if (branched)
        MaterialVec=Tau'*Coords;
        RotStat(seed,iT)=dot(MaterialVec,u1);
    else
        RotStat(seed,iT)=dot(Tau(1,:),Tau0);
    end
end
end
AllMeanDrifts(:,iTrial)=MeanDrift;
AllRotStats(iTrial,:) = mean(RotStat);
end
AllMeanDrifts=AllMeanDrifts/(nSteps*nSeed);
% TODO: rotational diffusion
%subplot(1,2,2)
ts=(1:nSteps)*dt;
errorbar(ts,mean(AllRotStats),2*std(AllRotStats)/sqrt(nTrial),'LineWidth',2.0)
hold on
plot(ts,exp(-alpha1*ts),'-k')