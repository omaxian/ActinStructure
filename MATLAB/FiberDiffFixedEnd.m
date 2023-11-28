nTrial = 10;
nSeed = 20;
nSteps = 2000;
RotStat = zeros(nSeed,nSteps);
AllRotStats = zeros(nTrial,nSteps);
for iTrial=1:nTrial
MeanDrift=zeros(3,1);
for seed=1:nSeed
a = 4e-3;
X0 = [0 0 0];
Tau = [1 0 0];
NMon = 10;
% Diffusion
kbT = 4.1e-3;
dt = 1e-4;
delta = 1e-4;
MInverse = 6*pi*a*eye(3*NMon);
for iT=1:nSteps
    X = BranchedFiberPoints(X0,NMon,[],[],Tau,a);
    Xc = mean(X);
    Xb = X(end,:);
    K = RigidRTMatrix(X);
    Nnow = pinv(K'*MInverse*K);
    if (iT==1) 
        XAll0=X;
        C=Nnow(1:3,1:3);
        [U,V]=eig(C);
        alpha1 = kbT*(V(2,2)+V(3,3));
        u1 = U(:,1);
        Tau0 = Tau';
    end
    vtilde = randn(3,1);
    %vtilde = load('WDrift.txt');
    %vtilde = vtilde(7:end);
    % Update X0 and rotate taus
    TausTilde = rotate(Tau,delta*vtilde);
    XTilde = BranchedFiberPoints([0 0 0],NMon,[],[],TausTilde,a);
    deltaX = XTilde(end,:)-Xb;
    XTilde = XTilde-deltaX;
    KTilde = RigidRTMatrix(XTilde);
    NTilde = pinv(KTilde'*MInverse*KTilde);
    Wn =randn(3,1);%load('W1.txt');
    %Wn=Wn(7:12);
    MeanDrift = MeanDrift+ kbT/delta*(NTilde-Nnow)*vtilde;
    v = sqrt(2*kbT/dt)*real(Nnow^(1/2))*Wn + kbT/delta*(NTilde-Nnow)*vtilde;
    % Update COM
    %X0 = Xc+dt*v(4:6)'+rotate(X0-Xc,dt*v(1:3)');
    Tau = rotate(Tau,dt*v(1:3)');
    X = BranchedFiberPoints([0 0 0],NMon,[],[],Tau,a);
    % Subtract so barbed end stays in place
    deltaX = X(end,:)-Xb;
    X = X-deltaX;
    X0 = X(1,:);
    RotStat(seed,iT)=dot(Tau,Tau0);
end
end
AllMeanDrifts(:,iTrial)=MeanDrift;
AllRotStats(iTrial,:) = mean(RotStat);
end
% TODO: rotational diffusion
%subplot(1,2,2)
errorbar(ts,mean(AllRotStats),2*std(AllRotStats)/sqrt(nTrial),'LineWidth',2.0)
hold on
plot(ts,exp(-alpha1*ts),'-k')