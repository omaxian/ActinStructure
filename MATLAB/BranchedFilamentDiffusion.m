nTrial = 10;
nSeed = 10;
nSteps = 1000;
SquareDisp = zeros(nSeed,nSteps);
RotStat = zeros(nSeed,nSteps);
AllSquareDisp = zeros(nTrial,nSteps);
AllRotStats = zeros(nTrial,nSteps);
for iTrial=1:nTrial
for seed=1:nSeed
a = 4e-3;
X0 = [0 0 0];
nFibers = 3;
nPerFiber = [10; 5; 3];
AttachedTo = [-1; 1; 2]; % 1 is the mother
branchpts = [-1; 5; 3];
% Taus = [1 0 0];
% rng(seed);
% for iBranch=2:nFibers
%     Taus(iBranch,:)=  RandomSeventyDegree(Taus(iBranch-1,:));
% end
Taus = [-0.195726 , -0.585448 , 0.786728;
0.804808 , -0.0597155 , 0.590524;
0.825568 , 0.206785 , -0.525051];
% Diffusion
kbT = 4.1e-3;
dt = 1e-4;
delta = 1e-4;
MInverse = 6*pi*a*eye(3*sum(nPerFiber));
for iT=1:nSteps
    X = BranchedFiberPoints(X0,nPerFiber,AttachedTo,branchpts,Taus,a);
    Xc = mean(X);
    K = RigidRTMatrix(X);
    Nnow = (K'*MInverse*K)^(-1);
    % Solve for the mobility center
    if (iT==1) 
        XAll0=X;
        Xmean0=mean(X);
        C=Nnow(1:3,1:3);
        B = Nnow(4:6,1:3);
        slope = trace(Nnow(4:6,4:6))*2*kbT;
        LMatrix = [-C(1,3) -C(2,3) C(1,1)+C(2,2);...
                   C(1,2) -C(3,3)-C(1,1) C(3,2);...
                   C(2,2)+C(3,3) -C(2,1) -C(3,1)];
        RSide = [B(1,2)-B(2,1); B(1,3)-B(3,1); B(2,3)-B(3,2)];
        deltaR = LMatrix \ RSide;
        [U,V]=eig(C);
        alpha1 = kbT*(V(2,2)+V(3,3));
        u1 = U(:,1);
        Coords = Taus' \ u1;
    end
    vtilde = randn(6,1);
    %vtilde = load('WDrift.txt');
    %vtilde = vtilde(7:end);
    % Update X0 and rotate taus
    TausTilde = rotate(Taus,delta*vtilde(1:3)');
    X0Tilde = X0 + delta*vtilde(4:6)';
    XTilde = BranchedFiberPoints(X0Tilde,nPerFiber,AttachedTo,branchpts,TausTilde,a);
    KTilde = RigidRTMatrix(XTilde);
    NTilde = (KTilde'*MInverse*KTilde)^(-1);
    Wn =randn(6,1);%load('W1.txt');
    %Wn=Wn(7:12);
    v = sqrt(2*kbT/dt)*real(Nnow^(1/2))*Wn;% + kbT/delta*(NTilde-Nnow)*vtilde;
    % Update COM
    X0 = Xc+dt*v(4:6)'+rotate(X0-Xc,dt*v(1:3)');
    Taus = rotate(Taus,dt*v(1:3)');
    X = BranchedFiberPoints(X0,nPerFiber,AttachedTo,branchpts,Taus,a);
    SquareDisp(seed,iT)=norm(mean(X)-Xmean0)^2;
    MaterialVec=Taus'*Coords;
    RotStat(seed,iT)=dot(MaterialVec,U(:,1)');
end
end
AllSquareDisp(iTrial,:) = mean(SquareDisp);
AllRotStats(iTrial,:) = mean(RotStat);
end
ts=(1:nSteps)*dt;
subplot(1,2,1)
errorbar(ts,mean(AllSquareDisp),2*std(AllSquareDisp)/sqrt(nTrial))
hold on
plot(ts,ts*slope,'-k')
% TODO: rotational diffusion
subplot(1,2,2)
errorbar(ts,mean(AllRotStats),2*std(AllRotStats)/sqrt(nTrial),'LineWidth',2.0)
hold on
plot(ts,exp(-alpha1*ts),'-k')