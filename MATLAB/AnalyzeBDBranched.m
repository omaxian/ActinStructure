% Theory
a = 4e-3;
kbT = 4.1e-3;
mu=1;
MInverse=6*pi*a*mu;
X0 = [0 0 0];
nFibers = 3;
nPerFiber = [10; 5; 3];
AttachedTo = [-1; 1; 2]; % 1 is the mother
branchpts = [-1; 5; 3];
Taus = [-0.195726 , -0.585448 , 0.786728;
0.804808 , -0.0597155 , 0.590524;
0.825568 , 0.206785 , -0.525051];
X = BranchedFiberPoints(X0,nPerFiber,AttachedTo,branchpts,Taus,a);
K = RigidRTMatrix(X);
Nnow = (K'*MInverse*K)^(-1);
XAll0=X;
[U,V]=eig(Nnow(1:3,1:3));
alpha1 = kbT*(V(2,2)+V(3,3));
slope = trace(Nnow(4:6,4:6))*2*kbT;
u1 = U(:,1);
Coords = Taus' \ u1;
Nf = 10;
Nt = 28;
nTri = 100;
nError = 5;
nSeed = nTri/nError;
dt = 1e-4;
for iError=1:nError
for seed=1:nSeed
loadseed = nSeed*(iError-1)+seed;
Locs=load(strcat('X_',num2str(loadseed),'.txt'));
XAll0 = Locs(Nf+1:Nt,:);
Tau0 = Locs(Nt,:)-Locs(Nf+1,:)/norm(Locs(Nt,:)-Locs(Nf+1,:));
numSteps=length(Locs)/Nt;
Coms = zeros(numSteps,3);
MaterialVecs = zeros(numSteps,3);
for i=1:numSteps
    Coms(i,:)=mean(Locs((i-1)*Nt+Nf+1:i*Nt,:));
    TausNow = [(Locs((i-1)*Nt+Nf+2,:)-Locs((i-1)*Nt+Nf+1,:))/(2*a);...
        (Locs((i-1)*Nt+Nf+12,:)-Locs((i-1)*Nt+Nf+11,:))/(2*a); ...
        (Locs((i-1)*Nt+Nf+17,:)-Locs((i-1)*Nt+Nf+16,:))/(2*a)];
    MaterialVecs(i,:)=TausNow'*Coords;
end
for i=1:numSteps
    disps = Coms(1:end-(i-1),:)-Coms(i:end,:);
    normdisp2 = sum(disps.*disps,2);
    meandisp2(i) = mean(normdisp2);
    %uncdisp2(i) = std(normdisp2)*2/sqrt(numSteps-(i-1));
    meandotu = sum(MaterialVecs.*u1',2)';
    %uncdotu(i) = std(dottanvecs)*2/sqrt(numSteps-(i-1));
end
allmeandisp2(seed,:)=meandisp2;
allmeandotu(seed,:)=meandotu;
end
L2Displacements(iError,:) = mean(allmeandisp2);
Rotations(iError,:) = mean(allmeandotu);
end
ts = (0:numSteps-1)*dt;
subplot(1,2,1)
errorbar(ts,mean(L2Displacements),2*std(L2Displacements)/sqrt(nError),'LineWidth',2.0)
hold on
plot(ts,slope*ts,'-k');
xlim([0 (numSteps-1)*dt])
subplot(1,2,2)
errorbar(ts,mean(Rotations),2*std(Rotations)/sqrt(nError),'LineWidth',2.0)
hold on
plot(ts,exp(-alpha1*ts),'-k')
xlim([0 (numSteps-1)*dt])
