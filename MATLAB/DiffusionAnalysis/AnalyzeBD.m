addpath(genpath('../../Python'))
Nt = 1000;
ThisFibInds = 971:980;
nTri = 20;
nError = 2;
nSeed = nTri/nError;
dt = 1e-4;
a = 4e-3;
mu = 0.01;
kbT = 4.1e-3;
for iError=1:nError
for seed=1:nSeed
loadseed = nSeed*(iError-1)+seed;
Locs=load(strcat('X_',num2str(loadseed),'.txt'));
XAll0 = Locs(ThisFibInds,:);
MInverse = 6*pi*a*mu;
K0 = RigidRTMatrix(XAll0);
N0 = pinv(K0'*MInverse*K0);
[U,V]=eig(N0(1:3,1:3));
evals = diag(V);
ind=find(evals<1e-10);
alpha1 = kbT*sum(evals);
u1 = U(:,ind);
slope = trace(N0(4:6,4:6))*2*kbT;
numSteps=length(Locs)/Nt;
Coms = zeros(numSteps,3);
tanvecs = zeros(numSteps,3);
for i=1:numSteps
    ThisStep = Locs((i-1)*Nt+1:i*Nt,:);
    FibThisStep = ThisStep(ThisFibInds,:);
    Coms(i,:)=mean(FibThisStep);
    tanvecs(i,:)=(FibThisStep(end,:)-FibThisStep(1,:))/norm(FibThisStep(end,:)-FibThisStep(1,:));
end
for i=1:numSteps
    disps = Coms(1:end-(i-1),:)-Coms(i:end,:);
    normdisp2 = sum(disps.*disps,2);
    meandisp2(i) = mean(normdisp2);
    %uncdisp2(i) = std(normdisp2)*2/sqrt(numSteps-(i-1));
    dottanvecs = dot(tanvecs(1:end-(i-1),:),tanvecs(i:end,:),2);
    meandotu(i) = mean(dottanvecs);
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