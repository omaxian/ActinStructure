addpath(genpath('../../Python'))
Nf = 8;
Nt = 8;
nTri = 50;
nError = 5;
nSeed = nTri/nError;
dt = 1e-4;
a = 4e-3;
mu = 0.068;
kbT = 4.1e-3;
for iError=1:nError
for seed=1:nSeed
loadseed = nSeed*(iError-1)+seed;
Locs=load(strcat('X_',num2str(loadseed),'.txt'));
XAll0 = Locs(1:Nf,:);
MInverse = 6*pi*a*mu;
K0 = RigidRTMatrix(XAll0);
N0 = pinv(K0'*MInverse*K0);
[U,V]=eig(N0(1:3,1:3));
alpha1 = kbT*(V(2,2)+V(3,3));
u1 = U(:,1);
slope = trace(N0(4:6,4:6))*2*kbT;
numSteps=length(Locs)/Nt;
Coms = zeros(numSteps,3);
tanvecs = zeros(numSteps,3);
for i=1:numSteps
    Coms(i,:)=mean(Locs((i-1)*Nt+1:(i-1)*Nt+Nf,:));
    tanvecs(i,:)=(Locs((i-1)*Nt+Nf,:)-Locs((i-1)*Nt+1,:))/norm(Locs((i-1)*Nt+Nf,:)-Locs((i-1)*Nt+1,:));
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