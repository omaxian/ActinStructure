Nmon=1000;
nSeed=10;
for seed=1:nSeed
X=load(strcat('X_',num2str(seed),'.txt'));
[nSteps,~]=size(X);
nSteps = nSteps/Nmon;
MeanDisps = zeros(1,nSteps);
X0 = X(1:Nmon,:);
for iT=1:nSteps
    disp = X((iT-1)*Nmon+1:iT*Nmon,:) - X0;
    MeanDisps(iT)=mean(sum(disp.*disp,2));
end
AllMeanDisps(seed,:)=MeanDisps;
end
dt = 1e-2;
a = 4e-3;
kbT = 4.1e-3;
mu=15;
slope = kbT/(pi*mu*a);
ts=(0:nSteps-1)*dt;
errorbar(ts,mean(AllMeanDisps),2*std(AllMeanDisps)/sqrt(nSeed));
hold on
%plot(ts,slope*ts,'-k')