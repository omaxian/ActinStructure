% Parameters
Nmon=1000;
a = 4e-3;
Raa = 10*a;
lambdaMon = 10;
lambdaOff = 12;
mu = 0.068;
kbT = 4.1e-3;
R=lambdaMon*Raa^2;
D=2*kbT/(6*pi*mu*a);
r=R/D % << 1 is reaction limited, >> 1 is diffusion limited
% Determine effective on rate
kplus = 1/2*4*pi/3*Raa^3*lambdaMon;
kplus_D = 2*pi*D*Raa*(1-sqrt(D/(lambdaMon*Raa^2))*tanh(sqrt(lambdaMon*Raa^2/D)));
kminus = lambdaOff;
% Steady state quadratic equation for concentration
for A=kplus_D
B = kminus/2;
V = 1;
C = -Nmon/V*kminus/2;
cFree =((-B+sqrt(B^2-4*A*C))/(2*A)); % average number of 2's
NumFree = cFree*V;
NumPoly = (Nmon-NumFree)/2;
end

% Simulation data
logdt=2;
dt=10^(-logdt);
tf = 1;
nTs = floor(tf/dt+1e-6)+1;
nTrial=5;
MeanNum = zeros(nTrial,nTs);
nSeed=10;
for iTrial=1:nTrial
    TrialNum = zeros(nSeed,nTs);
    for iseed=1:nSeed
        thisSeed = (iTrial-1)*nSeed+iseed;
        Fibs = load(strcat('Dt',num2str(logdt),'Fibers_',num2str(thisSeed),'.txt'));
        TrialNum(iseed,:)=Fibs;
    end
    MeanNum(iTrial,:)=mean(TrialNum);
end
% Plot the mean and standard deviation and compare to theory
ts=(0:nTs-1)*dt*lambdaMon;
set(gca,'ColorOrderIndex',logdt-1)
plot(ts,mean(MeanNum))
hold on
set(gca,'ColorOrderIndex',logdt-1)
skip = (nTs-1)/100;
start = 1;
errorbar(ts(start:skip:end),mean(MeanNum(:,start:skip:end)),2*std(MeanNum(:,start:skip:end))/sqrt(nTrial),...
    'o','MarkerSize',0.1,'LineWidth',1.0)
% Confidence interval
mean(mean(MeanNum(:,(end-1)/2+1:end)'))
2*std(mean(MeanNum(:,(end-1)/2+1:end)'))/sqrt(nTrial)