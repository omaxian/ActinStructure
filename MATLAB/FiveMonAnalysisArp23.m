% Parameters
Conc = 2; % in uM
LBox = 2;
ArpConc = 0.5;

% Parameters 
kplusDimer = 3.5e-3; % uM^(-1)*s^(-1) 
kminusDimer = 0;%0.041; %s^(-1)
kplusTrimer = 13e-1; % uM^(-1)*s^(-1) 
kminusTrimer = 0;%22; %s^(-1)
kplusBarbed = 11.6; % uM^(-1)*s^(-1) 
kminusBarbed = 0;%1.4; %s^(-1)
kplusPointed = 1.3; %uM^(-1)*s^(-1)
kminusPointed = 0;%0.8; %s^(-1)

kplusARF = 5.2e-3;

% Convert to microscopic assuming well-mixed system
Volume = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; % everything will be in s^(-1)
RxnRates=[kplusDimer*ConversionFactor; kminusDimer; kplusTrimer*ConversionFactor; ...
    kminusTrimer; kplusBarbed*ConversionFactor; kminusBarbed; ...
    kplusPointed*ConversionFactor; kminusPointed; kplusARF*ConversionFactor^2];

Nmon = floor(Conc*Volume/uMInvToMicron3);
ArpMon = floor(ArpConc*Volume/uMInvToMicron3);
nMax=5;
% Solve the ODEs
RHSFcn = @(t,y) RHS(t,y,RxnRates,nMax);
y0 = [Nmon;zeros(nMax-1,1);zeros(nMax,1);ArpMon];
tf=200;
[tvals,yvals] = ode45(RHSFcn,[0 tf],y0);

% Import the data
load('StructInfo.txt');
load('NumFibs.txt');
load('BoundFormins.txt');
load('BranchedOrLinear.txt');
nT = length(NumFibs);
NumOfEach = zeros(11,nT);
TotalEntries = NumFibs+3;
StartIndex = [0;cumsum(TotalEntries)];
FibStarts = [0;cumsum(NumFibs)];
for iT=1:nT
    StartInd = StartIndex(iT)+1;
    EndInd = StartIndex(iT)+TotalEntries(iT);
    for iS=0:2
        NumOfEach(1+iS,iT)=StructInfo(StartInd+iS);
    end
    rest = StructInfo(StartInd+3:EndInd);
    Branched = BranchedOrLinear(FibStarts(iT)+1:FibStarts(iT+1));
    NumOfEach(4,iT)=sum(rest==4 & ~Branched);
    NumOfEach(5,iT)=sum(rest==5 & ~Branched);
    NumOfEach(6,iT)=sum(rest==1 & Branched);
    NumOfEach(7,iT)=sum(rest==2 & Branched);
    NumOfEach(8,iT)=sum(rest==3 & Branched);
    NumOfEach(9,iT)=sum(rest==4 & Branched);
    NumOfEach(10,iT)=sum(rest==5 & Branched);
    NumOfEach(11,iT)=ArpMon-sum(Branched);
end
ts=0.5:0.5:tf;
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
nexttile
set(gca,'ColorOrderIndex',1)
plot(tvals,yvals(:,1:nMax)/LBox^3,'-.')
hold on
%plot(tf*ones(5,1),Nums/LBox^3,'ko')
plot(ts,NumOfEach(1:nMax,:)/LBox^3)
nexttile
set(gca,'ColorOrderIndex',1)
plot(tvals,yvals(:,nMax+1:2*nMax+1)/LBox^3,'-.')
hold on
%plot(tf*ones(5,1),Nums/LBox^3,'ko')
plot(ts,NumOfEach(nMax+1:2*nMax+1,:)/LBox^3)


function dydt = RHS(t,y,RxnRates,nMax)
    dydt=zeros(length(y),1);
    BarbedPlus = RxnRates(5);
    BarbedMinus = RxnRates(6);
    PolyOn = RxnRates(5)+RxnRates(7);
    PolyOff = RxnRates(6)+RxnRates(8);
    ARFPlus = RxnRates(9);
    BranchElig=0;
    for iD=4:nMax
        BranchElig=BranchElig+y(iD)+y(nMax+iD);
    end
    dydt(1) = -2*RxnRates(1)*y(1).^2 + 2*RxnRates(2)*y(2) ...
        -RxnRates(3)*y(1)*y(2) + RxnRates(4)*y(3) ...
        -ARFPlus*y(1)*y(2*nMax+1)*BranchElig;
    for iD = 4:nMax
        dydt(1) = dydt(1) - PolyOn*y(iD-1)*y(1)+ PolyOff*y(iD);
    end
    for iD = 1:nMax-1
        dydt(1) = dydt(1) - BarbedPlus*y(nMax+iD)*y(1);
    end
    % Dimers: form, unform, form trimers, unform trimers
    dydt(2) = RxnRates(1)*y(1).^2 + RxnRates(4)*y(3) ...
        - RxnRates(2)*y(2) - RxnRates(3)*y(1)*y(2);
    % Trimers: form, unform, form tetramers, unform tetramers
    dydt(3) = RxnRates(3)*y(1)*y(2) - RxnRates(4)*y(3) ...
        - PolyOn*y(1)*y(3) + PolyOff*y(4); % trimers
    for iD = 4:nMax-1
        dydt(iD)=PolyOn*y(iD-1)*y(1) - PolyOff*y(iD) - PolyOn*y(iD)*y(1) + PolyOff*y(iD+1);
    end
    iD = nMax;
    dydt(iD)=PolyOn*y(iD-1)*y(1) - PolyOff*y(iD);
    % Arp 2/3 ones
    dydt(nMax+1) = ARFPlus*y(1)*y(2*nMax+1)*BranchElig - BarbedPlus*y(nMax+1)*y(1);
    for iD=2:nMax-1
        dydt(nMax+iD) = BarbedPlus*y(1)*(y(nMax+iD-1)-y(nMax+iD));
    end
    dydt(2*nMax)= BarbedPlus*y(1)*y(2*nMax-1);
    dydt(2*nMax+1)=- ARFPlus*y(1)*y(2*nMax+1)*BranchElig;
end