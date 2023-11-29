% Parameters
Conc = 2; % in uM
LBox = 2;

% Parameters 
kplusDimer = 3.5e-3; % uM^(-1)*s^(-1) 
kminusDimer = 0.041; %s^(-1)
kplusTrimer = 13e-1; % uM^(-1)*s^(-1) 
kminusTrimer = 22; %s^(-1)
kplusBarbed = 11.6; % uM^(-1)*s^(-1) 
kminusBarbed = 1.4; %s^(-1)
kplusPointed = 1.3; %uM^(-1)*s^(-1)
kminusPointed = 0.8; %s^(-1)

% Convert to microscopic assuming well-mixed system
Volume = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; % everything will be in s^(-1)
RxnRates=[kplusDimer*ConversionFactor; kminusDimer; kplusTrimer*ConversionFactor; ...
    kminusTrimer; kplusBarbed*ConversionFactor; kminusBarbed; ...
    kplusPointed*ConversionFactor; kminusPointed];

Nmon = floor(Conc*Volume/uMInvToMicron3);
V = 4^3;
nMax=5;
% Steady state equations
alphafun = @(alp) OneDFcn(alp,RxnRates,Nmon,nMax);
alpha = fsolve(alphafun,1);
Nums = zeros(nMax,1);
PolyOn = RxnRates(5)+RxnRates(7);
PolyOff = RxnRates(6)+RxnRates(8);
Nums(1) = alpha*(PolyOff/PolyOn);
Nums(2) = Nums(1).^2*RxnRates(1)/RxnRates(2);
Nums(3) = Nums(2).*Nums(1)*RxnRates(3)/RxnRates(4);
for iD=4:nMax
    Nums(iD)=Nums(iD-1)*Nums(1)*PolyOn/PolyOff;
end
% Solve the ODEs
RHSFcn = @(t,y) RHS(t,y,RxnRates);
y0 = [Nmon;zeros(nMax-1,1)];
tf=200;
[tvals,yvals] = ode45(RHSFcn,[0 tf],y0);
plot(tvals,yvals/LBox^3,'-.')
hold on
%plot(tf*ones(5,1),Nums/LBox^3,'ko')

% Import the data
load('StructInfo.txt');
load('NumFibs.txt');
nT = length(NumFibs);
NumOfEach = zeros(5,nT);
TotalEntries = NumFibs+3;
StartIndex = [0;cumsum(TotalEntries)];
for iT=1:nT
    StartInd = StartIndex(iT)+1;
    EndInd = StartIndex(iT)+TotalEntries(iT);
    for iS=0:2
        NumOfEach(1+iS,iT)=StructInfo(StartInd+iS);
    end
    rest = StructInfo(StartInd+3:EndInd);
    NumOfEach(4,iT)=sum(rest==4);
    NumOfEach(5,iT)=sum(rest==5);
end
ts=5:5:tf;
set(gca,'ColorOrderIndex',1)
plot(ts,NumOfEach/LBox^3)


%nMons=[nMons;Nmon];
%alphas=[alphas;alpha];

% lb = zeros(nMax,1);
% myfunC = @(c) EqFcnC(c,kplus_p,kminus_p,kplus_d,kminus_d,Nmon);
% opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
% x3 = fmincon(@(x)0,c0,[],[],[],[],lb,[],myfunC,opts);

% for pmon=1:4
% imon=pmon;
% if (imon>1)
%     imon=pmon+1;
% end
% MeanOfAll=MeanOfAll_3;
% StdOfAll=StdOfAll_3;
% ind=2;
% ts=(0:nTs-1)/(nTs-1)*2000;
% subplot(2,2,pmon)
% set(gca,'ColorOrderIndex',ind)
% plot(ts,MeanOfAll(imon,:))
% hold on
% set(gca,'ColorOrderIndex',ind)
% skip = 50;
% start = 25;
% errorbar(ts(start:skip:end),MeanOfAll(imon,start:skip:end),...
%     2*StdOfAll(imon,start:skip:end)/sqrt(nTrial),'o','MarkerSize',0.1,'LineWidth',2.0)
% plot(xlim,[x(imon) x(imon)],':k')
% title(strcat(num2str(imon),' monomers'))
% end
% % Confidence interval
% mean(MeanOfAll(:,(end-1)/2+1:end)')
% StdOfAll(:,(end-1)/2+1:end)'

function dydt = RHS(t,y,RxnRates)
    nMax = length(y);
    dydt=zeros(length(y),1);
    PolyOn = RxnRates(5)+RxnRates(7);
    PolyOff = RxnRates(6)+RxnRates(8);
    dydt(1) = -2*RxnRates(1)*y(1).^2 + 2*RxnRates(2)*y(2) ...
        -RxnRates(3)*y(1)*y(2) + RxnRates(4)*y(3);
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
    for iD = 4:nMax
        dydt(1) = dydt(1) - PolyOn*y(iD-1)*y(1) + PolyOff*y(iD);
    end
end

function val = OneDFcn(alpha,RxnRates,Nmon,nMax)
    PolyOn = RxnRates(5)+RxnRates(7);
    PolyOff = RxnRates(6)+RxnRates(8);
    Nums = zeros(nMax,1);
    Nums(1) = alpha*(PolyOff/PolyOn);
    Nums(2) = Nums(1).^2*RxnRates(1)/RxnRates(2);
    Nums(3) = Nums(2).*Nums(1)*RxnRates(3)/RxnRates(4);
    for iD=4:nMax
        Nums(iD)=Nums(iD-1)*Nums(1)*PolyOn/PolyOff;
    end
    val = sum((1:nMax)'.*Nums)-Nmon;
    return
    Numerator = 2*alpha-alpha.^2-(1+Nmon)*alpha.^(Nmon)+Nmon*alpha.^(Nmon+1);
    Num1 = Numerator./(1-alpha).^2;
    if (alpha==1)
        Num1=1/2*(-2 + Nmon + Nmon^2);
    end
    val = -Nmon+kminus_p/kplus_p*alpha ...
        + (kplus_d/kminus_d)*(kminus_p/kplus_p)^2*alpha.*Num1;
end