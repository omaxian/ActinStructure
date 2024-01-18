% This validates dynamics of formin without branching
% Parameters
Conc = 2; % in uM
LBox = 3;
ForminConc = 0.1;

% Parameters 
kplusDimer = 3.5e-3; % uM^(-1)*s^(-1) 
kminusDimer = 0.041; %s^(-1)
kplusTrimer = 13e-2; % uM^(-1)*s^(-1) 
kminusTrimer = 22; %s^(-1)
kplusBarbed = 1.6; % uM^(-1)*s^(-1) 
kminusBarbed = 1.4; %s^(-1)
kplusPointed = 1.3; %uM^(-1)*s^(-1)
kminusPointed = 0.8; %s^(-1)
kForNuc = 2e-3; % uM^(-2)*s^(-1)
kplusFor = 29.1; %uM^(-1)*s^(-1)
kminusFor = 8.1e-2; %s^(-1)
ForEnhance = 2;

% Convert to microscopic assuming well-mixed system
Volume = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; % everything will be in s^(-1)
RxnRates=[kplusDimer*ConversionFactor; kminusDimer; kplusTrimer*ConversionFactor; ...
    kminusTrimer; kplusBarbed*ConversionFactor; kminusBarbed; ...
    kplusPointed*ConversionFactor; kminusPointed; ...
    kForNuc*ConversionFactor^2; kplusFor*ConversionFactor; kminusFor; ForEnhance];

Nmon = floor(Conc*Volume/uMInvToMicron3);
Nfor = floor(ForminConc*Volume/uMInvToMicron3);
nMax=5;
% Solve the ODEs
RHSFcn = @(t,y) RHS(t,y,RxnRates,nMax);
y0 = [Nmon;zeros(nMax-1,1);Nfor;zeros(nMax-1,1)];
tf=40;
[tvals,yvals] = ode45(RHSFcn,[0 tf],y0);

% Import the data
nError=1;
nTrial=1;
NumFibs=load(strcat('ForminOnly_NumFibs1.txt'));
nT = length(NumFibs);
MeanNumOfEach = zeros(2*nMax,nT,nError);
for iError=1:nError
NumOfEach = zeros(2*nMax,nT,nTrial);
for iTrial=1:nTrial
TriIndex = (iError-1)*nTrial+iTrial;
StructInfo=load(strcat('ForminOnly_StructInfo',num2str(TriIndex),'.txt'));
NumFibs=load(strcat('ForminOnly_NumFibs',num2str(TriIndex),'.txt'));
BoundFormins=load(strcat('ForminOnly_BoundFormins',num2str(TriIndex),'.txt'));
BranchedOrLinear=load(strcat('ForminOnly_BranchedOrLinear',num2str(TriIndex),'.txt'));
nT = length(NumFibs);
TotalEntries = NumFibs+3;
StartIndex = [0;cumsum(TotalEntries)];
FibStarts = [0;cumsum(NumFibs)];
for iT=1:nT
    StartInd = StartIndex(iT)+1;
    EndInd = StartIndex(iT)+TotalEntries(iT);
    for iS=0:2
        NumOfEach(1+iS,iT,iTrial)=StructInfo(StartInd+iS);
    end
    rest = StructInfo(StartInd+3:EndInd);
    Branched = BranchedOrLinear(FibStarts(iT)+1:FibStarts(iT+1));
    HasFormin = BoundFormins(FibStarts(iT)+1:FibStarts(iT+1));
    NumOfEach(4,iT,iTrial)=sum(rest==4 & ~Branched & ~HasFormin);
    NumOfEach(5,iT,iTrial)=sum(rest==5 & ~Branched & ~HasFormin);
    NumOfEach(6,iT,iTrial)=Nfor-sum(HasFormin);
    NumOfEach(7,iT,iTrial)=sum(rest==2 & ~Branched & HasFormin);
    NumOfEach(8,iT,iTrial)=sum(rest==3 & ~Branched & HasFormin);
    NumOfEach(9,iT,iTrial)=sum(rest==4 & ~Branched & HasFormin);
    NumOfEach(10,iT,iTrial)=sum(rest==5 & ~Branched & HasFormin);
end
end
MeanNumOfEach(:,:,iError) = mean(NumOfEach,3);
end
MeanOfMean =  mean(MeanNumOfEach,3);
StdOfMean = std(MeanNumOfEach,0,3);
ts=0.5:0.5:tf;
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');
nexttile
set(gca,'ColorOrderIndex',1)
plot(tvals,yvals(:,1:nMax)/LBox^3,'-.')
hold on
%plot(tf*ones(5,1),Nums/LBox^3,'ko')
errorbar(ts,MeanOfMean(1:nMax,:)/LBox^3,StdOfMean(1:nMax,:)/LBox^3)
nexttile
set(gca,'ColorOrderIndex',1)
plot(tvals,yvals(:,nMax+1:2*nMax)/LBox^3,'-.')
hold on
%plot(tf*ones(5,1),Nums/LBox^3,'ko')
errorbar(ts,MeanOfMean(nMax+1:2*nMax,:)/LBox^3,StdOfMean(nMax+1:2*nMax,:)/LBox^3)

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

function dydt = RHS(t,y,RxnRates,nMax)
    dydt=zeros(length(y),1);
    PolyOn = RxnRates(5)+RxnRates(7);
    PolyOff = RxnRates(6)+RxnRates(8);
    % Stuff without attached formin
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
    % Formin bound things
    % nMax+1 = formin
    % nMax+2 = formin + dimer
    % nMax+n = formin + n-mer
    ForminNuc = RxnRates(9);
    ForminBind = RxnRates(10);
    Forminunbind = RxnRates(11);
    ForminEnhance = RxnRates(12);
    EnhancedPolyOn = ForminEnhance*RxnRates(5)+RxnRates(7);
    dydt(1) = dydt(1) - 2*ForminNuc*y(nMax+1)*y(1)*y(1);
    dydt(nMax+1) = -ForminNuc*y(nMax+1)*y(1)*y(1);
    dydt(nMax+2) = ForminNuc*y(nMax+1)*y(1)*y(1) -EnhancedPolyOn*y(nMax+2)*y(1)...
        + PolyOff*y(nMax+3);
    for iD=3:nMax-1
        dydt(nMax+iD) = EnhancedPolyOn*y(nMax+iD-1)*y(1)-PolyOff*y(nMax+iD) + ...
            PolyOff*y(nMax+iD+1)-EnhancedPolyOn*y(nMax+iD)*y(1);
    end
    iD=nMax;
    dydt(nMax+iD) = EnhancedPolyOn*y(nMax+iD-1)*y(1)-PolyOff*y(nMax+iD);
    for iD = 3:nMax
        dydt(1) = dydt(1) - EnhancedPolyOn*y(nMax+iD-1)*y(1) + PolyOff*y(nMax+iD);
    end
    % Formin binding and unbinding (from fibers 4 or larger)
    for iD=4:nMax
        dydt(iD)=dydt(iD)-ForminBind*y(iD)*y(nMax+1)+Forminunbind*y(nMax+iD);
        dydt(nMax+1)=dydt(nMax+1)-ForminBind*y(iD)*y(nMax+1)+Forminunbind*y(nMax+iD); % formin
        dydt(nMax+iD)=dydt(nMax+iD)+ForminBind*y(iD)*y(nMax+1)-Forminunbind*y(nMax+iD);
    end
end