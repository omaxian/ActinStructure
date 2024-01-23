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
kplusFor = 0*29.1; %uM^(-1)*s^(-1)
kminusFor = 0*8.1e-2; %s^(-1)
ForminEnhance = 1;

% Convert to microscopic assuming well-mixed system
Volume = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; % everything will be in s^(-1)
RxnRates=[kplusDimer*ConversionFactor; kminusDimer; kplusTrimer*ConversionFactor; ...
    kminusTrimer; kplusBarbed*ConversionFactor; kminusBarbed; ...
    kplusPointed*ConversionFactor; kminusPointed];
BarbedOnOff = [kplusFor*ConversionFactor kminusFor];
AlphaDimerMinus = [1 0];
AlphaDimerPlus = [1 kForNuc*ConversionFactor/kplusDimer];
AlphaTrimerMinus = [1 0];
AlphaTrimerPlus = [1 10];%(kplusBarbed+kplusPointed)/kplusTrimer];
AlphaBarbed = [1 ForminEnhance];
Nmon = floor(Conc*Volume/uMInvToMicron3);
Nfor = floor(ForminConc*Volume/uMInvToMicron3);
nMax=5;
% Solve the ODEs
RHSFcn = @(t,y) RHS(t,y,RxnRates,nMax,2,0,BarbedOnOff, ...
    AlphaDimerPlus,AlphaDimerMinus,AlphaTrimerPlus,AlphaTrimerMinus,AlphaBarbed);
y0 = [Nmon;zeros(nMax-1,1);Nfor;zeros(nMax-1,1)];
tf=40;
[tvals,yvals] = ode45(RHSFcn,[0 tf],y0);

% Import the data
nError=2;
nTrial=5;
NumFibs=load(strcat('ForminOnly_NumFibs1.txt'));
nT = length(NumFibs);
MeanNumOfEach = zeros(2*nMax,nT,nError);
for iError=1:nError
NumOfEach = zeros(2*nMax,nT,nTrial);
for iTrial=1:nTrial
TriIndex = (iError-1)*nTrial+iTrial;
FreeMons=load(strcat('ForminOnly_FreeMons',num2str(TriIndex),'.txt'));
StructInfo=load(strcat('ForminOnly_StructInfo',num2str(TriIndex),'.txt'));
NumFibs=load(strcat('ForminOnly_NumFibs',num2str(TriIndex),'.txt'));
BoundBarbed=load(strcat('ForminOnly_BoundFormins',num2str(TriIndex),'.txt'));
BranchedOrLinear=load(strcat('ForminOnly_BranchedOrLinear',num2str(TriIndex),'.txt'));
nT = length(NumFibs);
FibStarts = [0;cumsum(NumFibs)];
for iT=1:nT
    NumOfEach(1,:,iTrial)=FreeMons;
    FibNumbers = StructInfo(FibStarts(iT)+1:FibStarts(iT+1));
    Branched = BranchedOrLinear(FibStarts(iT)+1:FibStarts(iT+1));
    for iB=0:1
        MatchBarbed = BoundBarbed(FibStarts(iT)+1:FibStarts(iT+1))==iB;
        if (iB == 1)
            NumOfEach(5*iB+1,iT,iTrial) = Nfor-sum(MatchBarbed);
        end
        for Num=2:5
            NumOfEach(5*iB+Num,iT,iTrial)=sum(FibNumbers==Num & ~Branched & MatchBarbed);
        end
    end
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

function dydt = RHS(t,y,RxnRates,nMax,nBarbedProts,nMonProts,BarbedOnOff, ...
    AlphaDimerPlus,AlphaDimerMinus,AlphaTrimerPlus,AlphaTrimerMinus,AlphaBarbed)
    dydt=zeros(length(y),1);
    % 1 = free monomers
    % 2--5 = nothing on them
    % 5*nBarbed+1 = free barbed protein
    % 5*nBarbed+j = number with protein j on barbed end
    DimerPlus = RxnRates(1);
    DimerMinus = RxnRates(2);
    TrimerPlus = RxnRates(3);
    TrimerMinus = RxnRates(4);
    BrbPlus = RxnRates(5);
    BrbMinus = RxnRates(6);
    PtPlus = RxnRates(7);
    PtMinus = RxnRates(8);
    PolyOff = BrbMinus+PtMinus;
    FreeMon = y(1);
    for iB=1:nBarbedProts
        % Nucleation and polymerization
        PolyOn = BrbPlus*AlphaBarbed(iB)+PtPlus;
        if (iB==1)
            dydt((iB-1)*nMax+2)=dydt((iB-1)*nMax+2)+FreeMon*FreeMon*DimerPlus;
            dydt(1)=dydt(1)-2*FreeMon*FreeMon*DimerPlus;
        else
            dydt((iB-1)*nMax+2) = dydt((iB-1)*nMax+2)+FreeMon*FreeMon*DimerPlus...
                *AlphaDimerPlus(iB)*y((iB-1)*nMax+1);
            dydt(1) = dydt(1)-2*FreeMon*FreeMon*DimerPlus*AlphaDimerPlus(iB)*y((iB-1)*nMax+1);
            dydt((iB-1)*nMax+1) = dydt((iB-1)*nMax+1)-FreeMon*FreeMon*DimerPlus...
                *AlphaDimerPlus(iB)*y((iB-1)*nMax+1);
        end
        dydt((iB-1)*nMax+2) =  dydt((iB-1)*nMax+2) ...
            - DimerMinus*AlphaDimerMinus(iB)*y((iB-1)*nMax+2) ...
            - TrimerPlus*AlphaTrimerPlus(iB)*y((iB-1)*nMax+2)*FreeMon ...
            + TrimerMinus*AlphaTrimerMinus(iB)*y((iB-1)*nMax+3);
        dydt(1) = dydt(1)+2*DimerMinus*AlphaDimerMinus(iB)*y((iB-1)*nMax+2) ...
            -TrimerPlus*AlphaTrimerPlus(iB)*y((iB-1)*nMax+2)*FreeMon;
        dydt((iB-1)*nMax+3) = dydt((iB-1)*nMax+3) ...
            + TrimerPlus*AlphaTrimerPlus(iB)*y((iB-1)*nMax+2)*FreeMon ...
            - TrimerMinus*AlphaTrimerMinus(iB)*y((iB-1)*nMax+3) ...
            - PolyOn*y((iB-1)*nMax+3)*FreeMon + PolyOff*y((iB-1)*nMax+4);
        dydt(1) = dydt(1) + TrimerMinus*AlphaTrimerMinus(iB)*y((iB-1)*nMax+3);
        for iM=4:nMax-1
            dydt((iB-1)*nMax+iM) = dydt((iB-1)*nMax+iM) ...
                + PolyOn*FreeMon*(y((iB-1)*nMax+iM-1)-y((iB-1)*nMax+iM)) ...
                + PolyOff*(y((iB-1)*nMax+iM+1)-y((iB-1)*nMax+iM));
            dydt(1) = dydt(1) - PolyOn*FreeMon*y((iB-1)*nMax+iM-1)+ PolyOff*y((iB-1)*nMax+iM);
        end
        dydt((iB-1)*nMax+nMax) = dydt((iB-1)*nMax+nMax) ...
            + PolyOn*FreeMon*(y((iB-1)*nMax+nMax-1)) ...
            - PolyOff*y((iB-1)*nMax+nMax);
        dydt(1) = dydt(1) - PolyOn*FreeMon*y((iB-1)*nMax+nMax-1)+ PolyOff*y((iB-1)*nMax+nMax);
    end
end