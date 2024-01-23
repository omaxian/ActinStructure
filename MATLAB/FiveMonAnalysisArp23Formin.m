% This validates dynamics of growing branches with formins. The only
% things NOT validated by this is depolymerization of branched networks.
% Parameters
Conc = 2; % in uM
LBox = 3;
ForminConc = 0.1;
ArpConc = 0;

% Parameters 
kplusDimer = 3.5e-3; % uM^(-1)*s^(-1) 
kminusDimer = 0.041; %s^(-1)
kplusTrimer = 13e-2; % uM^(-1)*s^(-1) 
kminusTrimer = 22; %s^(-1)
kplusBarbed = 1.6; % uM^(-1)*s^(-1) 
kminusBarbed = 0; %s^(-1)
kplusPointed = 1.3; %uM^(-1)*s^(-1)
kminusPointed = 0; %s^(-1)
kForNuc = 2e-3; % uM^(-2)*s^(-1)
kplusFor = 0*29.1; %uM^(-1)*s^(-1)
kminusFor = 0*8.1e-2; %s^(-1)
ForEnhance = 1;

kplusARF = 5.2e-2;
kminusARF = 0.5;

% Convert to microscopic assuming well-mixed system
Volume = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; % everything will be in s^(-1)
RxnRates=[kplusDimer*ConversionFactor; kminusDimer; kplusTrimer*ConversionFactor; ...
    kminusTrimer; kplusBarbed*ConversionFactor; kminusBarbed; ...
    kplusPointed*ConversionFactor; kminusPointed; ...
    kForNuc*ConversionFactor^2; kplusFor*ConversionFactor; kminusFor; ForEnhance; ...
    kplusARF*ConversionFactor^2;kminusARF];

Nmon = floor(Conc*Volume/uMInvToMicron3);
Nfor = floor(ForminConc*Volume/uMInvToMicron3);
ArpMon = floor(ArpConc*Volume/uMInvToMicron3);
nMax=5;
% Solve the ODEs
RHSFcn = @(t,y) RHS(t,y,RxnRates,nMax);
% 1-5 = unadorned actin
% 6-10 = formin attached 1-5
% 11-15 = branched actin 1-5
% 16 = arp 2/3 
% 17-18 = branched actin and formin 4+5
y0 = [Nmon;zeros(nMax-1,1);Nfor; zeros(nMax-1,1);zeros(nMax,1);ArpMon;0;0];
tf=40;
[tvals,yvals] = ode45(RHSFcn,[0 tf],y0);

% Import the data
nError=3;
nTrial=10;
NumFibs=load(strcat('ForminArpNoUnbind_NumFibs1.txt'));
nT = length(NumFibs);
MeanNumOfEach = zeros(18,nT,nError);
for iError=1:nError
NumOfEach = zeros(18,nT,nTrial);
for iTrial=1:nTrial
TriIndex = (iError-1)*nTrial+iTrial;
StructInfo=load(strcat('ForminArpNoUnbind_StructInfo',num2str(TriIndex),'.txt'));
NumFibs=load(strcat('ForminArpNoUnbind_NumFibs',num2str(TriIndex),'.txt'));
BoundFormins=load(strcat('ForminArpNoUnbind_BoundFormins',num2str(TriIndex),'.txt'));
BranchedOrLinear=load(strcat('ForminArpNoUnbind_BranchedOrLinear',num2str(TriIndex),'.txt'));
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
    Branched = BranchedOrLinear(FibStarts(iT)+1:FibStarts(iT+1))==1;
    HasFormin = BoundFormins(FibStarts(iT)+1:FibStarts(iT+1))>0;
    NumOfEach(4,iT,iTrial)=sum(rest==4 & ~Branched & ~HasFormin);
    NumOfEach(5,iT,iTrial)=sum(rest==5 & ~Branched & ~HasFormin);
    NumOfEach(6,iT,iTrial)=Nfor-sum(HasFormin);
    NumOfEach(7,iT,iTrial)=sum(rest==2 & ~Branched & HasFormin);
    NumOfEach(8,iT,iTrial)=sum(rest==3 & ~Branched & HasFormin);
    NumOfEach(9,iT,iTrial)=sum(rest==4 & ~Branched & HasFormin);
    NumOfEach(10,iT,iTrial)=sum(rest==5 & ~Branched & HasFormin);
    NumOfEach(11,iT,iTrial)=sum(rest==1 & Branched & ~HasFormin);
    NumOfEach(12,iT,iTrial)=sum(rest==2 & Branched & ~HasFormin);
    NumOfEach(13,iT,iTrial)=sum(rest==3 & Branched & ~HasFormin);
    NumOfEach(14,iT,iTrial)=sum(rest==4 & Branched & ~HasFormin);
    NumOfEach(15,iT,iTrial)=sum(rest==5 & Branched & ~HasFormin);
    NumOfEach(16,iT,iTrial)=ArpMon-sum(Branched);
    NumOfEach(17,iT,iTrial)=sum(rest==4 & Branched & HasFormin);
    NumOfEach(18,iT,iTrial)=sum(rest==5 & Branched & HasFormin);
end
end
MeanNumOfEach(:,:,iError) = mean(NumOfEach,3);
end
MeanOfMean =  mean(MeanNumOfEach,3);
StdOfMean = std(MeanNumOfEach,0,3);
ts=0.5:0.5:tf;
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
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
nexttile
set(gca,'ColorOrderIndex',1)
plot(tvals,yvals(:,2*nMax+1:end)/LBox^3,'-.')
hold on
%plot(tf*ones(5,1),Nums/LBox^3,'ko')
errorbar(ts,MeanOfMean(2*nMax+1:end,:)/LBox^3,StdOfMean(2*nMax+1:end,:)/LBox^3)


function dydt = RHS(t,y,RxnRates,nMax)
    % Attachment ONLY 
    dydt=zeros(length(y),1);
    BarbedPlus = RxnRates(5);
    PolyOn = RxnRates(5)+RxnRates(7);
    
    % Step 1: actin by itself
    dydt(1) = -2*RxnRates(1)*y(1).^2 - RxnRates(3)*y(1)*y(2)...
        + 2*RxnRates(2)*y(2)+ RxnRates(4)*y(3); % monomers
    for iD = 3:nMax-1
        dydt(1) = dydt(1) - PolyOn*y(iD)*y(1);
    end
    dydt(2) = RxnRates(1)*y(1).^2 - RxnRates(3)*y(1)*y(2) ...
        - RxnRates(2)*y(2) + RxnRates(4)*y(3);% dimers
    % Trimers: form, unform, form tetramers, unform tetramers
    dydt(3) = RxnRates(3)*y(1)*y(2) - PolyOn*y(1)*y(3) ...
        - RxnRates(4)*y(3); % trimers
    dydt(4) = PolyOn*y(1)*(y(3)-y(4));
    dydt(5) = PolyOn*y(1)*y(4);

    % Step 2: adding formin
    ForminNuc = RxnRates(9);
    ForminBind = RxnRates(10);
    Forminunbind = RxnRates(11);
    ForminEnhance = RxnRates(12);
    EnhancedPolyOn = ForminEnhance*RxnRates(5)+RxnRates(7);
    dydt(1) = dydt(1) - 2*ForminNuc*y(nMax+1)*y(1)*y(1) ...
        - EnhancedPolyOn*(y(nMax+2)+y(nMax+3)+y(nMax+4))*y(1);
    dydt(nMax+1) = -ForminNuc*y(nMax+1)*y(1)*y(1);
    dydt(nMax+2) = ForminNuc*y(nMax+1)*y(1)*y(1) -EnhancedPolyOn*y(nMax+2)*y(1);
    dydt(nMax+3) = EnhancedPolyOn*y(1)*(y(nMax+2)-y(nMax+3));
    dydt(nMax+4) = EnhancedPolyOn*y(1)*(y(nMax+3)-y(nMax+4));
    dydt(nMax+5) = EnhancedPolyOn*y(1)*y(nMax+4);
    for iD=4:nMax
        dydt(iD)=dydt(iD)-ForminBind*y(iD)*y(nMax+1)+Forminunbind*y(nMax+iD);
        dydt(nMax+1)=dydt(nMax+1)-ForminBind*y(iD)*y(nMax+1)+Forminunbind*y(nMax+iD);
        dydt(nMax+iD)=dydt(nMax+iD)+ForminBind*y(iD)*y(nMax+1)-Forminunbind*y(nMax+iD);
    end

    % Arp 2/3 dynamics
    ARFPlus = RxnRates(13);
    ARFMinus = RxnRates(14);
    BranchElig=y(4)+y(5)+y(nMax+4)+y(nMax+5)+y(2*nMax+4)+y(2*nMax+5)+y(3*nMax+2)+y(3*nMax+3);
    % Branch formation
    dydt(1) = dydt(1) - ARFPlus*y(1)*y(3*nMax+1)*BranchElig + ARFMinus*y(2*nMax+1) ...
        - BarbedPlus*y(1)*sum(y(2*nMax+1:2*nMax+4));
    dydt(3*nMax+1) =  - ARFPlus*y(1)*y(3*nMax+1)*BranchElig + ARFMinus*y(2*nMax+1);
    % Branch polymerization
    dydt(2*nMax+1) = ARFPlus*y(1)*y(3*nMax+1)*BranchElig - BarbedPlus*y(2*nMax+1)*y(1) ...
         - ARFMinus*y(2*nMax+1);
    dydt(2*nMax+2) = BarbedPlus*y(1)*(y(2*nMax+1)-y(2*nMax+2));
    dydt(2*nMax+3) = BarbedPlus*y(1)*(y(2*nMax+2)-y(2*nMax+3));
    dydt(2*nMax+4) = BarbedPlus*y(1)*(y(2*nMax+3)-y(2*nMax+4));
    dydt(2*nMax+5) = BarbedPlus*y(1)*y(2*nMax+4);


    % Branches and formins (come back to this later)
    % 3*nMax+1 = Arp 2/3
    % 3*nMax+2 = formin branched fiber length 4
    % 3*nMax+3 = formin branched fiber length 5
    dydt(2*nMax+4)=dydt(2*nMax+4)-ForminBind*y(nMax+1)*y(2*nMax+4)+Forminunbind*y(3*nMax+2);
    dydt(2*nMax+5)=dydt(2*nMax+5)-ForminBind*y(nMax+1)*y(2*nMax+5)+Forminunbind*y(3*nMax+3);
    dydt(3*nMax+2) = dydt(3*nMax+2) - ForminEnhance*BarbedPlus*y(1)*y(3*nMax+2) ...
        -Forminunbind*y(3*nMax+2) + ForminBind*y(nMax+1)*y(2*nMax+4); %R4
    dydt(3*nMax+3) = dydt(3*nMax+3) + ForminEnhance*BarbedPlus*y(1)*y(3*nMax+2) ...
        -Forminunbind*y(3*nMax+3) + ForminBind*y(nMax+1)*y(2*nMax+5); %R5
    dydt(1) = dydt(1) -  ForminEnhance*BarbedPlus*y(1)*y(3*nMax+2);
    dydt(nMax+1) = dydt(nMax+1)+Forminunbind*(y(3*nMax+2)+y(3*nMax+3))...
        - ForminBind*y(nMax+1)*(y(2*nMax+4)+y(2*nMax+5));
end