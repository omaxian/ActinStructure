% This validates dynamics of formin without branching
% Parameters
Conc = 2; % in uM
LBox = 3;
ForminConc = 0.1;
ProfConc = 0.05;
ArpConc = 0.1;

% Parameters 
kplusDimer = 3.5e-3; % uM^(-1)*s^(-1) 
kminusDimer = 0.041; %s^(-1)
kplusTrimer = 13e-2; % uM^(-1)*s^(-1) 
kminusTrimer = 22; %s^(-1)
kplusBarbed = 1.6; % uM^(-1)*s^(-1) 
kminusBarbed = 0*1.4; %s^(-1)
kplusPointed = 1.3; %uM^(-1)*s^(-1)
kminusPointed = 0*0.8; %s^(-1)
kForNuc = 2e-3; % uM^(-2)*s^(-1)
kplusFor = 5; %uM^(-1)*s^(-1)
kminusFor = 8.1e-2; %s^(-1)
ForminEnhance = 2;

% Convert to microscopic assuming well-mixed system
Volume = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
ConversionFactor = uMInvToMicron3/Volume; % everything will be in s^(-1)
RxnRates=[kplusDimer*ConversionFactor; kminusDimer; kplusTrimer*ConversionFactor; ...
    kminusTrimer; kplusBarbed*ConversionFactor; kminusBarbed; ...
    kplusPointed*ConversionFactor; kminusPointed];
BarbedOnOff = [kplusFor*ConversionFactor kminusFor];

AlphaBarbedMinus = [1 0.1];
AlphaDimerMinus = [1 0.2];
AlphaDimerPlus = [1 kForNuc*ConversionFactor/kplusDimer; ...
    1.7 0.5*kForNuc*ConversionFactor/kplusDimer];
AlphaTrimerMinus = [1 0.4];
AlphaTrimerPlus = [1 10; 0.4 2];
AlphaBarbedPlus = [1 ForminEnhance; 2 2.5];
AlphaPointed = [1;0.5];
MonEqs = 5*ConversionFactor;
Nmon = floor(Conc/ConversionFactor);
Nfor = floor(ForminConc/ConversionFactor);
Narp = floor(ArpConc/ConversionFactor);
NBarbed = Nfor;
nMax=5;
BranchRates = [10e-2*ConversionFactor^2 0.5];
AlphaBranchOn = [1; 1.3];
% Solve the ODEs
nBarbedProts = ForminConc>0;
nMonProts = floor(ProfConc*Volume/uMInvToMicron3);
RHSFcn = @(t,y) RHS(t,y,RxnRates,nMax,nBarbedProts,nMonProts,BarbedOnOff,MonEqs, ...
    AlphaDimerPlus,AlphaDimerMinus,AlphaTrimerPlus,AlphaTrimerMinus,AlphaBarbedPlus,...
    AlphaBarbedMinus,AlphaPointed,BranchRates,AlphaBranchOn);
y0 = zeros((nBarbedProts+1)*nMax,1);
y0(1) = Nmon;
for iB=1:nBarbedProts
    y0(iB*nMax+1)=NBarbed(iB);
end
if (ArpConc > 0)
    y0 = [y0; Narp; zeros((nBarbedProts+1)*nMax,1)];
end
tf=40;
[tvals,yvals] = ode45(RHSFcn,[0 tf],y0);

% Import the data
nError=2;
nTrial=5;
NumFibs=load(strcat('ForminOnly_NumFibs1.txt'));
nT = length(NumFibs);
MeanNumOfEach = zeros((nBarbedProts+1)*nMax,nT,nError);
if (ArpConc > 0)
    MeanNumOfEach = zeros(2*(nBarbedProts+1)*nMax+1,nT,nError);
end
for iError=1:nError
NumOfEach = zeros((nBarbedProts+1)*nMax,nT,nTrial);
if (ArpConc > 0)
    NumOfEach = zeros(2*(nBarbedProts+1)*nMax+1,nT,nTrial);
end
for iTrial=1:nTrial
TriIndex = (iError-1)*nTrial+iTrial
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
    Branched = BranchedOrLinear(FibStarts(iT)+1:FibStarts(iT+1))==1;
    for iB=0:nBarbedProts
        MatchBarbed = BoundBarbed(FibStarts(iT)+1:FibStarts(iT+1))==iB;
        if (iB > 0)
            NumOfEach(5*iB+1,iT,iTrial) = NBarbed(iB)-sum(MatchBarbed);
        end
        for Num=2:5
            NumOfEach(5*iB+Num,iT,iTrial)=sum(FibNumbers==Num & ~Branched & MatchBarbed);
        end
    end
    if (ArpConc > 0)
        FirstIndex = nMax*(nBarbedProts+1)+1;
        for iB=0:nBarbedProts
            MatchBarbed = BoundBarbed(FibStarts(iT)+1:FibStarts(iT+1))==iB;
            for Num=1:5
                NumOfEach(FirstIndex+5*iB+Num,iT,iTrial)=...
                    sum(FibNumbers==Num & Branched & MatchBarbed);
            end
        end
        NumOfEach(FirstIndex,iT,iTrial)=Narp - sum(Branched);
    end
end
end
MeanNumOfEach(:,:,iError) = mean(NumOfEach,3);
end
MeanOfMean =  mean(MeanNumOfEach,3);
StdOfMean = std(MeanNumOfEach,0,3);
ts=0.5:0.5:tf;
if (ArpConc==0)
tiledlayout(1,nBarbedProts+1,'Padding', 'none', 'TileSpacing', 'compact');
else
tiledlayout(2,nBarbedProts+1,'Padding', 'none', 'TileSpacing', 'compact');
end
for iB=0:nBarbedProts
    nexttile
    set(gca,'ColorOrderIndex',1)
    inds = iB*nMax+1:(iB+1)*nMax;
    x=get(gca,'ColorOrderIndex');
    plot(tvals,yvals(:,inds)/LBox^3,'-.')
    hold on
    %plot(tf*ones(5,1),Nums/LBox^3,'ko')
    set(gca,'ColorOrderIndex',x)
    errorbar(ts,MeanOfMean(inds,:)/LBox^3,StdOfMean(inds,:)/LBox^3)
end
if (ArpConc > 0)
    for iB=0:nBarbedProts
    nexttile
    set(gca,'ColorOrderIndex',1)
    StartIndex = (nBarbedProts+1)*nMax+1;
    inds = StartIndex+iB*nMax+1:StartIndex+(iB+1)*nMax;
    if (iB==0)
        inds=[inds StartIndex];
    end
    x=get(gca,'ColorOrderIndex');
    plot(tvals,yvals(:,inds)/LBox^3,'-.')
    hold on
    %plot(tf*ones(5,1),Nums/LBox^3,'ko')
    set(gca,'ColorOrderIndex',x)
    errorbar(ts,MeanOfMean(inds,:)/LBox^3,StdOfMean(inds,:)/LBox^3)
    end
end


function dydt = RHS(t,y,RxnRates,nMax,nBarbedProts,NumMonProts,BarbedOnOff,MonEqConsts, ...
    AlphaDimerPlus,AlphaDimerMinus,AlphaTrimerPlus,AlphaTrimerMinus,AlphaBarbedPlus, ...
    AlphaBarbedMinus, AlphaPointed,BranchRates,AlphaBranchOn)
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
    FreeMon = y(1);
    BoundToEach = FreeMon;
    % Determine percentage of monomer bound to different proteins
    if (~isempty(NumMonProts))
        EqSum = sum(MonEqConsts.*NumMonProts);
        BoundToEach=[FreeMon; FreeMon.*NumMonProts.*MonEqConsts]/(1+EqSum);
    end
    for iB=1:nBarbedProts+1
        % Nucleation and polymerization
        DimerPlusRate = DimerPlus*sum(BoundToEach.*BoundToEach.*AlphaDimerPlus(:,iB));
        DimerMinusRate = DimerMinus*AlphaDimerMinus(iB)*y((iB-1)*nMax+2);
        TrimerPlusRate = TrimerPlus*sum(BoundToEach.*AlphaTrimerPlus(:,iB))*y((iB-1)*nMax+2);
        TrimerMinusRate = TrimerMinus*AlphaTrimerMinus(iB)*y((iB-1)*nMax+3);
        PolyOn = sum((PtPlus*AlphaPointed+BrbPlus*AlphaBarbedPlus(:,iB)).*BoundToEach);
        PolyOff = BrbMinus*AlphaBarbedMinus(:,iB)+PtMinus;
        if (iB > 1)
            DimerPlusRate = DimerPlusRate*y((iB-1)*nMax+1);
            dydt((iB-1)*nMax+1) = dydt((iB-1)*nMax+1)-DimerPlusRate + DimerMinusRate;
        end
        dydt((iB-1)*nMax+2) =  dydt((iB-1)*nMax+2) + DimerPlusRate - DimerMinusRate ...
            - TrimerPlusRate + TrimerMinusRate;
        dydt(1) = dydt(1)-2*DimerPlusRate+2*DimerMinusRate -TrimerPlusRate + TrimerMinusRate;
        dydt((iB-1)*nMax+3) = dydt((iB-1)*nMax+3) + TrimerPlusRate  - TrimerMinusRate ...
            - PolyOn*y((iB-1)*nMax+3) + PolyOff*y((iB-1)*nMax+4);
        for iM=4:nMax-1
            dydt((iB-1)*nMax+iM) = dydt((iB-1)*nMax+iM) ...
                + PolyOn*(y((iB-1)*nMax+iM-1)-y((iB-1)*nMax+iM)) ...
                + PolyOff*(y((iB-1)*nMax+iM+1)-y((iB-1)*nMax+iM));
            dydt(1) = dydt(1) - PolyOn*y((iB-1)*nMax+iM-1)+ PolyOff*y((iB-1)*nMax+iM);
        end
        dydt((iB-1)*nMax+nMax) = dydt((iB-1)*nMax+nMax) ...
            + PolyOn*(y((iB-1)*nMax+nMax-1)) - PolyOff*y((iB-1)*nMax+nMax);
        dydt(1) = dydt(1) - PolyOn*y((iB-1)*nMax+nMax-1)+ PolyOff*y((iB-1)*nMax+nMax);
    end
    for iB=2:nBarbedProts+1
        for j = 2:nMax
            dBind = -BarbedOnOff(2*(iB-1))*y((iB-1)*nMax+j) ...
                + BarbedOnOff(2*(iB-1)-1)*y(j)*y((iB-1)*nMax+1);
            dydt(j)=dydt(j)- dBind;
            dydt((iB-1)*nMax+1) = dydt((iB-1)*nMax+1) - dBind;
            dydt((iB-1)*nMax+j) = dydt((iB-1)*nMax+j) + dBind;
        end
    end
    % Adjustments for branching
    if (length(y) > (nBarbedProts+1)*nMax)
        % There are 2 more entries: one for free branchers and nBarbedProts+1 for
        % branches of length 1
        ArpInd = (nBarbedProts+1)*nMax+1;
        nArp = y(ArpInd);
        IndsForBranching4 = [4:nMax:(nBarbedProts+1)*nMax ArpInd+4:nMax:length(y)];
        IndsForBranching5 = [5:nMax:(nBarbedProts+1)*nMax ArpInd+5:nMax:length(y)];
        NewBranchRate = BranchRates(1)*sum(AlphaBranchOn.*BoundToEach)*nArp...
            *(4*sum(y(IndsForBranching4))+5*sum(y(IndsForBranching5)))...
            - BranchRates(2)*y(ArpInd+1);
        dydt(1) = dydt(1) - NewBranchRate;
        dydt(ArpInd) = dydt(ArpInd) - NewBranchRate;
        dydt(ArpInd+1) = dydt(ArpInd+1) + NewBranchRate;
        % Bulk unbinding of protein-bound
%         for iB=2:nBarbedProts+1
%             UnbindRate = BranchRates(2)*y(ArpInd+(iB-1)*nMax+1);
%             dydt(1) = dydt(1) + UnbindRate;
%             dydt(ArpInd) = dydt(ArpInd) + UnbindRate;
%             dydt(ArpInd+(iB-1)*nMax+1) = dydt(ArpInd+(iB-1)*nMax+1) - UnbindRate;
%         end
        % Binding of barbed proteins to branches
        for iB=2:nBarbedProts+1
            for j = 1:nMax
                dBind = -BarbedOnOff(2*(iB-1))*y(ArpInd+(iB-1)*nMax+j) ...
                    + BarbedOnOff(2*(iB-1)-1)*y(ArpInd+j)*y((iB-1)*nMax+1);
                dydt(ArpInd+j) = dydt(ArpInd+j) - dBind;
                dydt((iB-1)*nMax+1) = dydt((iB-1)*nMax+1) - dBind;
                dydt(ArpInd+(iB-1)*nMax+j) = dydt(ArpInd+(iB-1)*nMax+j) + dBind;
            end
        end
        % Growth of branches
        for iB=1:nBarbedProts+1
            PolyOn = sum((BrbPlus*AlphaBarbedPlus(:,iB)).*BoundToEach);
            PolyOff = BrbMinus*AlphaBarbedMinus(:,iB);
            for iM=1:nMax-1
                index = ArpInd+(iB-1)*nMax+iM;
                flux = PolyOn*y(index) - PolyOff*y(index+1);
                dydt(index) = dydt(index) - flux;
                dydt(index+1) = dydt(index+1) + flux;
                dydt(1) = dydt(1) - flux;
            end
        end   
    end 
end