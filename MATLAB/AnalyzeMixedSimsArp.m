% Generate the following plots
% (1) Plot of free actin concentration over time
% (2) Histogram of length of fibers (in monomers)
ForminConcs = [0];
spacing = 2e-3;
for iC=1:length(ForminConcs)
if (iC==1)
    PConcs=0;
else
    PConcs=[2.5 5];
end
for iP=1:length(PConcs)
Conc = 1.5;
ConcProf =PConcs(iP);
ForminConc = ForminConcs(iC);
ArpConc= 20e-3;
LBox = 10;
MinForFiber = 4;
Vol = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
Nmon = floor(Conc*Vol/uMInvToMicron3);
NArp = floor(ArpConc*Vol/uMInvToMicron3);
tf = 10800;
dt = 10;
nT = tf/dt;
nError=3;
nTrial=10;
NumOfEach = zeros(6,nT,nError);
PercentBranched = zeros(nError,nT);
for iError=1:nError
LinFibLens = cell(nT,1);
BranchFibLens = cell(nT,1);
for iTrial=1:nTrial
index = nTrial*(iError-1)+iTrial;
FileName = strcat('Tf',num2str(tf),'_Box',num2str(LBox),'_Actin',num2str(Conc),...
    'uM_Prof',num2str(ConcProf),'uM_Arp',...
     num2str(ArpConc*1000),'nM_Formin',num2str(ForminConc*1e4),...
     'em4uM_',num2str(index),'.txt');
FreeMons = load(strcat('FreeMons',FileName));
StructInfo=load(strcat('StructInfo',FileName));
NumFibs=load(strcat('NumFibs',FileName));
Branched=load(strcat('BranchedOrLinear',FileName));
BoundFormins=load(strcat('BoundProteins',FileName));
FibStartIndex = [0;cumsum(NumFibs)];
for iT=1:nT
    NumOfEach(1,iT,iError)=NumOfEach(1,iT,iError)+1/nTrial*FreeMons(iT)/Vol*uMInvToMicron3;
    TheseLens = StructInfo(FibStartIndex(iT)+1:FibStartIndex(iT+1));
    if (sum(TheseLens)+FreeMons(iT)-Nmon ~=0)
        keyboard
    end
    TheseFormins = BoundFormins(FibStartIndex(iT)+1:FibStartIndex(iT+1));
    TheseBranched = Branched(FibStartIndex(iT)+1:FibStartIndex(iT+1));
    nFibs = NumFibs(iT);
    for iS=2:MinForFiber-1
        inds = TheseLens==iS;
        TheseLens(inds)=[];
        TheseFormins(inds)=[];
        TheseBranched(inds)=[];
        nFibs = nFibs - sum(inds);
    end
    NumOfEach(2,iT,iError)=NumOfEach(2,iT,iError)+1/nTrial*nFibs/Vol;
    % Number of linear fibers, branched fibers, and free arps 
    if (isempty(TheseBranched))
        TheseBranched=0;
        pBranch=0;
        NLin=0;
    else
        LinFibLens{iT} = [LinFibLens{iT};TheseLens(TheseBranched==0 | TheseBranched==2)];
        BranchFibLens{iT} = [BranchFibLens{iT};TheseLens(TheseBranched==1)];
        pBranch = sum(TheseBranched > 0)/length(TheseBranched);
        TotalLength = sum(TheseLens)*spacing;
        NLin=sum(TheseBranched==0 | TheseBranched==2);
        NumOfEach(3,iT,iError)=NumOfEach(3,iT,iError)+1/nTrial*(nFibs-NLin)/TotalLength;
        NumOfEach(6,iT,iError)=NumOfEach(6,iT,iError)+1/nTrial*sum(TheseFormins)/nFibs;
    end
    NumOfEach(4,iT,iError)=NumOfEach(4,iT,iError)+1/nTrial*NLin/Vol;
    FreeArps = (NArp-sum(TheseBranched==1))/NArp;
    NumOfEach(5,iT,iError)=NumOfEach(5,iT,iError)+1/nTrial*FreeArps;
    PercentBranched(iError,iT)=PercentBranched(iError,iT)+1/nTrial*pBranch;
end
end
if (iError==1)
    CountsByTime_Lin = cell(nT,1);
    EdgesByTime_Lin = cell(nT,1);
    CountsByTime_Br = cell(nT,1);
    EdgesByTime_Br = cell(nT,1);
    for iT=1:nT
        [CountsByTime_Lin{iT}(iError,:),EdgesByTime_Lin{iT}]=histcounts(LinFibLens{iT});
        dx = EdgesByTime_Lin{iT}(2)-EdgesByTime_Lin{iT}(1);
        CountsByTime_Lin{iT}(iError,:)=CountsByTime_Lin{iT}(iError,:)/(length(LinFibLens{iT})*dx);
        [CountsByTime_Br{iT}(iError,:),EdgesByTime_Br{iT}]=histcounts(BranchFibLens{iT});
        dx = EdgesByTime_Br{iT}(2)-EdgesByTime_Br{iT}(1);
        CountsByTime_Br{iT}(iError,:)=CountsByTime_Br{iT}(iError,:)/(length(BranchFibLens{iT})*dx);
    end
else
    for iT=1:nT
        CountsByTime_Lin{iT}(iError,:)=histcounts(LinFibLens{iT},EdgesByTime_Lin{iT});
        dx = EdgesByTime_Lin{iT}(2)-EdgesByTime_Lin{iT}(1);
        CountsByTime_Lin{iT}(iError,:)=CountsByTime_Lin{iT}(iError,:)/(length(LinFibLens{iT})*dx);
        CountsByTime_Br{iT}(iError,:)=histcounts(BranchFibLens{iT},EdgesByTime_Br{iT});
        dx = EdgesByTime_Br{iT}(2)-EdgesByTime_Br{iT}(1);
        CountsByTime_Br{iT}(iError,:)=CountsByTime_Br{iT}(iError,:)/(length(BranchFibLens{iT})*dx);
    end
end
end
% Take means and standard deviations
MeanOfEach = mean(NumOfEach,3);
StdOfEach = std(NumOfEach,0,3);
ts=(1:nT)*dt;
figure(1)
subplot(2,2,1)
set(gca,'ColorOrderIndex',iC)
LineSt='-';
plot(ts,Conc-MeanOfEach(1,:),LineSt)
hold on
skip=floor(nT/20);
try
set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1)
catch
set(gca,'ColorOrderIndex',1)
end
set(gca,'ColorOrderIndex',iC)
errorbar(ts(skip:skip:end),Conc-MeanOfEach(1,skip:skip:end),...
    2*StdOfEach(1,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.5)
title('Polymerized actin ($\mu$M)')
xlabel('$t$ (s)')
xlim([0 tf])
subplot(2,2,2)
plot(ts,MeanOfEach(6,:),LineSt)
hold on
try
set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1)
catch
set(gca,'ColorOrderIndex',1)
end
set(gca,'ColorOrderIndex',iC)
errorbar(ts(skip:skip:end),MeanOfEach(6,skip:skip:end),...
    2*StdOfEach(6,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.5)
xlabel('$t$ (s)')
xlim([0 tf])
title('Fraction formin-bound barbed ends')
subplot(2,2,3)
plot(ts,MeanOfEach(4,:),LineSt)
hold on
try
set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1)
catch
set(gca,'ColorOrderIndex',1)
end
set(gca,'ColorOrderIndex',iC)
errorbar(ts(skip:skip:end),MeanOfEach(4,skip:skip:end),...
    2*StdOfEach(4,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.5)
xlabel('$t$ (s)')
xlim([0 tf])
title('Linear fibers per volume')
subplot(2,2,4)
plot(ts,MeanOfEach(3,:),LineSt)
hold on
try
set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1)
catch
set(gca,'ColorOrderIndex',1)
end
set(gca,'ColorOrderIndex',iC)
errorbar(ts(skip:skip:end),MeanOfEach(3,skip:skip:end),...
    2*StdOfEach(3,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.5)
title('Branches per 1 $\mu$m')
xlabel('$t$ (s)')
xlim([0 tf])
figure(2)
plot(ts,MeanOfEach(5,:))
hold on
set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')-1)
errorbar(ts(skip:skip:end),MeanOfEach(5,skip:skip:end),...
    2*StdOfEach(5,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.5)
ylabel('Fraction free arp 2/3')
% Histogram
MeanCountsByTime = cell(nT,1);
StdCountsByTime = cell(nT,1);
xpl = cell(nT,1);
for iT=1:nT
    MeanCountsByTime{iT}=mean(CountsByTime_Lin{iT});
    StdCountsByTime{iT}=std(CountsByTime_Lin{iT});
    xpl{iT}=(EdgesByTime_Lin{iT}(1:end-1)+EdgesByTime_Lin{iT}(2:end))*0.5;
end
figure(3)
tspl =[100 300 600 1200 3600 10800]/dt;
subplot(1,length(ForminConcs),iC)
for iT=tspl
    errorbar(xpl{iT},MeanCountsByTime{iT},2*StdCountsByTime{iT}/sqrt(nError),'LineWidth',2.0);
    hold on
end
legend(strcat('$t=$',num2str(tspl(1)*dt)),strcat('$t=$',num2str(tspl(2)*dt)),...
    strcat('$t=$',num2str(tspl(3)*dt)),strcat('$t=$',num2str(tspl(4)*dt)))
title(strcat(num2str(ArpConc*1000),' nM'))
figure(4)
MeanCountsByTime = cell(nT,1);
StdCountsByTime = cell(nT,1);
xpl = cell(nT,1);
for iT=1:nT
    MeanCountsByTime{iT}=mean(CountsByTime_Br{iT});
    StdCountsByTime{iT}=std(CountsByTime_Br{iT});
    xpl{iT}=(EdgesByTime_Br{iT}(1:end-1)+EdgesByTime_Br{iT}(2:end))*0.5;
end
subplot(1,length(ForminConcs),iC)
for iT=tspl
    errorbar(xpl{iT},MeanCountsByTime{iT},2*StdCountsByTime{iT}/sqrt(nError),'LineWidth',2.0);
    hold on
end
legend(strcat('$t=$',num2str(tspl(1)*dt)),strcat('$t=$',num2str(tspl(2)*dt)),...
    strcat('$t=$',num2str(tspl(3)*dt)),strcat('$t=$',num2str(tspl(4)*dt)))
title(strcat(num2str(ForminConc*1000),' nM Formin'))
end
end
