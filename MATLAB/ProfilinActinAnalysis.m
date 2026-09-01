% Generate the following plots
% (1) Plot of free actin concentration over time
% (2) Histogram of length of fibers (in monomers)
figure(1)
Conc = 1.5;
ProfConc = 1.5;
KProf=1;
MinForFiber = 4;
tf = 28800;
LBox=10;
Vol = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
Nmon = floor(Conc*Vol/uMInvToMicron3);
dt = 10;
nT = tf/dt+1;
nError=2;
nTrial=3;
NumOfEach = zeros(5,nT,nError);
for iError=1:nError
FibLens = cell(nT,1);
nSimWithFib=zeros(1,nT);
for iTrial=1:nTrial
index = nTrial*(iError-1)+iTrial;
FileName = strcat('Actin',num2str(Conc),...
'uM_KProf',num2str(KProf),...
'_Prof',num2str(ProfConc),'uM_',num2str(index),'.txt');
FreeMons = load(strcat('FreeMonConc',FileName));
ProfActin = load(strcat('ProfBoundConc',FileName));
NumFibs=load(strcat('NumFibs',FileName));
NumOnEach=load(strcat('NumPerFib',FileName));
FibStartIndex = [0;cumsum(NumFibs)];
for iT=1:nT
    NumOfEach(1,iT,iError)=NumOfEach(1,iT,iError)+1/nTrial*FreeMons(iT)/Vol*uMInvToMicron3;
    NumOfEach(2,iT,iError)=NumOfEach(2,iT,iError)+1/nTrial*ProfActin(iT)/Vol*uMInvToMicron3;
    TheseLens = NumOnEach(FibStartIndex(iT)+1:FibStartIndex(iT+1));
    if (sum(TheseLens)+FreeMons(iT)+ProfActin(iT)-Nmon ~=0)
        keyboard
    end
    nFibs = NumFibs(iT);
    for iS=2:MinForFiber-1
        inds = TheseLens==iS;
        TheseLens(inds)=[];
        nFibs = nFibs - sum(inds);
    end
    pActin = sum(TheseLens);
    NumOfEach(3,iT,iError)=NumOfEach(3,iT,iError)+1/nTrial*pActin/Vol*uMInvToMicron3;
    NumOfEach(4,iT,iError)=NumOfEach(4,iT,iError)+1/nTrial*nFibs/Vol;
    FibLens{iT} = [FibLens{iT};TheseLens];
end
end
if (iError==1)
    CountsByTime = cell(nT,1);
    EdgesByTime = cell(nT,1);
    for iT=1:nT
        [CountsByTime{iT}(iError,:),EdgesByTime{iT}]=histcounts(FibLens{iT});
        dx = EdgesByTime{iT}(2)-EdgesByTime{iT}(1);
        CountsByTime{iT}(iError,:)=CountsByTime{iT}(iError,:)/(length(FibLens{iT})*dx);
    end
else
    for iT=1:nT
        CountsByTime{iT}(iError,:)=histcounts(FibLens{iT},EdgesByTime{iT});
        dx = EdgesByTime{iT}(2)-EdgesByTime{iT}(1);
        CountsByTime{iT}(iError,:)=CountsByTime{iT}(iError,:)/(length(FibLens{iT})*dx);
    end
end
end
% Take means and standard deviations
MeanOfEach = mean(NumOfEach,3);
StdOfEach = std(NumOfEach,0,3);
ts=(0:nT-1)*dt;
figure(1);
set(gca,'ColorOrderIndex',1)
for j=1:3
pltop=plot(ts,MeanOfEach(j,:));
hold on
skip=floor(nT/30);
skip=skip+mod(skip,2);
set(gca,'ColorOrderIndex',j)
pltop2=errorbar(ts(skip:skip:end),MeanOfEach(j,skip:skip:end),...
    2*StdOfEach(j,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.1);
end
legend('Free','','Profilin bound','','Polymerized')
xlim([0 28800])
xticks(0:7200:28800)
xticklabels(xticks/3600)
xlabel('Time (hr)')
ylabel('Concentration ($\mu$M)')
figure(2);
% Histogram
MeanCountsByTime = cell(nT,1);
StdCountsByTime = cell(nT,1);
xpl = cell(nT,1);
for iT=1:nT
    MeanCountsByTime{iT}=mean(CountsByTime{iT});
    StdCountsByTime{iT}=std(CountsByTime{iT});
    xpl{iT}=(EdgesByTime{iT}(1:end-1)+EdgesByTime{iT}(2:end))*0.5;
end
%subplot(2,3,RunIndex)
tspl =[3600 10800 21600 28800]/dt;
for iT=tspl
    errorbar(xpl{iT},MeanCountsByTime{iT},2*StdCountsByTime{iT}/sqrt(nError),'LineWidth',2.0);
    hold on
    xlabel('Number of monomers')
end
legend(strcat('$t=$',num2str(tspl(1)*dt)),strcat('$t=$',num2str(tspl(2)*dt)),...
    strcat('$t=$',num2str(tspl(3)*dt)),strcat('$t=$',num2str(tspl(4)*dt)))
title('Filament lengths (\# monomers)')