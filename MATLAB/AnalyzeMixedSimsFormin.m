% Generate the following plots
% (1) Plot of free actin concentration over time
% (2) Histogram of length of fibers (in monomers)
Alphas = [0.5 1 2];
for iFc=1:length(Alphas)
alpha = Alphas(iFc);
Conc = 5; %uM
ForminConc = 0.12;
LBox = 5;
Vol = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
Nmon = floor(Conc*Vol/uMInvToMicron3);
NFormin = floor(ForminConc*1e-3*Vol/uMInvToMicron3);
tf = 2000;
dt = 1;
nT = tf/dt;
nError=5;
nTrial=10;
NumOfEach = zeros(6,nT,nError);
for iError=1:nError
FibLens = cell(nT,1);
for iTrial=1:nTrial
index = nTrial*(iError-1)+iTrial;
%if (alpha==1)
%StructInfo=load(strcat('StructInfo',num2str(Conc),'uM_Formin',num2str(ForminConc),...
%    'nM_',num2str(index),'.txt'));
%NumFibs=load(strcat('NumFibs',num2str(Conc),'uM_Formin',num2str(ForminConc),...
%    'nM_',num2str(index),'.txt'));
%BoundFormins=load(strcat('BoundFormin',num2str(Conc),'uM_Formin',num2str(ForminConc),...
%    'nM_',num2str(index),'.txt'));
%else
StructInfo=load(strcat('StructInfo',num2str(Conc),'uM_Formin0.12000000000000001',...
    'nM_Alpha',num2str(alpha),'_',num2str(index),'.txt'));
NumFibs=load(strcat('NumFibs',num2str(Conc),'uM_Formin0.12000000000000001',...
    'nM_Alpha',num2str(alpha),'_',num2str(index),'.txt'));
BoundFormins=load(strcat('BoundFormin',num2str(Conc),'uM_Formin0.12000000000000001',...
    'nM_Alpha',num2str(alpha),'_',num2str(index),'.txt'));
%end
TotalEntries = NumFibs+3;
StartIndex = [0;cumsum(TotalEntries)];
FibStartIndex = [0;cumsum(NumFibs)];
for iT=1:nT
    StartInd = StartIndex(iT)+1;
    EndInd = StartIndex(iT)+TotalEntries(iT);
    for iS=0:2
        NumOfEach(1+iS,iT,iError)=NumOfEach(1+iS,iT,iError)...
            +1/nTrial*StructInfo(StartInd+iS)/Vol*uMInvToMicron3;
    end
    NumOfEach(4,iT,iError)=NumOfEach(4,iT,iError)+1/nTrial*NumFibs(iT)/Vol;
    % Number of fibers with bound formins
    nBoundFormins = sum(BoundFormins(FibStartIndex(iT)+1:FibStartIndex(iT+1)));
    NumOfEach(5,iT,iError)=NumOfEach(5,iT,iError)+...
        1/nTrial*nBoundFormins/Vol;
    NumOfEach(6,iT,iError)=NumOfEach(6,iT,iError)...
        +1/nTrial*(NFormin-nBoundFormins)/Vol*uMInvToMicron3*1e3/ForminConc;
    FibLens{iT} = [FibLens{iT};StructInfo(StartInd+3:EndInd)];
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
ts=(1:nT)*dt;
figure(1)
set(gca,'ColorOrderIndex',iFc)
plot(ts,MeanOfEach(1,:))
hold on
skip=floor(nT/20);
set(gca,'ColorOrderIndex',iFc)
errorbar(ts(skip:skip:end),MeanOfEach(1,skip:skip:end),...
    2*StdOfEach(1,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.1)
ylabel('Concentration free monomers ($\mu$M)')
figure(2)
set(gca,'ColorOrderIndex',iFc)
plot(ts,MeanOfEach(6,:))
hold on
skip=floor(nT/20);
set(gca,'ColorOrderIndex',iFc)
errorbar(ts(skip/2:skip:end),MeanOfEach(6,skip/2:skip:end),...
    2*StdOfEach(6,skip/2:skip:end)/sqrt(nError),'o','MarkerSize',0.1)
ylabel('Fraction free formin')
figure(3)
set(gca,'ColorOrderIndex',iFc)
plot(ts,MeanOfEach(4,:))
hold on
skip=floor(nT/20);
set(gca,'ColorOrderIndex',iFc)
errorbar(ts(skip/2:skip:end),MeanOfEach(4,skip/2:skip:end),...
    2*StdOfEach(4,skip/2:skip:end)/sqrt(nError),'o','MarkerSize',0.1)
set(gca,'ColorOrderIndex',iFc)
plot(ts,MeanOfEach(5,:),':')
hold on
skip=floor(nT/20);
set(gca,'ColorOrderIndex',iFc)
errorbar(ts(skip:skip:end),MeanOfEach(5,skip:skip:end),...
    2*StdOfEach(5,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.1)
ylabel('Fibers per volume')
legend('All fibers','','Formin at barbed end')
% Histogram
MeanCountsByTime = cell(nT,1);
StdCountsByTime = cell(nT,1);
xpl = cell(nT,1);
for iT=1:nT
    MeanCountsByTime{iT}=mean(CountsByTime{iT});
    StdCountsByTime{iT}=std(CountsByTime{iT});
    xpl{iT}=(EdgesByTime{iT}(1:end-1)+EdgesByTime{iT}(2:end))*0.5;
end
figure(4)
subplot(1,3,iFc)
tspl =[500 1000 2000]/dt;
for iT=tspl
    errorbar(xpl{iT},MeanCountsByTime{iT},2*StdCountsByTime{iT}/sqrt(nError),'LineWidth',2.0);
    hold on
end
%legend(strcat('$t=$',num2str(tspl(1)*dt)),strcat('$t=$',num2str(tspl(2)*dt)),...
%    strcat('$t=$',num2str(tspl(3)*dt)))
%title(strcat(num2str(ForminConcs(iFc)),' nM Formin'))
title(strcat('$\alpha_\textrm{for}=$',num2str(alpha)))
end

