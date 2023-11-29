% Generate the following plots
% (1) Plot of free actin concentration over time
% (2) Histogram of length of fibers (in monomers)
Conc = 10;
LBox = 5;
Vol = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
Nmon = floor(Conc*Vol/uMInvToMicron3);
tf = 10800*(Conc==2)+5400*(Conc==5)+2000*(Conc==10);
dt = 2;
nT = tf/dt;
nError=5;
nTrial=10;
NumOfEach = zeros(4,nT,nError);
for iError=1:nError
FibLens = cell(nT,1);
for iTrial=1:nTrial
index = nTrial*(iError-1)+iTrial;
StructInfo=load(strcat('StructInfo',num2str(Conc),'uM_',num2str(index),'.txt'));
NumFibs=load(strcat('NumFibs',num2str(Conc),'uM_',num2str(index),'.txt'));
TotalEntries = NumFibs+3;
StartIndex = [0;cumsum(TotalEntries)];
for iT=1:nT
    StartInd = StartIndex(iT)+1;
    EndInd = StartIndex(iT)+TotalEntries(iT);
    for iS=0:2
        NumOfEach(1+iS,iT,iError)=NumOfEach(1+iS,iT,iError)...
            +1/nTrial*StructInfo(StartInd+iS)/Vol*uMInvToMicron3;
    end
    NumOfEach(4,iT,iError)=NumOfEach(4,iT,iError)+1/nTrial*NumFibs(iT)/Vol;
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
plot(ts,MeanOfEach(1,:))
hold on
skip=floor(nT/20);
errorbar(ts(skip:skip:end),MeanOfEach(1,skip:skip:end),...
    2*StdOfEach(1,skip:skip:end)/sqrt(nError),'o')
figure(2)
plot(ts,MeanOfEach(4,:))
hold on
skip=floor(nT/20);
errorbar(ts(skip:skip:end),MeanOfEach(4,skip:skip:end),...
    2*StdOfEach(4,skip:skip:end)/sqrt(nError),'o')
% Histogram
MeanCountsByTime = cell(nT,1);
StdCountsByTime = cell(nT,1);
xpl = cell(nT,1);
for iT=1:nT
    MeanCountsByTime{iT}=mean(CountsByTime{iT});
    StdCountsByTime{iT}=std(CountsByTime{iT});
    xpl{iT}=(EdgesByTime{iT}(1:end-1)+EdgesByTime{iT}(2:end))*0.5;
end
figure(3)
if (Conc==10)
tspl =[100 500 2000]/dt;
subplot(1,3,3)
elseif (Conc==5)
tspl =[500 1000 5000]/dt;
subplot(1,3,2)  
elseif (Conc==2)
tspl =[2000 5400 10800]/dt;
subplot(1,3,1)  
end
for iT=tspl
    errorbar(xpl{iT},MeanCountsByTime{iT},2*StdCountsByTime{iT}/sqrt(nError),'LineWidth',2.0);
    hold on
end
legend(strcat('$t=$',num2str(tspl(1)*dt)),strcat('$t=$',num2str(tspl(2)*dt)),...
    strcat('$t=$',num2str(tspl(3)*dt)))

