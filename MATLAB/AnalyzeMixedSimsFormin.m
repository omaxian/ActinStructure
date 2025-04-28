% Generate the following plots
% (1) Plot of free actin concentration over time
% (2) Histogram of length of fibers (in monomers)
figure(1)
tiledlayout(1,4,'Padding', 'none', 'TileSpacing', 'compact');
ProfConcs = [0 0.75 1.5 6]; % Do not change this; change the indices
SeedConcs= [0 0.001 0.002];
ConcFormins = [0];
LineStys = ["--","--","-."];
for iP=1:length(ProfConcs)
figure(1)
nexttile
for iK=1:length(SeedConcs)
figure(1)
RunIndex = 1;%length(ConcFormins)*(iC-1)+iFc;
Conc = 1.5; %uM
SeedConc=SeedConcs(iK);
ConcProf = ProfConcs(iP);
KProf=1;%KProfs(iK);
ArpConc = 160e-3;
ForminConc = 0;
MinForFiber = 4;
tf = 28800;
if (ConcProf<1.5)
    LBox =5;
else
    LBox=10;
end
Vol = LBox^3;
uMInvToMicron3 = 1.0e15/(6.022e17);
Nmon = floor(Conc*Vol/uMInvToMicron3);
NFormin = floor(ForminConc*Vol/uMInvToMicron3);
dt = 10;
nT = tf/dt;
nError=3;
nTrial=5;
NumOfEach = zeros(5,nT,nError);
for iError=1:nError
FibLens = cell(nT,1);
nSimWithFib=zeros(1,nT);
for iTrial=1:nTrial
index = nTrial*(iError-1)+iTrial;
if (ConcProf>-1)
    FileName = strcat('Tf',num2str(tf),'_Box',num2str(LBox),'_Actin',num2str(Conc),...
    'uM_Seed_',num2str(SeedConc),'_KProf',num2str(KProf),...
    '_Prof',num2str(ConcProf),'uM_Arp',...
     num2str(ArpConc*1000),'nM_Formin',num2str(ForminConc*1e4),...
     'em4uM_',num2str(index),'.txt');
else
    FileName = strcat('Tf10800_Box10_Actin1.5uM_Prof0uM_Arp0nM_Formin0em4uM_',...
        num2str(index),'.txt');
end
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
    nFibs = NumFibs(iT);
    for iS=2:MinForFiber-1
        inds = TheseLens==iS;
        TheseLens(inds)=[];
        TheseFormins(inds)=[];
        nFibs = nFibs - sum(inds);
    end
    NumOfEach(2,iT,iError)=NumOfEach(2,iT,iError)+1/nTrial*nFibs/Vol;
    % Number of fibers with bound formins
    nBoundFormins = sum(TheseFormins);
    NumOfEach(3,iT,iError)=NumOfEach(3,iT,iError)+1/nTrial*nBoundFormins/Vol;
    if (ForminConc>0)
        NumOfEach(4,iT,iError)=NumOfEach(4,iT,iError)...
            +1/nTrial*(NFormin-nBoundFormins)/NFormin;
    end
    if (sum(TheseLens)==0)
        NumOfEach(5,iT,iError)=NumOfEach(5,iT,iError)+0;
    else
        nSimWithFib(iT)=nSimWithFib(iT)+1;
        NumOfEach(5,iT,iError)=NumOfEach(5,iT,iError)+...
         sum(TheseLens(TheseFormins==1))/sum(TheseLens);
    end
    FibLens{iT} = [FibLens{iT};TheseLens];
end
end
NumOfEach(5,:,iError)=NumOfEach(5,:,iError)./nSimWithFib;
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
%subplot(1,4,iP)
set(gca,'ColorOrderIndex',1)
LineSty = LineStys(iK);
pltop=plot(ts,Conc-MeanOfEach(1,:),LineSty);
hold on
skip=floor(nT/30);
skip=skip+mod(skip,2);
set(gca,'ColorOrderIndex',1)
pltop2=errorbar(ts(skip:skip:end),Conc-MeanOfEach(1,skip:skip:end),...
    2*StdOfEach(1,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.1);
title(strcat(num2str(ProfConcs(iP)), ' $\mu$M profilin'))
if (iK==1)
if (ArpConc==0)
if (ConcProf==0)
    plot(Times1,AllData_Actin1p5Prof0)
elseif (ConcProf==0.75)
    plot(Times1,AllData_Actin1p5Prof0p75)
elseif (ConcProf==1.5)
    plot(Times1,AllData_Actin1p5Prof1p5)
elseif (ConcProf==6)
    plot(Times1,AllData_Actin1p5Prof6)
end
elseif (ArpConc==0.08)
if (ConcProf==0)
    plot(Times1,ConvertedData_June5(:,5))
    plot(Times1,ConvertedData_June12_Int(:,5))
elseif (ConcProf==0.75)
    plot(Times1,ConvertedData_June5(:,6))
    plot(Times1,ConvertedData_June12_Int(:,6))
elseif (ConcProf==1.5)
    plot(Times1,ConvertedData_June5(:,7))
    plot(Times1,ConvertedData_June12_Int(:,7))
elseif (ConcProf==6)
    plot(Times1,ConvertedData_June5(:,8))
    plot(Times1,ConvertedData_June12_Int(:,8))
end
elseif (ArpConc==0.16)
if (ConcProf==0)
    plot(Times1,ConvertedData_June5(:,9))
    plot(Times1,ConvertedData_June12_Int(:,9))
elseif (ConcProf==0.75)
    plot(Times1,ConvertedData_June5(:,10))
    plot(Times1,ConvertedData_June12_Int(:,10))
elseif (ConcProf==1.5)
    plot(Times1,ConvertedData_June5(:,11))
    plot(Times1,ConvertedData_June12_Int(:,11))
elseif (ConcProf==6)
    plot(Times1,ConvertedData_June5(:,12))
    plot(Times1,ConvertedData_June12_Int(:,12))
end
end
end
ylim([0 1.5])
xlabel('$t$ (hr)')
xlim([0 28800])
uistack(pltop2,'top')
uistack(pltop,'top')
xticks(0:7200:28800)
xticklabels(xticks/3600)
if (iP==1)
    ylabel('Polymerized actin ($\mu$M)')
end
%ylabel('Poly ($\mu$M)')
% figure(2)
% set(gca,'ColorOrderIndex',iFc)
% plot(ts,MeanOfEach(4,:),LineSty)
% hold on
% set(gca,'ColorOrderIndex',iFc)
% errorbar(ts(skip/2:skip:end),MeanOfEach(4,skip/2:skip:end),...
%     2*StdOfEach(4,skip/2:skip:end)/sqrt(nError),'o','MarkerSize',0.1)
% ylabel('Fraction free formin')
figure(3)
set(gca,'ColorOrderIndex',iP)
plot(ts,MeanOfEach(2,:),LineSty)
hold on
set(gca,'ColorOrderIndex',iP)
errorbar(ts(skip/2:skip:end),MeanOfEach(2,skip/2:skip:end),...
    2*StdOfEach(2,skip/2:skip:end)/sqrt(nError),'o','MarkerSize',0.1)
% set(gca,'ColorOrderIndex',iP)
% plot(ts,MeanOfEach(3,:),':')
% hold on
% set(gca,'ColorOrderIndex',iFc)
% errorbar(ts(skip:skip:end),MeanOfEach(3,skip:skip:end),...
%    2*StdOfEach(3,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.1)
ylabel('Fibers per volume')
% %legend('All fibers','','Formin at barbed end')
% figure(4)
% plot(ts,MeanOfEach(5,:),LineSty)
% hold on
% skip=floor(nT/20);
% set(gca,'ColorOrderIndex',iFc)
% errorbar(ts(skip:skip:end),MeanOfEach(5,skip:skip:end),...
%     2*StdOfEach(5,skip:skip:end)/sqrt(nError),'o','MarkerSize',0.1)
% ylabel('\% in formin-bound ')
% % Histogram
% % MeanCountsByTime = cell(nT,1);
% % StdCountsByTime = cell(nT,1);
% % xpl = cell(nT,1);
% % for iT=1:nT
% %     MeanCountsByTime{iT}=mean(CountsByTime{iT});
% %     StdCountsByTime{iT}=std(CountsByTime{iT});
% %     xpl{iT}=(EdgesByTime{iT}(1:end-1)+EdgesByTime{iT}(2:end))*0.5;
% % end
% % figure(5)
% % %subplot(2,3,RunIndex)
% % tspl =[100 300 600 1200 3600 10800]/dt;
% % for iT=tspl
% %     errorbar(xpl{iT},MeanCountsByTime{iT},2*StdCountsByTime{iT}/sqrt(nError),'LineWidth',2.0);
% %     hold on
% %     xlabel('Number of monomers')
% % end
% % %legend(strcat('$t=$',num2str(tspl(1)*dt)),strcat('$t=$',num2str(tspl(2)*dt)),...
% % %    strcat('$t=$',num2str(tspl(3)*dt)))
% % %title(strcat(num2str(ForminConcs(iFc)),' nM Formin'))
% % title(strcat(num2str(ConcProf),' $\mu$M profilin'))
end
end
