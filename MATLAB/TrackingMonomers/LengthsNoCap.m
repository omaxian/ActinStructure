% Simulation data
Nmon=8000;
nTs = 1001;
nTrial=2;
nSeed=5;
MeanNumMonomers = zeros(nTrial,nTs);
MeanFiberLengths = zeros(10,nTs,nTrial);
maxFibs=0;
for iTrial=1:nTrial
    NumMonomers = zeros(nTs,nSeed);
    FiberLengths = zeros(10,nTs,nSeed);
    AllLengthsTrial = [];
    for iseed=1:nSeed
        thisSeed = (iTrial-1)*nSeed+iseed;
        IDs = load(strcat('8000Ld2Dt4IDs_',num2str(thisSeed),'.txt'));
        % Find all the monomers
        for iT=1:nTs
            inds = (iT-1)*Nmon+1:iT*Nmon;
            IDThis = IDs(inds);
            uIDs = unique(IDThis);  
            uIDs = sort(uIDs,'ascend');
            if (length(uIDs)-1 > maxFibs)
                maxFibs=length(uIDs)-1;
            end
            nMonsbyID = zeros(length(uIDs)-1,1);
            for iiD=1:length(uIDs)
                IDnum=uIDs(iiD);
                if (IDnum==-1)
                    nMons = sum(IDThis==IDnum);
                    NumMonomers(iT,iseed)=nMons;
                else
                    nMons = sum(IDThis==IDnum);
                    nMonsbyID(iiD-1) = nMons;
                end
            end
            FiberLengths(1:length(uIDs)-1,iT,iseed)=sort(nMonsbyID,'descend');
            if (iT > (nTs-1)/2) % Steady state length distribution
                AllLengthsTrial=[AllLengthsTrial;nMonsbyID];
            end
        end
    end
    MeanNumMonomers(iTrial,:)=mean(NumMonomers,2);
    MeanFiberLengths(:,:,iTrial)=mean(FiberLengths,3);
    AllLengthsByTrial{iTrial}=AllLengthsTrial;
end
MeanFiberLengths = MeanFiberLengths(1:maxFibs,:,:);
% Plot the mean and standard deviation and compare to theory
MeanOfAll = zeros(maxFibs,nTs);
StdOfAll = zeros(maxFibs,nTs);
for iFib=1:maxFibs
    MeanOfAll(iFib,:)=mean(MeanFiberLengths(iFib,:,:),3);
    arr = reshape(MeanFiberLengths(iFib,:,:),nTs,nTrial)';
    StdOfAll(iFib,:)=std(arr);
end

meanMon = mean(MeanNumMonomers);
stdMon = std(MeanNumMonomers);
%tiledlayout(1,maxFibs+1, 'Padding', 'none', 'TileSpacing', 'compact');
%nexttile
subplot(2,4,1)
% Plot monomers
tf=2000;
ts=(0:nTs-1)/(nTs-1)*tf;
plot(ts,meanMon)
hold on
skip = 50;
start = 25;
set(gca,'ColorOrderIndex',1)
errorbar(ts(start:skip:end),meanMon(start:skip:end),...
    2*stdMon(start:skip:end)/sqrt(nTrial),'o','MarkerSize',0.1,'LineWidth',2.0)
%plot(xlim,xBar*[1 1],':k')
title('Free monomers')
xlabel('$t$ (s)')
ylabel('Number of monomers')

% Plot fibers
for pmon=1:maxFibs
subplot(2,4,pmon+1)
%nexttile
plot(ts,MeanOfAll(pmon,:))
hold on
set(gca,'ColorOrderIndex',1)
skip = 50;
start = 25;
errorbar(ts(start:skip:end),MeanOfAll(pmon,start:skip:end),...
    2*StdOfAll(pmon,start:skip:end)/sqrt(nTrial),'o','MarkerSize',0.1,'LineWidth',2.0)
title(strcat('Fiber ',num2str(pmon)))
xlabel('$t$ (s)')
end