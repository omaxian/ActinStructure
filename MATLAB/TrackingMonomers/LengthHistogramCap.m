% Simulation data
Nmon=1000;
nTs = 2001;
nTrial=2;
nSeed=5;
MeanNumOfEach = zeros(5,nTs,nTrial);
for iTrial=1:nTrial
    NumOfEach = zeros(5,nTs,nSeed);
    TrialNum = zeros(nSeed,nTs);
    for iseed=1:nSeed
        thisSeed = (iTrial-1)*nSeed+iseed;
        IDs = load(strcat('Max5Dt4IDs_',num2str(thisSeed),'.txt'));
        % Find all the monomers
        for iT=1:nTs
            inds = (iT-1)*Nmon+1:iT*Nmon;
            IDThis = IDs(inds);
            uIDs = unique(IDThis);       
            for iiD=1:length(uIDs)
                IDnum=uIDs(iiD);
                if (IDnum==-1)
                    nMons = sum(IDThis==IDnum);
                    NumOfEach(1,iT,iseed)=nMons;
                else
                    nMons = sum(IDThis==IDnum);
                    NumOfEach(nMons,iT,iseed)=NumOfEach(nMons,iT,iseed)+1;
                end
            end
        end
    end
    MeanNumOfEach(:,:,iTrial)=mean(NumOfEach,3);
end
% Plot the mean and standard deviation and compare to theory
MeanOfAll_3 = zeros(5,nTs);
StdOfAll_3 = zeros(5,nTs);
for nMon=1:5
    MeanOfAll_3(nMon,:)=mean(MeanNumOfEach(nMon,:,:),3);
    arr = reshape(MeanNumOfEach(nMon,:,:),nTs,nTrial)';
    StdOfAll_3(nMon,:)=std(arr);
end