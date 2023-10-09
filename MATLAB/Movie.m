X=load('6000Dt4X_3.txt');
IDs=load('6000Dt4IDs_3.txt');
Nmon=6000;
nT=1001;
Ld=1;
f = figure;
f.Position = [100 100 1000 1000];
NFree = zeros(nT,1);
for iT=1:nT
    inds = (iT-1)*Nmon+1:iT*Nmon;
    Xpl = X(inds,:);
    Xpl = mod(Xpl,[Ld Ld Ld]);
    IDpl = IDs(inds);
    mons = IDpl==-1;
    NFree(iT)=sum(mons);
    scatter3(Xpl(mons,1),Xpl(mons,2),Xpl(mons,3),10,'filled')
    hold on
    maxFib = max(IDs(inds));
    for iFib=0:maxFib
        plInds = IDpl==iFib;
        scatter3(Xpl(plInds,1),Xpl(plInds,2),Xpl(plInds,3),10,'filled')
    end
    %scatter3(Xpl(1,1),Xpl(1,2),Xpl(1,3),50,'k','filled')
    hold off
    view([ -14.2290   28.1680])
    zlim([0 1])
    xlim([0 1])
    ylim([0 1])
    view(3)
    PlotAspect
    movieframes(iT)=getframe(f);
end