X=load('X_1.txt');
IDs=load('IDs_1.txt');
Nmon=1000;
nT=100;
f = figure;
%f.Position = [100 100 1000 1000];
NFree = zeros(nT,1);
for iT=1:nT
    inds = (iT-1)*Nmon+1:iT*Nmon;
    Xpl = X(inds,:);
    Xpl = mod(Xpl,[1 1 1]);
    mons = IDs(inds)==-1;
    NFree(iT)=sum(mons);
    bis = IDs(inds) > -1;
    scatter3(Xpl(mons,1),Xpl(mons,2),Xpl(mons,3),'filled')
    hold on
    scatter3(Xpl(bis,1),Xpl(bis,2),Xpl(bis,3),'filled')
    hold off
    view([ -14.2290   28.1680])
    zlim([0 1])
    xlim([0 1])
    ylim([0 1])
    PlotAspect
    movieframes(iT)=getframe(f);
end