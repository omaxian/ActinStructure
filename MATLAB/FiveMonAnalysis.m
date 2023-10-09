% Parameters
Nmon=6000;
a = 4e-3;
Raa = 2*a;
dplus = 10;
on_m = 0.005420;
off_m = 0.041;
on_p = 1006.570112+8981.702539;
off_p = 2.2;
mu = 0.01;
kbT = 4.1e-3;
R=on_m*Raa^2;
D=2*kbT/(6*pi*mu*a);
r=R/D % << 1 is reaction limited, >> 1 is diffusion limited
% Determine effective on rate
uMInvToMicron3 = 1.0e15/(6.022e17);
kplus_d = 3.5e-6*uMInvToMicron3;%1/2*4*pi/3*Raa^3*on_m;
kminus_d = 0.041;%off_m;
kplus_p = 12.9*uMInvToMicron3;%4*pi/3*Raa^3*on_p;
kminus_p = 2.2;%off_p;
% Steady state equations
myfun = @(c) EqFcn(c,kplus_p,kminus_p,kplus_d,kminus_d,Nmon);
nMax=Nmon;
c0=ones(nMax,1);
x = fsolve(myfun,c0);
alpha=(kplus_p*x(1)/kminus_p);
%nMons=[nMons;Nmon];
%alphas=[alphas;alpha];

% lb = zeros(nMax,1);
% myfunC = @(c) EqFcnC(c,kplus_p,kminus_p,kplus_d,kminus_d,Nmon);
% opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
% x3 = fmincon(@(x)0,c0,[],[],[],[],lb,[],myfunC,opts);

% for pmon=1:4
% imon=pmon;
% if (imon>1)
%     imon=pmon+1;
% end
% MeanOfAll=MeanOfAll_3;
% StdOfAll=StdOfAll_3;
% ind=2;
% ts=(0:nTs-1)/(nTs-1)*2000;
% subplot(2,2,pmon)
% set(gca,'ColorOrderIndex',ind)
% plot(ts,MeanOfAll(imon,:))
% hold on
% set(gca,'ColorOrderIndex',ind)
% skip = 50;
% start = 25;
% errorbar(ts(start:skip:end),MeanOfAll(imon,start:skip:end),...
%     2*StdOfAll(imon,start:skip:end)/sqrt(nTrial),'o','MarkerSize',0.1,'LineWidth',2.0)
% plot(xlim,[x(imon) x(imon)],':k')
% title(strcat(num2str(imon),' monomers'))
% end
% % Confidence interval
% mean(MeanOfAll(:,(end-1)/2+1:end)')
% StdOfAll(:,(end-1)/2+1:end)'

function [val1,val] =  EqFcnC(c,kplus_p,kminus_p,kplus_d,kminus_d,Nmon)
    val1=0*c;
    nMax=length(c);
    val = zeros(nMax,1);
    val(1) = sum((1:nMax)'.*c)-Nmon;
    val(2) = kminus_d*c(2)-kplus_d.*c(1).^2;   % Dimers
    for iP=2:nMax-1
        val(iP+1) = kminus_p*c(iP+1)-kplus_p*c(iP).*c(1);
    end
end

function val = EqFcn(c,kplus_p,kminus_p,kplus_d,kminus_d,Nmon)
    nMax=length(c);
    val = zeros(nMax,1);
    val(1) = sum((1:nMax)'.*c)-Nmon;
    val(2) = kminus_d*c(2)-kplus_d.*c(1).^2;   % Dimers
    for iP=2:nMax-1
        val(iP+1) = kminus_p*c(iP+1)-kplus_p*c(iP).*c(1);
    end
end