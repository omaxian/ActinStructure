function K = RigidRTMatrix(X)
    [N,~]=size(X);
    TauCross = zeros(3*N,3);
    I = zeros(3*N,3);
    COM = X(1,:);%mean(X);
    for iR=1:N
        inds = (iR-1)*3+1:iR*3;
        I(inds,:)=eye(3);
        Disp = X(iR,:)-COM;
        TauCross(inds,1:3)=-CPMatrix(Disp);
    end
    K = TauCross;%[TauCross I];
end