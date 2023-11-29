function X = BranchedFiberPoints(X0,nPerFib,AttachedTo,AttachmentPts,taus,a)
    spacing = 2*a;
    % Mother filament
    X  = X0 + taus(1,:).*(0:nPerFib(1)-1)'*spacing;
    % Build the branches
    nFil = length(nPerFib);
    CumN = [0;cumsum(nPerFib)];
    for iB=2:nFil
        Mother = AttachedTo(iB);
        AttachmentNode = CumN(Mother)+AttachmentPts(iB);
        Xa = X(AttachmentNode,:);
        NewX = Xa+ taus(iB,:).*(1:nPerFib(iB))'*spacing;
        X = [X;NewX];
    end
end
    