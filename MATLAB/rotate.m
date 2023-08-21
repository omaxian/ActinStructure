% Rodriguez rotation formula
function rotated_x=rotate(x,Omega)
    [nX,~] = size(x);
    nOm = norm(Omega);
    Omhat = Omega/norm(Omega);
    Px = Omhat.*sum(Omhat.*x,2);
    Crosses = zeros(nX,3);
    for iP=1:nX
        Crosses(iP,:)=cross(Omhat,x(iP,:));
    end
    rotated_x = Px+cos(nOm)*(x-Px)+sin(nOm)*Crosses;
    if (nOm < 1e-10)
        rotated_x = x;
    end
end
