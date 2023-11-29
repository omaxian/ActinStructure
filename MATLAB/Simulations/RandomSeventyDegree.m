function tau2 = RandomSeventyDegree(tau)
    % Find a random vector an angle of 70 degrees 
    Seventy = 70*pi/180;
    v = randn(1,3);
    v = v - dot(v,tau)*tau;
    v = v/norm(v);
    tau2 = sin(Seventy)*v + cos(Seventy)*tau;
end