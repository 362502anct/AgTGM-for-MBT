function [lambda,v,x,y,p] = findp_d2(A,R,C)
    [l,~] = size(A);

    L = chol(A,'lower');
    Linv = inv(L);
    LRL = Linv * R * Linv';
    LCL = Linv * C * Linv';
    [l,~] = size(A);
    
    x0 = pi/4;
    y0 = pi/4;
    alpha = 1;

    F = @(x,y) eye(l) + LRL * exp(-1i * x) + LRL' * exp(1i * x) +  ...
        LCL * exp(-1i * y) + LCL' * exp(1i * y);
 
    F0 = F(x0,y0);
    [v02,lambda02] = eigs(F0,2,'sm');
    v0 = v02(:,2);
    Lp = v02(:,1);
    lambda0 = lambda02(2,2);
 
    dFx = @(x) 1i * LRL' * exp(1i * x)  - 1i * LRL * exp(-1i * x);
    dFy = @(y) 1i * LCL' * exp(1i * y)  - 1i * LCL * exp(-1i * y);

    lambda = 0;
    eps = 1e-13;
    v = v0;

while abs(lambda - lambda0)/abs(lambda) > eps
    lambda0 = lambda;
    v0 = v;
    
    D = eye(l) - Lp * Lp' / (Lp' * Lp);
    
    dFx1 = dFx(x0);
    dFy1 = dFy(y0);
    dlambdax = v0' * dFx1 * v0 / (v0' * D * v0);
    dlambday = v0' * dFy1 * v0 / (v0' * D * v0);
   
    x = x0 - alpha * dlambdax;
    y = y0 - alpha * dlambday;
   
    F1 = F(x,y);
    [v2,lambda2] = eigs(F1,2,'sm');
    lambda = lambda2(2,2);
    v = v2(:,2);
    Lp = v2(:,1);
    
    x0 = x;
    y0 = y;
end

    p = L\v2(:,1);
    
end