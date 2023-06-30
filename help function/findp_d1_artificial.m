function [lambda,v,x,p] = findp_d1_artifical(D,A,R)
    
    L = chol(D,'lower');
    Linv = inv(L);
    LRL = Linv * R * Linv';
    LAL = Linv * A * Linv';
    [l,~] = size(A);
    
    x0 = pi/4;
    alpha = 1;

    F = @(x) LAL + LRL * exp(-1i * x) + LRL' * exp(1i * x) ;
 
    F0 = F(x0);
    [v02,lambda02] = eigs(F0,2,'sm');
    v0 = v02(:,2);
    Lp = v02(:,1);
    lambda0 = lambda02(2,2);
 
    dFx = @(x) 1i * LRL' * exp(1i * x)  - 1i * LRL * exp(-1i * x) ;

    lambda = 0;
    eps = 1e-15;
    v = v0;

while abs(lambda - lambda0)/abs(lambda) > eps
    lambda0 = lambda;
    v0 = v;
    
    Ds = eye(l) - Lp * Lp' / (Lp' * Lp);
    
    dFx1 = dFx(x0);
    dlambdax = v0' * dFx1 * v0 / (v0' * Ds * v0);
   
    x = x0 - alpha * dlambdax;
   
    F1 = F(x);
    [v2,lambda2] = eigs(F1,2,'sm');
    lambda = lambda2(2,2);
    v = v2(:,2);
    Lp = v2(:,1);
    
    x0 = x;
   
end

    p = L\v2(:,1);
    
end
