function [lambda,v,x] = grad_mf_d1(A,R,Rt,B,p)
if nargin < 5
    Rt = R';
end

x0 = pi;
alpha = 0.05;

F = @(x) A + R * exp(-1i * x) + Rt * exp(1i * x) ;
 
 F0 = F(x0);
 if p
      [v0,lambda0] = eigs(F0,B,1,'sm');
 else 
      [v0,lambda0] = eigs(F0,B,1);
 end

 
 dFx = @(x) 1i * Rt * exp(1i * x)  - 1i * R * exp(-1i * x) ;

lambda = 0;
eps = 1e-12;
v = v0;

while abs(lambda0 - lambda) > eps
    lambda0 = lambda;
    v0 = v;
    
   dFx1 = dFx(x0);
   dlambdax = v0' * dFx1 * v0/(v0' * B * v0);
   
   x = x0 - alpha * dlambdax;
   
    F1 = F(x);
     if p
          [v,lambda] = eigs(F1,B,1,'sm');
     else 
          [v,lambda] = eigs(F1,B,1);
     end
    
    x0 = x;
   
end
end

