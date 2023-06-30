function [lambda,v,x,y] = grad_mf_d2(A,R,C,B,p)
%% test the eigenvalue by gradient descent method
x0 = pi/4;
y0 = pi/4;
alpha = 1;

F = @(x,y) A + R * exp(-1i * x) + R' * exp(1i * x) + ...
     C * exp(-1i * y) + C' * exp(1i * y);
 
 F0 = F(x0,y0);
 if p
      [v0,lambda0] = eigs(F0,B,1,'sm');
 else 
      [v0,lambda0] = eigs(F0,B,1);
 end
 
 dFx = @(x,y) 1i * R' * exp(1i * x)  - 1i * R * exp(-1i * x) ;
 dFy = @(x,y) 1i * C' * exp(1i * y)  - 1i * C * exp(-1i * y) ;

lambda = 0;
eps = 1e-9;
v = v0;

while abs(lambda - lambda0) > eps
    lambda0 = lambda;
    v0 = v;
    
   dFx1 = dFx(x0,y0);
   dFy1 = dFy(x0,y0);
   
   dlambdax = v0' * dFx1 * v0/(v0' * B * v0);
   dlambday = v0' * dFy1 * v0/(v0' * B * v0);
   
   x = x0 - alpha * dlambdax;
   y = y0 - alpha * dlambday;
   
    F1 = F(x,y);
     if p
          [v,lambda] = eigs(F1,B,1,'sm');
     else 
          [v,lambda] = eigs(F1,B,1);
     end
    
    x0 = x;
    y0 = y;
   
end

end

