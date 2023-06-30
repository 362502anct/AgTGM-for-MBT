function x = Toeptwogrid_d2(A,A0inv,A1Zn1,Ap,P,b,omega1,omega2,Lp)

    step = 1;
    if nargin < 8
        omega1 = 7/8;
        omega2 = 7/12;
    end
    
    x = zeros(size(b));
    
    for i = 1:step
        B = (b - A1Zn1 * x);
        x1 = A0inv * B;
%         x1 = reshape(X1,N,1);
        x = x + omega1 * (x1 - x);
    end
    
    r = b - A * x;
    
    bP = P' * r;
    if(nargin < 9)   
        xP = Ap \ bP;
    else
        xP = Lp'\(Lp\bP);
    end
    x2 = P * xP;
    x = x + x2;
    
    for i = 1:step
%         x3 = A0inv * (b - A1Zn1 * x);
        B = (b - A1Zn1 * x);
        x3 = A0inv * B;
%         x3 = reshape(X3,N,1);
        x = x + omega2 * (x3 - x);
    end
end