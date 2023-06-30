function x = Toeptwogrid_d1(A,A0inv,A1Zn1,Ap,P,b)

    step = 1;
    omega1 = 7/8;
    omega2 = 7/12;
    
    x = zeros(size(b));
    
    for i = 1:step
        x1 = A0inv * (b - A1Zn1 * x);
        x = x + omega1 * (x1 - x);
    end
    
    r = b - A * x;
    
    bP = P' * r;
    xP = Ap \ bP;
    x0 = P * xP;
    x = x + x0;
    
    for i = 1:step
        x1 = A0inv * (b - A1Zn1 * x);
        x = x + omega2 * (x1 - x);
    end
end