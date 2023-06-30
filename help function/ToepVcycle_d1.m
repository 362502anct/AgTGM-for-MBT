function x = ToepVcycle_d1(A,b,n,l,nsize,omega1,omega2)

    step = 1;
    
    P = prolong_d1(n,l);
    x = zeros(size(b));
    A00 = A(1:l,1:l);
    A0inv = kron(speye(n),inv(A00));
    A1Zn1 = A - kron(speye(n),A00);
    
    for i = 1:step
        x1 = A0inv * (b - A1Zn1 * x);
        x = x + omega1 * (x1 - x);
    end
    
    r = b - A * x;
    bP = P' * r;
    Ap = P' * A * P;
    [np,~] = size(Ap);
    
    if np > nsize
       xP = ToepVcycle_d1(Ap,bP,np/l,l,nsize,omega1,omega2); 
    else
       xP = Ap \ bP;
    end

    x0 = P * xP;
    x = x + x0;
    
    for i = 1:step
        x1 = A0inv * (b - A1Zn1 * x);
        x = x + omega2 * (x1 - x);
    end
end