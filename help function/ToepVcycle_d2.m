function x = ToepVcycle_d2(A,P,b,n,m,l,nsize,omega1,omega2,level)

    if nargin < 10
        level = 0;
    end

    step = 1;
    
    x = zeros(size(b));
    A00 = A(1:l,1:l);
    A0inv = kron(speye(n),kron(speye(m),inv(A00)));
    A1Zn1 = A - kron(speye(n),kron(speye(m),A00));
    
    for i = 1:step
        x1 = A0inv * (b - A1Zn1 * x);
        x = x + omega1 * (x1 - x);
    end
    
    r = b - A * x;
    bP = P' * r;
    Ap = P' * A * P;
    if level <= 0
        n1 = (n - 1)/2;
        m1 = (m - 1)/2;
    else
        n1 = n;
        m1 = n;
    end
    [np,~] = size(Ap);
    
    if np > nsize
       l = np / (n1*m1);
       P1 = prolong_d2(n1,m1,l);
       xP = ToepVcycle_d2(Ap,P1,bP,n1,m1,l,nsize,omega1,omega2); 
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