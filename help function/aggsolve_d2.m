function x = aggsolve_d2(A,Chol,A0inv,A1Zn1,P,b,Ap,n,m,choice,uors)
    if nargin < 11
       uors = 1; 
    end
    
    step = 1;
    x = zeros(size(b));
    for i = 1:step       
        B = (b - A1Zn1 * x);
        x1 = A0inv * B;
        x = x1;
    end
    
    r = b - A * x;
    
    bP = P' * r;
    
    if choice == 'T'
        xP = Chol'\(Chol\bP);
    elseif choice == 'V'
        omega1 = 7/8;
        omega2 = 7/12;
        nsize = 1000;
        xP = ToepVcycle_d2(Ap,bP,n,m,1,nsize,omega1,omega2);
    end
    
    x0 = P * xP;
    x = x + x0;
    
    if uors == 1
        for i = 1:step
            B = (b - A1Zn1 * x);
            x1 = A0inv * B;
            x = x1;
        end
    end
end