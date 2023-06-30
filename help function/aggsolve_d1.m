function x = aggsolve_d1(A,para,Ap,A0inv,A1Zn1,A0inv_post, A1Zn1_post,P,b,choice,uors)

    if nargin < 11
       uors = 1; 
    end
    
    step = 1;
    
    x = zeros(size(b));
    for i = 1:step
        x1 = A0inv * (b - A1Zn1 * x);
        x = x1;
    end
    
    r = b - A * x;
    
    bP = P' * r;
    
    if choice == 'T'
        xP = Ap\bP;
     %   xP = dst_d1(bP,para);
    elseif choice == 'V'
        omega1 = 7/8;
        omega2 = 7/12;
        nsize = 1000;
        [n,~] = size(Ap);
        xP = ToepVcycle_d1(Ap,bP,n,1,nsize,omega1,omega2);
    end
    x0 = P * xP;
    x = x + x0;
    
    if uors == 1
        for i = 1:step
            b_smooth = b - A1Zn1_post * x;
            x1 = A0inv_post * b_smooth;
            x = x1;
        end
    end
end

