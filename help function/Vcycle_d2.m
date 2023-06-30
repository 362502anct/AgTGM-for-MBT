

function x = Vcycle_d2(Acell,Pcell,b,ncell,mcell,lcell,nsize,omega1cell,omega2cell,level,Lp)

    step = 1;
    
    x = zeros(size(b));
    A = Acell{level};
    P = Pcell{level + 1};
    l = lcell(level);
    n = ncell(level);
    m = mcell(level);
    omega1 = omega1cell(level);
    omega2 = omega2cell(level);
    totallev = length(ncell) - 1;
    
    A00 = A(1:l,1:l);
    A0inv = kron(speye(n),kron(speye(m),inv(A00)));
    A1Zn1 = A - kron(speye(n),kron(speye(m),A00));
    
    for i = 1:step
        x1 = A0inv * (b - A1Zn1 * x);
        x = x + omega1 * (x1 - x);
    end
    
    r = b - A * x;
    bP = P' * r;
    
    if level < totallev
       level = level + 1;
       xP = Vcycle_d2(Acell,Pcell,bP,ncell,mcell,lcell,nsize,omega1cell,omega2cell,level,Lp);
    else
       xP = Lp'\(Lp\bP);
    end

    x0 = P * xP;
    x = x + x0;
    
    for i = 1:step
        x1 = A0inv * (b - A1Zn1 * x);
        x = x + omega2 * (x1 - x);
    end
end