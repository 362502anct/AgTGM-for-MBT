function [u,e,iter] = pcgagg_fpd(A,C,R,D1,D2,n,m,b,tol,omega,choice)
    % construct the matrix A
     [l,~] = size(A);
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    Zm1 = diag(ones(m - 1,1),1);
    Zm1 = sparse(Zm1);
    
    A0inv = kron(speye(n),kron(speye(m), sparse(inv(full(A)))));
    A0 = kron(speye(m),A) + kron(Zm1,C) + kron(Zm1',C');
    A1 = kron(speye(m),R) + kron(Zm1',D1);
    AG = kron(speye(n),A0) + kron(Zn1,A1) + kron(Zn1',A1');
    A1Z1 = kron(speye(n),kron(Zm1,C) + kron(Zm1',C')) + kron(Zn1,A1) + kron(Zn1',A1');

    P = ones(1,l)/sqrt(l);
    P = sparse(P);

    AP = P * A * P';
    CP = P * C * P';
    RP = P * R * P';
    D1P = P * D1 * P';
    D2P = P * D2 * P';
    
    G1p = kron(speye(m),AP) + kron(Zm1,CP) + kron(Zm1',CP') ;
    G2p = kron(speye(m),RP) + kron(Zm1', D1P);
    Ap = kron(speye(n),G1p) + kron(Zn1,G2p) + kron(Zn1',G2p');
    Lp = chol(Ap,'lower');
%     
    Pnm = kron(speye(n),kron(speye(m),P'));
%     Ap = Pnm' * A * Pnm;
%     
%     Lp = chol(Ap,'lower');
%     
    omega1 = omega;
    omega2 = omega;
    nsize = 10000;
   
    
    if choice == 'V'
        inix = ToepVcycle_d2(AG,Pnm,b,n,m,l,nsize,omega1,omega2,1);
    elseif choice == 'T'
        inix = Toeptwogrid_d2(AG,A0inv,A1Z1,Ap,Pnm,b,omega1,omega2,Lp);   
    end
     
    mmax = 10000; % the maximal number of iterations 
    u = inix; 
    
    r = b - AG * u;
    e(1) = norm(r)/norm(b);
    fprintf('\n Initial residual = %e',e(1));
    iter = 1;
    t1 = 1.0;
    d = zeros(n*m*l,1);
    
    while iter < mmax && e(iter) > tol
        if choice == 'V'
            z = ToepVcycle_d2(AG,Pnm,r,n,m,l,nsize,omega1,omega2,1);
        elseif choice == 'T'
            z = Toeptwogrid_d2(AG,A0inv,A1Z1,Ap,Pnm,r,omega1,omega2,Lp);   
        end
        
        t1old = t1;
        t1 = z' * r;
        beta = t1/t1old; 
        d = z + beta * d;
        s = AG * d;
        
        suma = d' * s;
        tau = t1/suma;
        u = u + tau * d;
        r = r - tau * s;
        iter = iter + 1;
         e(iter) = norm(r)/norm(b);

        
    end
            fprintf('\n at step %1.0f, relative residual = %e',... 
            iter-1,e(iter));
end

