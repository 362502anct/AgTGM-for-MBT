function [u,e,iter] = Vcycleagg_fpd(A,C,R,D1,D2,n,m,b,tol,omega,omega11,nsize)
    % construct the matrix A
     [l,~] = size(A);
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    Zm1 = diag(ones(m - 1,1),1);
    Zm1 = sparse(Zm1);
    
%     A0inv = kron(speye(n),kron(speye(m), sparse(inv(full(A)))));
    A0 = kron(speye(m),A) + kron(Zm1,C) + kron(Zm1',C');
    A1 = kron(speye(m),R) + kron(Zm1',D1);
    AG = kron(speye(n),A0) + kron(Zn1,A1) + kron(Zn1',A1');
%     A1Z1 = kron(speye(n),kron(Zm1,C) + kron(Zm1',C')) + kron(Zn1,A1) + kron(Zn1',A1');

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
%     
    Pnm = kron(speye(n),kron(speye(m),P'));
%     Ap = Pnm' * A * Pnm;
%     
%     Lp = chol(Ap,'lower');
%     
    omega1 = omega;
    omega2 = omega * 2 / 3;
   
    % setup
    level = 1;
    ncell = [n, n];
    mcell = [m, m];
    lcell = [l, 1];
    np = n * m;
    while np > nsize
        n1 = (ncell(end) - 1)/2;
        m1 = (mcell(end) - 1)/2;
        ncell = [ncell, n1];
        mcell = [mcell, m1];
        lcell = [lcell, 1];
        np = n1 * m1;
        level = level + 1;
    end
    
    Acell = cell(level+1,1);
    Pcell = cell(level+1,1);
    Acell{1} = AG;
    Acell{2} = Ap;
    Pcell{1} = speye(n *m *l);
    Pcell{2} = Pnm;
%     omega1 = [omega1, 7/8];
%     omega2 = [omega2, 7/12];
    omega1 = [omega1, omega11];
    omega2 = [omega2, omega11 * 2 / 3];
    
    
    if level+1 > 2
        for i = 3:level+1
            Pcell{i} = prolong_d2(ncell(i - 1), mcell(i - 1), lcell(i - 1));
            Acell{i} = Pcell{i}' * Acell{i - 1} * Pcell{i};
%         omega1 = [omega1, 7/8];
%         omega2 = [omega2, 7/12];
            omega1 = [omega1, omega11];
            omega2 = [omega2, omega11 * 2 / 3];

        end
    end
    
     Chol_Ap = chol(Acell{end},'lower');
     
     finelevel = 1;
     inix = Vcycle_d2(Acell,Pcell,b,ncell,mcell,lcell,nsize,omega1,omega2,finelevel,Chol_Ap);
    
     
    mmax = 10000; % the maximal number of iterations 
    u = inix; 
    
    r = b - AG * u;
    e(1) = norm(r)/norm(b);
    fprintf('\n Initial residual = %e',e(1));
    iter = 1;
    
    while iter < mmax && e(iter) > tol
        z = Vcycle_d2(Acell,Pcell,r,ncell,mcell,lcell,nsize,omega1,omega2,finelevel,Chol_Ap);
        
        u = u + z;
        r = r - AG * z;
        iter = iter + 1;
         e(iter) = norm(r)/norm(b);
%         fprintf('\n at step %1.0f, relative residual = %e',... 
%             iter-1,e(iter));
        
    end
    fprintf('\n agg Vcycle at step %1.0f, relative residual = %e',... 
            iter-1,e(iter));

end