function [u,e, iter, time] = ToepVcyclerun_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,omega,nsize)

    t1 = cputime;
  % construct the matrix A
     [l,~] = size(A00);
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    Zm1 = diag(ones(m - 1,1),1);
    Zm1 = sparse(Zm1);
    
    A0inv = kron(speye(n),kron(speye(m), inv(A00)));
    A0 = kron(speye(m),A00) + kron(Zm1,A10) + kron(Zm1',A_10);
    A1 = kron(speye(m),A01);
    A = kron(speye(n),A0) + kron(Zn1,A1) + kron(Zn1',A1');
    A1Z1 = A - kron(speye(n),kron(speye(m),A00));

    P = prolong_d2(n,m,l);
    Ap = P' * A * P;
    Lp = chol(Ap,'lower');
    
    omega1 = omega;
    omega2 = omega * 2 / 3;

     % setup
    level = 0;
    ncell = [n];
    mcell = [m];
    lcell = [l];
    np = n * m * l;
    while np > nsize
        n1 = (ncell(end) - 1)/2;
        m1 = (mcell(end) - 1)/2;
        ncell = [ncell, n1];
        mcell = [mcell, m1];
        lcell = [lcell, l];
        np = n1 * m1 * l;
        level = level + 1;
    end
    
    Acell = cell(level+1,1);
    Pcell = cell(level+1,1);
    Acell{1} = A;
    Pcell{1} = speye(n *m *l);
    Pcell{2} = prolong_d2(ncell(1), mcell(1), lcell(1));
    Acell{2} = Pcell{2}' * Acell{1} * Pcell{2};
    omega1 = [omega1, omega];
    omega2 = [omega2, omega * 2 / 3];
    
    if level+1 > 2
        for i = 3:level+1
            Pcell{i} = prolong_d2(ncell(i - 1), mcell(i - 1), lcell(i - 1));
            Acell{i} = Pcell{i}' * Acell{i - 1} * Pcell{i};
            omega1 = [omega1, omega];
            omega2 = [omega2, omega * 2 / 3];
        end
    end
    
     Chol_Ap = chol(Acell{end},'lower');
     
     finelevel = 1;
     inix = Vcycle_d2(Acell,Pcell,b,ncell,mcell,lcell,nsize,omega1,omega2,finelevel,Chol_Ap);
     
    mmax = 10000; % the maximal number of iterations 
    u = inix; 
    
    r = b - A * u;
    e(1) = norm(r)/norm(b);
%     fprintf('\n Initial residual = %e',e(1));
    iter = 1;
    
    while iter < mmax && e(iter) > tol
        z = Vcycle_d2(Acell,Pcell,r,ncell,mcell,lcell,nsize,omega1,omega2,finelevel,Chol_Ap);

        u = u + z;
        r = r - A * z;
        iter = iter + 1;
         e(iter) = norm(r)/norm(b);
%         fprintf('\n at step %1.0f, relative residual = %e',... 
%             iter-1,e(iter));
        
    end
        fprintf('\n Toep Vcycle at step %1.0f, relative residual = %e',... 
            iter-1,e(iter));
    time = cputime - t1;

end