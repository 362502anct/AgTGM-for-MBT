function [u,e,iter, time] = pcgaggVcycle_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,uors,p,omega,omega11,nsize)

    t1 = cputime;
    if nargin < 11
        [l,~] = size(A00);
        p = ones(l,1)/sqrt(l);
    end
    
    if nargin < 12
        omega = 1;
    end

    if nargin < 10
        nors = 1;
        [l,~] = size(A00);
        p = ones(l,1)/sqrt(l);
    end
        
     [l,~] = size(A00);
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    Zm1 = diag(ones(m - 1,1),1);
    Zm1 = sparse(Zm1);
    
    A0inv = kron(speye(n),kron(speye(m),omega * inv(A00)));
    A0 = kron(speye(m),A00) + kron(Zm1,A10) + kron(Zm1',A_10);
    A1 = kron(speye(m),A01);
    A = kron(speye(n),A0) + kron(Zn1,A1) + kron(Zn1',A1');
    A1Z1 = A - (1/omega) * kron(speye(n),kron(speye(m),A00));

     Pl = p;
     Pl = sparse(Pl);
     P = kron(speye(n),kron(speye(m),Pl));
     
     A00p = Pl' * A00 * Pl;
     A01p = Pl' * A01 * Pl;
     A10p = Pl' * A10 * Pl;
     A0_1p = Pl' * A0_1 * Pl;
     A_10p = Pl' * A_10 * Pl;


     A0p = kron(speye(m),A00p) + kron(Zm1,A10p) + kron(Zm1',A_10p);
     A1p = kron(speye(m),A01p);
     Ap = kron(speye(n),A0p) + kron(Zn1,A1p) + kron(Zn1',A1p');
     
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
    Acell{1} = A;
    Acell{2} = Ap;
    Pcell{1} = speye(n *m *l);
    Pcell{2} = P;
    omega1 = [omega1, omega11];
    omega2 = [omega2, omega11 * 2 / 3];
    
    if level+1 > 2
        for i = 3:level+1
            Pcell{i} = prolong_d2(ncell(i - 1), mcell(i - 1), lcell(i - 1));
            Acell{i} = Pcell{i}' * Acell{i - 1} * Pcell{i};
            omega1 = [omega1, omega11];
            omega2 = [omega2, omega11 * 2 / 3];
        end
    end
    
     Chol_Ap = chol(Acell{end},'lower');
     
     finelevel = 1;
     inix = Vcycle_d2(Acell,Pcell,b,ncell,mcell,lcell,nsize,omega1,omega2,finelevel,Chol_Ap);
        
    mmax = 10000; % the maximal number of iterations 
    u = inix; 
    
    r = b - A * u;
    e(1) = norm(r)/norm(b);
    fprintf('\n Initial residual = %e',e(1));
    iter = 1;
    t1 = 1.0;
    d = zeros(n*m*l,1);
    
    while iter < mmax && e(iter) > tol

        z = Vcycle_d2(Acell,Pcell,r,ncell,mcell,lcell,nsize,omega1,omega2,finelevel,Chol_Ap);

         t1old = t1;
        t1 = z' * r;
        beta = t1/t1old; 
        d = z + beta * d;
        s = A * d;
        
        suma = d' * s;
        tau = t1/suma;
        u = u + tau * d;
        r = r - tau * s;
        iter = iter + 1;
         e(iter) = norm(r)/norm(b);
%         fprintf('\n at step %1.0f, relative residual = %e',... 
%             iter-1,e(iter));
        
    end
        fprintf('\n PCG agg Vcycle at step %1.0f, relative residual = %e',... 
            iter-1,e(iter));
        time = cputime - t1;

end