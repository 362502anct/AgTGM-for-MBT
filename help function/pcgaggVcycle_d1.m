function [u,e,iter] = pcgaggVcycle_d1(A0,A1,n,b,inix,tol,uors,p,omega,nsize)

    if nargin < 9
        omega = 1;
        choice = 'T';
    end

    
    if nargin < 8
        [l,~] = size(A0);
        p = ones(l,1)/sqrt(l);
        omega = 0.7;
        choice = 'T';
    end
    
    if nargin < 7
        uors = 1;
        [l,~] = size(A0);
        p = ones(l,1)/sqrt(l);
        omega = 0.7;
        choice = 'T';
    end
    
    [l,~] = size(A0);
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    
    A0inv = kron(speye(n),omega * inv(A0));
    A1Zn1 = (1 - 1/omega) * kron(speye(n),A0) + kron(Zn1,A1') + kron(Zn1',A1);
    A = (1/omega) * kron(speye(n),A0) + A1Zn1;

     Pl = p;
     Pl = sparse(Pl);
     P = kron(speye(n),Pl);
     
     A0p = Pl' * A0 * Pl;
     A1p = Pl' * A1 * Pl;
     para = [A0p,A1p];
     
     Ap = kron(speye(n),A0p) + kron(Zn1,A1p) + kron(Zn1',A1p');
     
    omega1 = omega;
    omega2 = omega;
    
    % setup
    level = 1;
    ncell = [n, n];
    lcell = [l, 1];
    np = n;
    while np > nsize
        n1 = (ncell(end) - 1)/2;
        ncell = [ncell, n1];
        lcell = [lcell, 1];
        np = n1;
        level = level + 1;
    end
    
    Acell = cell(level+1,1);
    Pcell = cell(level+1,1);
    Acell{1} = A;
    Acell{2} = Ap;
    Pcell{1} = speye(n *l);
    Pcell{2} = P;
    omega1 = [omega1, 7/8];
    omega2 = [omega2, 7/12];
    
    if level+1 > 1
        for i = 3:level+1
            Pcell{i} = prolong_d1(ncell(i - 1),  lcell(i - 1));
            Acell{i} = Pcell{i}' * Acell{i - 1} * Pcell{i};
                omega1 = [omega1, 7/8];
                omega2 = [omega2, 7/12];
        end
    end
    
     Chol_Ap = chol(Acell{end},'lower');
     
     finelevel = 1;
     inix = Vcycle_d1(Acell,Pcell,b,ncell,lcell,nsize,omega1,omega2,finelevel,Chol_Ap);
     
    mmax = 10000; % the maximal number of iterations 
    u = inix; 
    
    r = b - A * u;
    e(1) = norm(r)/norm(b);
    fprintf('\n Initial residual = %e',e(1));
    iter = 1;
    t1 = 1.0;
    d = zeros(n*l,1);
    
    while iter < mmax && e(iter) > tol
        z = Vcycle_d1(Acell,Pcell,r,ncell,lcell,nsize,omega1,omega2,finelevel,Chol_Ap);
        
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

end

