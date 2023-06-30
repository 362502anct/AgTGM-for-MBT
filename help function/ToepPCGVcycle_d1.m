function [u,e,iter] = ToepPCGVcycle_d1(A0,A1,n,b,tol,omega,nsize)

    % construct the matrix
     [l,~] = size(A0);
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    
    A0inv = kron(speye(n), inv(A0));
    A1Zn1 =  kron(Zn1,A1') + kron(Zn1',A1);
    A = kron(speye(n),A0) + A1Zn1;
    
    % construct the prolongation matrix
    P = prolong_d1(n,l);
    Ap = P' * A * P;
    
    
    omega1 = omega;
    omega2 = omega * 2 /3;
    
     % setup
    level = 0;
    ncell = [n];
    lcell = [l];
    np = n * l;
    while np > nsize
        n1 = (ncell(end) - 1)/2;
        ncell = [ncell, n1];
        lcell = [lcell, l];
        np = n1 * l;
        level = level + 1;
    end
    
    Acell = cell(level+1,1);
    Pcell = cell(level+1,1);
    Acell{1} = A;
    Pcell{1} = speye(n *l);
         Pcell{2} = prolong_d1(ncell(1), lcell(1));
    Acell{2} = Pcell{2}' * Acell{1} * Pcell{2};
    omega1 = [omega1, omega];
    omega2 = [omega2, omega * 2 / 3];
    
    if level+1 > 1
        for i = 3:level+1
            Pcell{i} = prolong_d1(ncell(i - 1),  lcell(i - 1));
            Acell{i} = Pcell{i}' * Acell{i - 1} * Pcell{i};
    omega1 = [omega1, omega];
    omega2 = [omega2, omega * 2 / 3];
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
 fprintf('\n PCG Toep Vcycle at step %1.0f, relative residual = %e',... 
            iter-1,e(iter));

end