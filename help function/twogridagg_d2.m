function [u,e, iter, time] = twogridagg_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,uors,p,omega,choice)

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
    omega2 = omega *2 /3;
    nsize = 10000;
     
     Chol_Ap = chol(Ap,'lower');
        if choice == 'V'
            inix = ToepVcycle_d2(A,P,b,n,m,l,nsize,omega1,omega2,1);
        elseif choice == 'T'
            inix = Toeptwogrid_d2(A,A0inv,A1Z1,Ap,P,b,omega1,omega2,Chol_Ap);   
        end
%      inix = aggsolve_d2(A,Chol_Ap,A0inv,A1Z1,P,b,Ap,n,m,choice,uors);
     
    mmax = 10000; % the maximal number of iterations 
    u = inix; 
    
    r = b - A * u;
    e(1) = norm(r)/norm(b);
    fprintf('\n Initial residual = %e',e(1));
    iter = 1;
    t1 = 1.0;
    d = zeros(n*m*l,1);

    
    while iter < mmax && e(iter) > tol
        if choice == 'V'
            z = ToepVcycle_d2(A,P,r,n,m,l,nsize,omega1,omega2,1);
        elseif choice == 'T'
            z = Toeptwogrid_d2(A,A0inv,A1Z1,Ap,P,r,omega1,omega2,Chol_Ap);   
        end
%         z = aggsolve_d2(A,Chol_Ap,A0inv,A1Z1,P,r,Ap,n,m,choice,uors);
        
        u = u + z;
        r = r - A * z;
        iter = iter + 1;
         e(iter) = norm(r)/norm(b);
%         fprintf('\n at step %1.0f, relative residual = %e',... 
%             iter-1,e(iter));
        
    end
        fprintf('\n agg TGM at step %1.0f, relative residual = %e',... 
            iter-1,e(iter));
    time = cputime - t1;

end