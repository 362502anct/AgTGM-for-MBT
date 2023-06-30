function [u,e,iter] = Toeptgm_fpd(A,C,R,D1,D2,n,m,b,tol,omega,choice)

    % construct the matrix A
     [l,~] = size(A);
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    Zm1 = diag(ones(m - 1,1),1);
    Zm1 = sparse(Zm1);
    
    A0inv = kron(speye(n),kron(speye(m), inv(A)));
    A0 = kron(speye(m),A) + kron(Zm1,C) + kron(Zm1',C');
    A1 = kron(speye(m),R) + kron(Zm1',D1) + kron(Zm1,D2);
    A = kron(speye(n),A0) + kron(Zn1,A1) + kron(Zn1',A1');
    A1Z1 = kron(speye(n),kron(Zm1,C) + kron(Zm1',C')) + kron(Zn1,A1) + kron(Zn1',A1');

    P = prolong_d2(n,m,l);
    Ap = P' * A * P;
    
    Lp = chol(Ap,'lower');
    
    omega1 = omega;
    omega2 = omega * 2 / 3;
    nsize = 10000;
    
    
    if choice == 'V'
        inix = ToepVcycle_d2(A,P,b,n,m,l,nsize,omega1,omega2);
    elseif choice == 'T'
        inix = Toeptwogrid_d2(A,A0inv,A1Z1,Ap,P,b,omega1,omega2,Lp);   
    end
     
    mmax = 10000; % the maximal number of iterations 
    u = inix; 
    
    r = b - A * u;
    e(1) = norm(r)/norm(b);
    fprintf('\n Initial residual = %e',e(1));
    iter = 1;
    
    while iter < mmax && e(iter) > tol
        if choice == 'V'
            z = ToepVcycle_d2(A,P,r,n,m,l,nsize,omega1,omega2);
        elseif choice == 'T'
            z = Toeptwogrid_d2(A,A0inv,A1Z1,Ap,P,r,omega1,omega2,Lp);   
        end
       
        u = u + z;
        r = r - A * z;
        iter = iter + 1;
         e(iter) = norm(r)/norm(b);
%         fprintf('\n at step %1.0f, relative residual = %e',... 
%             iter-1,e(iter));
        
    end

        fprintf('\n Toep TGM at step %1.0f, relative residual = %e',... 
            iter-1,e(iter));

end