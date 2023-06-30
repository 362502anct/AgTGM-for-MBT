function [u,e] = Toeptgm_d1(A0,A1,n,b,tol,choice)

    
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
    
    
    omega1 = 7/8;
    omega2 = 7/12;
    nsize = 1000;
    
    if choice == 'V'
        inix = ToepVcycle_d1(A,b,n,l,nsize,omega1,omega2);
    elseif choice == 'T'
        inix = Toeptwogrid_d1(A,A0inv,A1Zn1,Ap,P,b);   
    end
     
     
    mmax = 10000; % the maximal number of iterations 
    u = inix; 
    
    r = b - A * u;
    e(1) = norm(r)/norm(b);
    fprintf('\n Initial residual = %e',e(1));
    iter = 1;
    t1 = 1.0;
    d = zeros(n*l,1);
    
    while iter < mmax && e(iter) > tol
            
    if choice == 'V'
        z = ToepVcycle_d1(A,r,n,l,nsize,omega1,omega2);
    elseif choice == 'T'
        z = Toeptwogrid_d1(A,A0inv,A1Zn1,Ap,P,r);   
    end
        
    
        u = u + z;
        r = r - A * z;

        iter = iter+1;
         e(iter) = norm(r)/norm(b);
%         fprintf('\n at step %1.0f, relative residual = %e',... 
%             iter-1,e(iter));
        
    end
    fprintf('\n Toep TGM at step %1.0f, relative residual = %e',... 
            iter-1,e(iter));
    

end

