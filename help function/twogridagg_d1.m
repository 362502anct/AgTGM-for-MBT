function [u,e] = twogridagg_d1(A0,A1,n,b,inix,tol,uors,p,omega)

    if nargin < 8
            [l,~] = size(A0);
            p = ones(l,1)/sqrt(l);
            omega = 0.7;
    end
    
    choice = 'T';
    [l,~] = size(A0);
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    
    A0inv = kron(speye(n),omega * inv(A0));
    A1Zn1 = (1 - 1/omega) * kron(speye(n),A0) + kron(Zn1,A1') + kron(Zn1',A1);
    omega_post =  omega;
    A0inv_post = kron(speye(n),omega_post * inv(A0));
    A1Zn1_post = (1 - 1/omega_post) * kron(speye(n),A0) + kron(Zn1,A1') + kron(Zn1',A1);
    A = (1/omega) * kron(speye(n),A0) + A1Zn1;

     Pl = p;
     Pl = sparse(Pl);
     P = kron(speye(n),Pl);
     
     A0p = Pl' * A0 * Pl;
     A1p = Pl' * A1 * Pl;
     para = [A0p,A1p];
     
     Ap = kron(speye(n),A0p) + kron(Zn1,A1p) + kron(Zn1',A1p');
     
     inix = aggsolve_d1(A,para,Ap,A0inv,A1Zn1,A0inv_post,A1Zn1_post,P,b,choice,uors);
     
    mmax = 1000; % the maximal number of iterations 
    u = inix; 
    
    r = b - A * u;
    e(1) = norm(r)/norm(b);
    fprintf('\n Initial residual = %e',e(1));
    iter = 1;
    t1 = 1.0;
    d = zeros(n*l,1);
    
    while iter < mmax && e(iter) > tol
        z = aggsolve_d1(A,para,Ap,A0inv,A1Zn1,A0inv_post,A1Zn1_post,P,r,choice);
        
        u = u + z;
        r = r - A * z;
        iter = iter + 1;
         e(iter) = norm(r)/norm(b);
%         fprintf('\n at step %1.0f, relative residual = %e',... 
%            iter-1,e(iter));
        
    end
    
    fprintf('\n aggTGM at step %1.0f, relative residual = %e',... 
           iter-1,e(iter));

end