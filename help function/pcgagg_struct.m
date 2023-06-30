function [U,e] = pcgagg_struct(A,C,R,D1,D2,n,m,b,tol,step,P,omega,choice)
%--------------------------------------------------------------------------
% This function is used to solve the structured linear system Gx = b where
% G has the form 
%     G = kron(eye(n),kron(eye(m),A) + kron(Zm1,C) + kron(Zm1',C')) + 
%         kron(Zn1, kron(eye(m),R) + kron(Zm1', D1)) + 
%         kron(Zn1', kron(eye(m),R')+ kron(Zm1,D1')),
% by the aggregation based multigrid preconditioner conjugate gradient method.
% 
% n: the number of pixels in one row;
% m: the number of pixels in one column;
% tol: the tolerance;
% step: the number of pre- and post- smoother step;

% Copyright Chengtao An, 2020-02-09
% v1 Copyright Chengtao An, 2020-02-10: add the pre-smooth step and the
% recommanded step is 1;
% v2 Copyright Chengtao An, 2020-02-16: use the matrix form to replace the
% vector form;
%--------------------------------------------------------------------------
    
    [l,~] = size(A);
    Zn1 = diag(ones(n - 1,1),1);
    Zm1 = diag(ones(m - 1,1),1);
    Zn1 = sparse(Zn1);
    Zm1 = sparse(Zm1);
    ml = m*l;
    
    G1 = kron(speye(m),A) + kron(Zm1,C) + kron(Zm1',C') ;
    G2 = kron(speye(m),R) + kron(Zm1', D1);
    AG = kron(speye(n),G1) + kron(Zn1,G2) + kron(Zn1',G2');
    
    P = ones(1,l)/sqrt(l);
    P = sparse(P);

    AP = P * A * P';
    CP = P * C * P';
    RP = P * R * P';
    D1P = P * D1 * P';
    D2P = P * D2 * P';
    Para = [AP CP D2P RP D1P];
    
    Pmt = kron(eye(m),P');
    Pm = kron(eye(m),P);
    
    B = reshape(b,ml,n);
    
    G1p = kron(speye(m),AP) + kron(Zm1,CP) + kron(Zm1',CP') ;
    G2p = kron(speye(m),RP) + kron(Zm1', D1P);
    Ap = kron(speye(n),G1p) + kron(Zn1,G2p) + kron(Zn1',G2p');
    Lp = chol(Ap,'lower');
    %     %
%     Adiag = sparse(diag(diag(A)));
%     A1 = kron(eye(m),Adiag);
%     GA1 = G1 - A1;
%     Adiaginv = sparse(diag(1./diag(A)));
%     A1inv = kron(eye(m),Adiaginv);
%     E = C + D1 + R + C' + D1' + R';
%     e1 = sum(E);
%     Ae = sparse(A + diag(e1) - max(e1) * eye(l));
    Aja = A;
    A1 = kron(speye(m),Aja);
    GA1 = G1 - A1;
    Ainv = sparse(inv(full(Aja)));
    A1inv = kron(speye(m),Ainv);
   iniX = aggsolve_jacobi(Para,Pm,Pmt,GA1,G1,G2,A1inv,Zn1,B,step,omega,Ap,choice,Lp);
   % iniX = aggsolve(Para,Pm,Pmt,G1,G2,LG1,LG1t,Zn1,B,step);
    
    mmax = 10000; % the maximal number of iterations 
    U = iniX; 
    
    R = B - Gmultiple(G1,G2,Zn1,iniX);
    e(1) = norm(R,'fro')/norm(B,'fro');
    fprintf('\n Initial residual = %e',e(1));
    iter = 1;
    t1 = 1.0;
    D = zeros(ml,n);
    
    while iter < mmax && e(iter) > tol
       Z = aggsolve_jacobi(Para,Pm,Pmt,GA1,G1,G2,A1inv,Zn1,R,step,omega,Ap,choice,Lp);
       % Z = aggsolve(Para,Pm,Pmt,G1,G2,LG1,LG1t,Zn1,R,step);
        
        Z = real(Z);
        t1old = t1;
        t1 = innerpro(Z,R);
        beta = t1/t1old; 
        D = Z+beta*D; 
        
        S = Gmultiple(G1,G2,Zn1,D);
        
        suma = innerpro(D,S);
        tau = t1/suma;
        U = U+tau*D;
        R = R-tau*S; 
        iter = iter+1; 
        e(iter) = norm(R,'fro')/norm(B,'fro');
        fprintf('\n at step %1.0f, relative residual = %e',... 
            iter-1,e(iter));
        
    end
end