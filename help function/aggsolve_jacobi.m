function X = aggsolve_jacobi(Para,Pm,Pmt,GA1,G1,G2,A1inv,Zn1,B,step,omega,Ap,choice,Lp)
%--------------------------------------------------------------------------
% This function can be seen as a v-cycle of the aggregation-based  
% multigrid method. Naturally, this function can be used as a 
% preconditioner.
%
% b: the right hand side term
% Para: the parameters which represent the aggregated linear system;
% Pm, Pmt: the prolongation matrix and its transpose;
% We use the function  to solve the aggregated linear system;
% 
% As for the smooth procedure, we use the block jacobi iteration. Since
% the original matrix can be written as 
%   G = kron(eye(n),G1) + kron(Zn1,G2) + kron(Zn1',G2')
% And the smooth step is as follow
%   x = kron(eye(n),G1)\(b - (kron(Zn1,G2) + kron(Zn1',G2')) * x0)
% Also, to accerlate the speed, we compute the cholesky decomposition 
% of G1 as G1 = LG1 * LG1t.
%
% ml: the size of G1 and G2;
% step: the number of pre- and post- smoother step;

% Copyright Chengtao An, 2020-02-09
% v1 Copyright Chengtao An, 2020-02-10: add the pre-smooth step and the
% recommanded step is 1;
% v2 Copyright Chengtao An, 2020-02-16: use the matrix form to replace the
% vector form;
%--------------------------------------------------------------------------
    
   [ml,n] = size(B);
    X = zeros(ml,n);
    for i = 1:step
        RB = B - Gmultiple(GA1,G2,Zn1,X);
        X1 = A1inv * RB;
        X = X + omega * (X1 - X);
    end
    
    R1 = B - Gmultiple(G1,G2,Zn1,X);

    tol = 1e-13;
    R1P = Pm * R1;
    [nm,~] = size(Ap);
    m = nm/n;
    
    if choice == 'T'
        r1p = reshape(R1P,m *n,1);
        x2p = Lp'\(Lp\r1p);
        X2P = reshape(x2p,m,n);
        %[X2P, ~] = pcgsylv(R1P, tol, Para);
    elseif choice == 'V'
        omega1 = 7/8;
        omega2 = 7/12;
        nsize = 1000;
        r1P = reshape(R1P,n*m,1);
        x2P = ToepVcycle_d2(Ap,r1P,n,m,1,nsize,omega1,omega2);
        X2P = reshape(x2P,m,n);
    end
    
%     r1P = reshape(R1P,n*m,1);
%     x2P = GP\r1P;
%     X2P = reshape(x2P,m,n);
    X2 = Pmt * X2P;  
    
    X = X + X2;
    
    for i = 1:step
        RB2 = B - Gmultiple(GA1,G2,Zn1,X);
        X3 = A1inv * RB2;
        X = X + omega * (X3 - X);
    end  
end
