function [X, resi] = pcgsylv(C, tol, para)
%--------------------------------------------------------------------------
%   PCG preconditioned conjugate gradient method for linear systems
%                       Wx=c
%   which is equivalent to a matrix equation
%               AX+BXJ+B'*X*J'=C
%   matrix A and B are tridiagonal Toeplitz matrices:
%   A=tridiag(para(2),para(1),para(2));
%   B=tridiag(para(3), para(4), para(5)); is unsymmetric
%   J=diag(ones(n-1,1),1)
%
%   In most case q=2*sum(para); such that it is harmonic, however, q 
%   can be any real number > 2*sum(para) such that W is SPD.
%
% inputs: 
%   C: the right-hand side
%   tol: tolerance for iteration
%   maxit: the maximum number of iterations
%   para: the parameters for the matrix, para>0

% Ref: Y. Saad. Iterative method for sparse linear systems, pp 263

% Copyright Yangfeng Su, 2020-02-15
% v1 Copyright Chengtao An, 2020-02-16: adapt this function to the FPD
% problem
%--------------------------------------------------------------------------

    %setup the matrices
    if length(para) ~= 5
        error('the length of the parameter vector p should be 5');
    end

    [m,n]=size(C);
    e = ones(m,1);
    A = spdiags([para(2)*e para(1)*e para(2)*e], -1:1,m,m);
    B = spdiags([para(3)*e para(4)*e para(5)*e],-1:1, m,m);
    BT = B';
    J = spdiags(ones(n,1), 1, n, n);
    JT = J';

    %setup the preconditioner parameter. In preconditioner, we symmetrize the
    %matrix B
    para1 = [para(1), para(2), (para(3)+para(5))/2, para(4)];
    
    maxit = 1000;
    k = 0;
    resi = zeros(maxit+1,1);
    X0 = TDFFTsol(C,para1);

    %Compute the initial residual
    R = C - ( A * X0 + B * X0 * J + BT * X0 * JT);
    resi(1) = norm(R,'fro');
    %fprintf('k=%d, resi=%.4e\n', k, resi(1))

    %preconditioner solver
    Z = TDFFTsol(R,para1);
    P = Z;
    X = X0;

    while(resi(k+1)>tol && k<maxit)
        % matrix-vector multiplication
        AP = A * P + B * P * J + BT * P * JT;

        alpha = innerpro(R,Z) /innerpro(AP,P);

        % update the iterate
        X = X + alpha * P;

        %update the residual
        R1 = R - alpha * AP;

        %precondition
        Z1 = TDFFTsol(R1,para1);

        %update the conjugate direction
        beta = innerpro(R1,Z1)/innerpro(R,Z);
        P = Z1 + beta * P;

        R = R1;
        Z = Z1;
        k = k+1;
        resi(k+1) = norm(R,'fro');
        %fprintf('k=%d, resi=%.4e\n', k, resi(k+1));    
    end
    resi = resi(1:k+1);
end

