function X = TDFFTsol(C,p)
%--------------------------------------------------------------------------
% solve 2D Poisson-like equation with FFT
% W=block-tridiag(B, A, B), B is symmetric
% A=tridiag(p(2),p(1),p(2)); B=tridiag(p(3), p(4), p(3));
% The equation can be written as:
%  AX+BXJ=C
%  J=tridiag(1,0,1).

% Copyright Yangfeng Su, 2020-02-15
% v1 Copyright Chengtao An, 2020-02-16: adapt this function to the FPD
% problem
% v2 Copyright Chengtao An, 2020-02-17: add the discrete sine transform
% function dst
%--------------------------------------------------------------------------

    [m,n] = size(C);

    %diagonalize A and B
    X = conj(dst(conj(dst(C)'))');
    %X = conj(Sint(conj(Sint(C)'))');
    
    %soler the diagonal system
    emJ = 2 * cos(pi/(m+1)*(1:m)');
    eA = p(1) + p(2) * emJ;
    eB = p(4) + p(3) * emJ;
    enJ = 2 * cos(pi/(n+1)*(1:n)');
    t = eA * ones(1,n) + eB * enJ';
    X = X ./ t;

    %backward substitution    
    %X = conj(Sint(conj(Sint(C)'))');
    X = conj(dst(conj(dst(X)'))');
    X = X / ((m + 1) * (n + 1) / 4);

end

