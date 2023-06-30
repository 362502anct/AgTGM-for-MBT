function X = Gmultiple(G1,G2,Zn1,B)
%--------------------------------------------------------------------------
% This function is used to compute the matrix-vector multiplication
% as 
%       x = [kron(eye(n),G1) + kron(Zn1,G2) + kron(Zn1',G2')] * b
% where n is the size of Zn1 and ml is the size of G2.

% Copyright Chengtao An, 2020-02-09
% v2 Copyright Chengtao An, 2020-02-16: use the matrix form to replace the
% vector form;
%--------------------------------------------------------------------------
    
    X = G1 * B + G2 * B * Zn1' + G2' * B * Zn1;    
end