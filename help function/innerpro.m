function alpha = innerpro(A,B)
%--------------------------------------------------------------------------
% This function computes the inner product of two matrix.
%--------------------------------------------------------------------------
    
    alpha = sum(sum(A .* B));
    
end