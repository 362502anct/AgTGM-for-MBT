function P = linearinter_prolong(n,l)

    [K0,K1,K_1] = linearmatrix(l);
    
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    
    T = kron(speye(n),K0) + kron(Zn1,K_1) + kron(Zn1',K1);
    
    Kn = sparse((n -1)/2,n);
    Kn(:,2:end) = kron(speye((n - 1)/2),[1,0]);
    
    P = T * kron(Kn',speye(l));
     
end