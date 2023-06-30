function P = basisprojprolong_d1(n,l)
       
    [K0,K1,K_1,K2,K_2] = basisproj_symbol(l);
    
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    Zn2 = diag(ones(n - 2,1),2);
    Zn2 = sparse(Zn2);
    
    T = kron(speye(n),K0) + kron(Zn1,K_1) + kron(Zn1',K1) + kron(Zn2,K_2) + kron(Zn2',K2);
    
    Kn = sparse((n -1)/2,n);
    Kn(:,2:end) = kron(speye((n - 1)/2),[1,0]);
    
    P = T * kron(Kn',speye(l));
     
end
