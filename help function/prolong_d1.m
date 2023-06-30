function P = prolong_d1(n,l)
       
    z = 5;
    e = ones(l,1);
    Tp = eye(l) + ((z - 1)/l) * (e * e');
    
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    
    T = kron(speye(n),Tp) + kron(Zn1,Tp/2) + kron(Zn1',Tp/2);
    
    Kn = sparse((n -1)/2,n);
    Kn(:,2:end) = kron(speye((n - 1)/2),[1,0]);
    
    P = T * kron(Kn',speye(l));
     
end