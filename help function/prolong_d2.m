function P = prolong_d2(n,m,l)
       
    z = 5;
    e = ones(l,1);
    Tp = eye(l) + ((z - 1)/l) * (e * e');
    
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    Zm1 = diag(ones(m - 1,1),1);
    Zm1 = sparse(Zm1);
    
    Tp1 = kron(speye(m),Tp) + kron(Zm1,Tp/2) + kron(Zm1',Tp/2);
    Tp2 = kron(speye(m),Tp/2) + kron(Zm1,Tp/4) + kron(Zm1',Tp/4);
    T = kron(speye(n),Tp1) + kron(Zn1,Tp2) + kron(Zn1',Tp2');
    
    Kn = sparse((n -1)/2,n);
    Kn(:,2:end) = kron(speye((n - 1)/2),[1,0]);
    Km = sparse((m -1)/2,m);
    Km(:,2:end) = kron(speye((m - 1)/2),[1,0]);
    
    P = T * kron(Kn',kron(Km',speye(l)));
     
end