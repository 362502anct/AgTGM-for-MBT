function  x = dst_d1(b,p)
    
    m = length(b);
    x = dst(b);
    
    %soler the diagonal system
    emJ = 2 * cos(pi/(m+1)*(1:m)');
    t = p(1) + p(2) * emJ;
    x = x ./ t;

    %backward substitution    
    x = dst(x);
    x = x / ((m + 1) / 2);

end