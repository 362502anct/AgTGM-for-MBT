function [K0,K1,K_1] = linearmatrix(l)
    K0 = zeros(l,l);
    K1 = K0;
    K_1 = K0;
    K1(1,end) = 1;
    
    if rem(l,2)
        A2 = zeros(l,(l+1)/2);
        A3 = zeros(l,(l+1)/2);
        for i = 1:(l - 1)/2
            A2(2 * i - 1,i) = 2;
            A2(2 * i,i:i+1) = [1,1];
            A3(2*i,i) = 2;
            A3(2*i+1,i:i+1) = [1,1];
        end
        A2(end,end) = 2;
        A3(1,1) = 1;
        K0(:,(l+1)/2:end) = A2;
        K_1(:,1:(l+1)/2) = A3;
    else
        A1 = zeros(l,l/2);
        A1(1,1) = 1;
        A1(end,end) = 2;
        if l/2 > 1
            for i = 1:l/2 - 1
                A1(2*i,i) = 2;
                A1(2*i+1,i:i+1) = [1,1];
            end
        end
        K_1(:,1:l/2) = A1;
        K0(:,l/2+1:end) = A1;
        K0(1,l/2) = 1;
    end
end

