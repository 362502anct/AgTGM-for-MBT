function  [K0,K1,K_1,K2,K_2] = basisproj_symbol(l)

    switch l
        case 2
            % Q2
            node = [0,1,2]/2;
            phi1 = @(x) (x - node(2)) .* (x - node(3))  ./ ...
                ((node(1) - node(2)) .* (node(1) - node(3)) );
            phi2 = @(x) (x - node(1)) .* (x - node(3))   ./ ...
                ((node(2) - node(1)) .* (node(2) - node(3)) );
            phi3 = @(x) (x - node(1)) .* (x - node(2))  ./ ...
                ((node(3) - node(1)) .* (node(3) - node(2)) );

            node2 = [1,2,3,4]/4;
            phi1_2 = phi1(node2);
            phi2_2 = phi2(node2);
            phi3_2 = phi3(node2);

            size = 2;
            P = zeros(7,3);
            P(1:size * 2,1) = phi2_2;
            P(1:size * 2,2) = phi3_2;
            P(2*size + 1:end,2) = phi1_2(1:size * 2 - 1);
            P(2*size + 1:end,3) = phi2_2(1:size * 2 - 1);


            K_1 = P(1:size,1:size);
            K0 = P(size + 1:size *2,1:size);
            K_2 = zeros(size,size);
            K1 = P(size*2 + 1:size * 3,1:size);
            K2 = P(size*3 + 1:end,1:size);
            K2(size,:) =0;
        case 3
            % Q3
            node = [0,1,2,3]/3;
            phi1 = @(x) (x - node(2)) .* (x - node(3)) .* (x - node(4)) ./ ...
                ((node(1) - node(2)) .* (node(1) - node(3)) .* (node(1) - node(4)));
            phi2 = @(x) (x - node(1)) .* (x - node(3)) .* (x - node(4))  ./ ...
                ((node(2) - node(1)) .* (node(2) - node(3)) .* (node(2) - node(4)));
            phi3 = @(x) (x - node(1)) .* (x - node(2)) .* (x - node(4))  ./ ...
                ((node(3) - node(1)) .* (node(3) - node(2)) .* (node(3) - node(4)));
            phi4 = @(x) (x - node(1)) .* (x - node(2)) .* (x - node(3))  ./ ...
                ((node(4) - node(1)) .* (node(4) - node(2)) .* (node(4) - node(3)));

            node2 = [1,2,3,4,5,6]/6;
            phi1_2 = phi1(node2);
            phi2_2 = phi2(node2);
            phi3_2 = phi3(node2);
            phi4_2 = phi4(node2);

            size = 3;
            P = zeros(11,5);
            P(1:size * 2,1) = phi2_2;
            P(1:size * 2,2) = phi3_2;
            P(1:size * 2,3) = phi4_2;
            P(2*size + 1:end,3) = phi1_2(1:size * 2 - 1);
            P(2*size + 1:end,4) = phi2_2(1:size * 2 - 1);
            P(2*size + 1:end,5) = phi3_2(1:size * 2 - 1);

            K_1 = P(1:size,1:size);
            K0 = P(size + 1:size *2,1:size);
            K_2 = zeros(size,size);
            K1 = P(size*2 + 1:size * 3,1:size);
            K2 = P(size*3 + 1:end,1:size);
            K2(size,:) =0;
        case 4
            % Q4
            node = [0,1,2,3,4]/4;
            phi1 = @(x) (x - node(2)) .* (x - node(3)) .* (x - node(4)) .* (x - node(5)) ./ ...
                ((node(1) - node(2)) .* (node(1) - node(3)) .* (node(1) - node(4)) .* (node(1) - node(5)));
            phi2 = @(x) (x - node(1)) .* (x - node(3)) .* (x - node(4)) .* (x - node(5)) ./ ...
                ((node(2) - node(1)) .* (node(2) - node(3)) .* (node(2) - node(4)) .* (node(2) - node(5)));
            phi3 = @(x) (x - node(1)) .* (x - node(2)) .* (x - node(4)) .* (x - node(5)) ./ ...
                ((node(3) - node(1)) .* (node(3) - node(2)) .* (node(3) - node(4)) .* (node(3) - node(5)));
            phi4 = @(x) (x - node(1)) .* (x - node(2)) .* (x - node(3)) .* (x - node(5)) ./ ...
                ((node(4) - node(1)) .* (node(4) - node(2)) .* (node(4) - node(3)) .* (node(4) - node(5)));
            phi5 = @(x) (x - node(1)) .* (x - node(2)) .* (x - node(3)) .* (x - node(4)) ./ ...
                ((node(5) - node(1)) .* (node(5) - node(2)) .* (node(5) - node(3)) .* (node(5) - node(4)));

            node2 = [1,2,3,4,5,6,7,8]/8;
            phi1_2 = phi1(node2);
            phi2_2 = phi2(node2);
            phi3_2 = phi3(node2);
            phi4_2 = phi4(node2);
            phi5_2 = phi5(node2);

            P = zeros(15,7);
            P(1:8,1) = phi2_2;
            P(1:8,2) = phi3_2;
            P(1:8,3) = phi4_2;
            P(1:8,4) = phi5_2;
            P(9:15,4) = phi1_2(1:7);
            P(9:15,5) = phi2_2(1:7);
            P(9:15,6) = phi3_2(1:7);
            P(9:15,7) = phi4_2(1:7);
            K_1 = P(1:4,1:4);
            K0 = P(5:8,1:4);
            K_2 = zeros(4,4);
            K1 = P(9:12,1:4);
            K2 = P(13:end,1:4);
            K2(4,:) =0;
    end

end

