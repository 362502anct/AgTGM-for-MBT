%% Qp lagrangian
clc
size_test = 5; % the size of the test case from 1 to size_test
iter_mat1 = zeros(size_test,3*4); % TGM
iter_mat2 = zeros(size_test,3*4); % Vcycle
iter_mat3 = zeros(size_test,3*4); % Vcycle PCG
iter_mat4 = zeros(size_test,3*4); % TGM PCG
time_mat1 = zeros(size_test,3*4); % TGM 
time_mat2 = zeros(size_test,3*4); % Vcycle
time_mat3 = zeros(size_test,3*4); % Vcycle PCG
time_mat4 = zeros(size_test,3*4); % TGM PCG

for case_num = 1:3
    switch case_num
        case 1
            l = 2;
            A0 = (1/3) * [16,-8;-8,14];
            A1 = (1/3) * [0,-8;0,1];
        case 2
            l = 3;
            A0 = (1/40) * [432,-297,54;
            -297,432 ,-189;
            54,-189,296];
            A1 = (1/40) * [0,0,-189;
            0, 0 ,54;
            0, 0 ,-13];
        case 3
            l = 4;
            A0 = (1/945) * [16640 ,-14208, 5888, -1472;
            -14208, 22320, -14208, 3048;
            5888 ,-14208, 16640 ,-6848;
            -1472 ,3048 ,-6848 ,9850];
            A1 = (1/945) * [0, 0, 0, -6848;
            0, 0, 0, 3048;
            0 ,0 ,0 ,-1472;
            0, 0, 0, 347];
    end

        for i = 1:size_test
        t =9+i;  
        n = 2^t - 1;
        b = rand(n*l,1);
        tol = 1e-7;
        smallsize = 1050;

        Zn1 = diag(ones(n - 1,1),1);
        Zn1 = sparse(Zn1);
        A1Zn1 = kron(Zn1,A1') + kron(Zn1',A1);
        A = kron(speye(n),A0) + A1Zn1;

        inix = kron(speye(n),inv(A0)) * b;
        p = ones(l,1)/sqrt(l);

        B = A0 - (A0 * p * p' * A0) / (p' * A0 * p);
        [lambdaf,vf,xf] = grad_mf_d1(A0,A1,A1',B,1);
        [U,V] = eig(A0);
        A0_5 = U * diag(sqrt(diag(V))) * U';
        invA0_5 = inv(A0_5);

        [lambdagmax,vgmax,xgmax] = grad_mf_d1(-eye(l),-invA0_5 * A1 * invA0_5,-invA0_5 * A1' * invA0_5,eye(l),0);
        lambdagmax = real(lambdagmax);
        omega = 2/( - lambdagmax);
        uors = 1;
        
        f = @(x) -eye(l) - invA0_5 * A1 * invA0_5 * exp(1i * x) - invA0_5 * A1' * invA0_5 * exp(-1i * x);
        g = @(x) eye(l) + invA0_5 * A1 * invA0_5 * exp(1i * x) + invA0_5 * A1' * invA0_5 * exp(-1i * x);
        
        fun_eig = @(x) compute_eig(f,x); 
        theta0 = [pi/2];
        [omega21,fval1,eflag,output] = patternsearch(fun_eig,theta0);
        gun_eig = @(x) compute_eig_2(g,x);
        [omega22,fval2,eflag,output] = patternsearch(gun_eig,theta0);

        omega = 2 / (-fval1);

        tic 
        [u,e,iter11] = Vcycleagg_d1(A0,A1,n,b,inix,tol,uors,p,omega,1000000);
        toc
        time11 = toc;

        tic 
        [u,e,iter12] = Vcycleagg_d1(A0,A1,n,b,inix,tol,uors,p,omega,smallsize);
        toc
        time12 = toc;
        
        tic
        [u,e,iter13] = pcgaggVcycle_d1(A0,A1,n,b,inix,tol,uors,p,omega,smallsize);
        toc
        time13 = toc;

        tic
        [u,e,iter14] = pcgaggVcycle_d1(A0,A1,n,b,inix,tol,uors,p,omega,1000000);
        toc
        time14 = toc;

        tic
        [u1,e1,iter21] = ToepVcyclerun_d1(A0,A1,n,b,tol,omega,1000000);
        toc
        time21 = toc;

        tic
        [u1,e1,iter22] = ToepVcyclerun_d1(A0,A1,n,b,tol,omega,smallsize);
        toc
        time22 = toc;
        
        tic
        [u1,e1,iter23] = ToepPCGVcycle_d1(A0,A1,n,b,tol,omega,smallsize);
        toc
        time23 = toc;
        
         tic
        [u1,e1,iter24] = ToepPCGVcycle_d1(A0,A1,n,b,tol,omega,1000000);
        toc
        time24 = toc;

        tic
        [u1,e1,iter31] = basisprojVcycle_d1(A0,A1,n,b,tol,omega,10000000);
        toc
        time31 = toc;
        
        tic
        [u1,e1,iter32] = basisprojVcycle_d1(A0,A1,n,b,tol,omega,smallsize);
        toc
        time32 = toc;
        
        tic
        [u1,e1,iter33] = basisprojVcycle_pcg(A0,A1,n,b,tol,omega,smallsize);
        toc
        time33 = toc;
        
        tic
        [u1,e1,iter34] = basisprojVcycle_pcg(A0,A1,n,b,tol,omega,10000000);
        toc
        time34 = toc;

        tic
        [u1,e1,iter41] = linearinterVcycle_d1(A0,A1,n,b,tol,omega,1000000);
        toc
        time41 = toc;

        tic
        [u1,e1,iter42] = linearinterVcycle_d1(A0,A1,n,b,tol,omega,smallsize);
        toc
        time42 = toc;
        
        tic
        [u1,e1,iter43] = linearinterVcycle_pcg(A0,A1,n,b,tol,omega,smallsize);
        toc
        time43 = toc;
        
        tic
        [u1,e1,iter44] = linearinterVcycle_pcg(A0,A1,n,b,tol,omega,10000000);
        toc
        time44 = toc;

        iter_mat1(i,case_num:3:end) = [iter11,iter21,iter31,iter41];        
        time_mat1(i,case_num:3:end) = [time11,time21,time31,time41];
        iter_mat2(i,case_num:3:end) = [iter12,iter22,iter32,iter42];        
        time_mat2(i,case_num:3:end) = [time12,time22,time32,time42];
        iter_mat3(i,case_num:3:end) = [iter13,iter23,iter33,iter43];        
        time_mat3(i,case_num:3:end) = [time13,time23,time33,time43];
        iter_mat4(i,case_num:3:end) = [iter14,iter24,iter34,iter44];        
        time_mat4(i,case_num:3:end) = [time14,time24,time34,time44];
        end
end

%% DG
%clc;

size_test = 5; % the size of the test case from 1 to size_test
result1_d2 = zeros(size_test,4); % TGM
result2_d2 = zeros(size_test,4); % Vcycle
result3_d2 = zeros(size_test,4); % TGM PCG
result4_d2 = zeros(size_test,4); % Vcycle PCG

A00 = [127/360,41/480,-43/320,41/480,-1/360,-2/45,-43/320,-2/45,13/288;
41/480,103/90,41/480,-1/360,5/24,-1/360,-2/45,-113/240,-2/45;
-43/320,41/480,127/360,-2/45,-1/360,41/480,13/288,-2/45,-43/320;
41/480,-1/360,-2/45,103/90,5/24,-113/240,41/480,-1/360,-2/45;
-1/360,5/24,-1/360,5/24,158/45,5/24,-1/360,5/24,-1/360;
-2/45,-1/360,41/480,-113/240,5/24,103/90,-2/45,-1/360,41/480;
-43/320,-2/45,13/288,41/480,-1/360,-2/45,127/360,41/480,-43/320;
-2/45,-113/240,-2/45,-1/360,5/24,-1/360,41/480,103/90,41/480;
13/288,-2/45,-43/320,-2/45,-1/360,41/480,-43/320,41/480,127/360];

A_10 = [5/288,5/576,-5/1152,23/720,23/1440,-23/2880,-11/1440,-11/2880,11/5760;
5/576,5/72,5/576,23/1440,23/180,23/1440,-11/2880,-11/360,-11/2880;
-5/1152,5/576,5/288,-23/2880,23/1440,23/720,11/5760,-11/2880,-11/1440;
-17/144,-17/288,17/576,-47/360,-47/720,47/1440,23/720,23/1440,-23/2880;
-17/288,-17/36,-17/288,-47/720,-47/90,-47/720,23/1440,23/180,23/1440;
17/576,-17/288,-17/144,47/1440,-47/720,-47/360,-23/2880,23/1440,23/720;
-7/288,-7/576,7/1152,-17/144,-17/288,17/576,5/288,5/576,-5/1152;
-7/576,-7/72,-7/576,-17/288,-17/36,-17/288,5/576,5/72,5/576;
7/1152,-7/576,-7/288,17/576,-17/288,-17/144,-5/1152,5/576,5/288];

A0_1 = [5/288,23/720,-11/1440,5/576,23/1440,-11/2880,-5/1152,-23/2880,11/5760;
-17/144,-47/360,23/720,-17/288,-47/720,23/1440,17/576,47/1440,-23/2880;
-7/288,-17/144,5/288,-7/576,-17/288,5/576,7/1152,17/576,-5/1152;
5/576,23/1440,-11/2880,5/72,23/180,-11/360,5/576,23/1440,-11/2880;
-17/288,-47/720,23/1440,-17/36,-47/90,23/180,-17/288,-47/720,23/1440;
-7/576,-17/288,5/576,-7/72,-17/36,5/72,-7/576,-17/288,5/576;
-5/1152,-23/2880,11/5760,5/576,23/1440,-11/2880,5/288,23/720,-11/1440;
17/576,47/1440,-23/2880,-17/288,-47/720,23/1440,-17/144,-47/360,23/720;
7/1152,17/576,-5/1152,-7/576,-17/288,5/576,-7/288,-17/144,5/288];

A10 = A_10';
A01 = A0_1';

l = 9;
A = eye(l);
[U,V] = eig(A00);
A0_5 = U * diag(sqrt(diag(V))) * U';
invA0_5 = inv(A0_5);
R = invA0_5 * A10 * invA0_5;
C = invA0_5 * A01' * invA0_5;

f = @(x) -A -R * exp(-1i * x(1)) - R' * exp(1i * x(1)) - ...
     C * exp(-1i * x(2)) - C' * exp(1i * x(2));
g = @(x) A + R * exp(-1i * x(1)) + R' * exp(1i * x(1)) + ...
     C * exp(-1i * x(2)) + C' * exp(1i * x(2));

fun_eig = @(x) compute_eig(f,x); 
theta0 = [pi, pi];
[omega21,fval,eflag,output] = patternsearch(fun_eig,theta0);
gun_eig = @(x) compute_eig_2(g,x);
[omega22,fval2,eflag,output] = patternsearch(gun_eig,theta0);

omega1 = 2 / (-fval );

Pl = ones(l,1)/sqrt(l);
A00p = Pl' * A00 * Pl;
A01p = Pl' * A01 * Pl;
A10p = Pl' * A10 * Pl;
A0_1p = Pl' * A0_1 * Pl;
A_10p = Pl' * A_10 * Pl;
     
f = @(x) - 1 - A10p/A00p * exp(-1i * x(1)) - A10p/A00p * exp(1i * x(1)) - ...
     A01p/A00p * exp(-1i * x(2)) - A01p/A00p * exp(1i * x(2));
g = @(x) 1 + A10p/A00p * exp(-1i * x(1)) + A10p/A00p * exp(1i * x(1)) + ...
     A01p/A00p * exp(-1i * x(2)) + A01p/A00p * exp(1i * x(2));

theta0 = [pi, pi];
[omega21,fval11,eflag,output] = patternsearch(f,theta0);
gun_eig = @(x) compute_eig_2(g,x);
[omega22,fval12,eflag,output] = patternsearch(g,theta0);

omega11 = 2 / (-fval11);

l = 9;


for i = 1:size_test
    t = 4+i;
    n = 2^t - 1;
    m = n;
    b = rand(l*n*m,1);
    tol = 1e-7;
    Zn1 = diag(ones(n - 1,1),1);
    Zn1 = sparse(Zn1);
    Zm1 = diag(ones(m - 1,1),1);
    Zm1 = sparse(Zm1);
    A0 = kron(speye(m),A00) + kron(Zm1,A10) + kron(Zm1',A_10);
    A1 = kron(speye(m),A01);

    p = ones(l,1)/sqrt(l);
    uors = 1;

    tic
    [u,e,iter11, time11] = Vcycleagg_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,uors,p,omega1,omega11,100000000);
    toc
    time11 = toc;
    tic
    [u1,e1,iter12,time12] = ToepVcyclerun_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,omega1,100000000);
    toc
    time12 = toc;
    result1_d2(i,:) = [iter11-1,time11,iter12-1,time12];
    
    tic
    [u,e,iter21,time21] = Vcycleagg_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,uors,p,omega1,omega11,10000);
    toc
    time21 = toc;
    tic
    [u1,e1,iter22,time22] = ToepVcyclerun_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,omega1,10000);
    toc
    time22 = toc;
    result2_d2(i,:) = [iter21-1,time21,iter22-1,time22];
    
    tic
    [u,e,iter31,time31] = pcgaggVcycle_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,uors,p,omega1,omega11,10000);
    toc
    time31 = toc;
    tic
    [u1,e1,iter32,time32] = ToepPCGVcycle_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,omega1,10000);
    toc
    time32 = toc;
    result3_d2(i,:) = [iter31 - 1, time31, iter32-1, time32];
    
    tic
    [u,e,iter41,time41] = pcgaggVcycle_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,uors,p,omega1,omega11,100000000);
    toc
    time41 = toc;
    tic
    [u1,e1,iter42,time42] = ToepPCGVcycle_d2(A00,A01,A10,A0_1,A_10,n,m,b,tol,omega1,100000000);
    toc
    time42 = toc;
    result4_d2(i,:) = [iter41 - 1, time41, iter42-1, time42];
end

%% Test FPD problem
%clc; clear;
size_test = 5;

result = cell(3,4);

for case_num = 1:3
   switch case_num
       case 1
           load('2k.mat');
            nvsource = 3240;
            interaction = 1;
            m = 1920;%row
            n = 1080;%col
            k = 7;
            omega = 0.8;
       case 2
           load('640_360.mat');
            nvsource = 1440;
            interaction = 2;
            m = 640;%row
            n = 360;%col
            k = 14;
            omega = 0.9;           
       case 3
            load('4k.mat');
            nvsource = 129602;
            interaction = 4;
            m = 3840;%row
            n = 2160;%col
            k = 28;
            omega = 1.1;
   end
   
    A = A11_cf;
    R = A12_cf;
    C = B11_cf;
    D1 = B21_cf;
    [l,~] = size(A);
    D2 = zeros(l,l);

    A0 = eye(l);
    [U,V] = eig(full(A));
    A0_5 = U * diag(sqrt(diag(V))) * U';
    invA0_5 = inv(A0_5);
    Chol_A = chol(A, 'lower');
    invChol_A = inv(full(Chol_A));
    R0 = invChol_A * R * invChol_A';
    C0 = invChol_A * C' * invChol_A';
    D10 = invChol_A * D1 * invChol_A';

    F = @(x,y) A + R * exp(-1i * x) + R' * exp(1i * x) + ...
         C * exp(-1i * y) + C' * exp(1i * y)  + ...
         D1 * exp(-1i * (x - y)) + D1' * exp(1i * (x - y));

    f = @(x) -A0 -R0 * exp(-1i * x(1)) - R0' * exp(1i * x(1)) - ...
         C0 * exp(-1i * x(2)) - C0' * exp(1i * x(2)) - ...
         D10 * exp(-1i * (x(1) - x(2))) - D10' * exp(1i * (x(1) - x(2)));
    g = @(x) A0 + R0 * exp(-1i * x(1)) + R0' * exp(1i * x(1)) + ...
         C0 * exp(-1i * x(2)) + C0' * exp(1i * x(2)) + ...
         D10 * exp(-1i * (x(1) - x(2))) + D10' * exp(1i * (x(1) - x(2)));
     
    fun_eig = @(x) compute_eig(f,x); 
    theta0 = [pi/2, 0];
    [omega21,fval,eflag,output] = patternsearch(fun_eig,theta0);
    gun_eig = @(x) compute_eig_2(g,x);
    [omega22,fval2,eflag,output] = patternsearch(gun_eig,theta0);

    omega1 = (-fval) / 2;    
    
     P = ones(1,l)/sqrt(l);
    AP = P * A * P';
    CP = P * C * P';
    RP = P * R * P';
    D1P = P * D1 * P';
    
    f = @(x) -1 -RP/AP * exp(-1i * x(1)) - RP/AP * exp(1i * x(1)) - ...
         CP/AP * exp(-1i * x(2)) - CP/AP * exp(1i * x(2)) - ...
         D1P/AP * exp(-1i * (x(1) - x(2))) - D1P/AP * exp(1i * (x(1) - x(2)));
    
    g = @(x) 1 + RP/AP * exp(-1i * x(1)) + RP/AP * exp(1i * x(1)) + ...
         CP/AP * exp(-1i * x(2)) + CP/AP * exp(1i * x(2)) + ...
         D1P/AP * exp(-1i * (x(1) - x(2))) + D1P/AP * exp(1i * (x(1) - x(2)));

    theta0 = [pi, pi/2];
    [omega21,fval11,eflag,output] = patternsearch(f,theta0);
    [omega22,fval12,eflag,output] = patternsearch(g,theta0);

    omega11 = (-fval11) / 2;

    result1 = zeros(size_test,4);
    result2 = zeros(size_test,4);
    result3 = zeros(size_test,4);
    result4 = zeros(size_test,4);


    for i = 1:size_test
    t = 4+i;
    m = 2^t - 1;
    n = m;

    Zm1 = spdiags(ones(m,1),1,m,m);
    Zn1 = spdiags(ones(n,1),1,n,n);

    bmk = rand(m*k,1);
    b = kron(ones(n,1),bmk);
    tol = 1e-7;
    step = 1;

    p = ones(1,l)/sqrt(l);

    tic
    [U,e,iter11] = Vcycleagg_fpd(A,C,R,D1,D2,n,m,b,tol,omega1,omega11,100000000);
    toc
    time11 = toc;

    tic
    [u,e,iter12] = ToepVcyclerun_fpd(A,C,R,D1,D2,n,m,b,tol,omega1,10000000);
    toc
    time12 = toc;
    result1(i,:) = [iter11 - 1,time11,iter12 - 1,time12];


    tic
    [U,e,iter21] = Vcycleagg_fpd(A,C,R,D1,D2,n,m,b,tol,omega1,omega11,10000);
    toc
    time21 = toc;

    tic
    [u,e,iter22] = ToepVcyclerun_fpd(A,C,R,D1,D2,n,m,b,tol,omega1,10000);
    toc
    time22 = toc;
    result2(i,:) = [iter21 - 1,time21,iter22 - 1,time22];


     tic
    [U,e,iter31] = pcgaggVcycle_fpd(A,C,R,D1,D2,n,m,b,tol,omega1,omega11,10000);
    toc
    time31 = toc;
     tic
    [u,e,iter32] = ToepPCGVcycle_fpd(A,C,R,D1,D2,n,m,b,tol,omega1,10000);
    toc
    time32 = toc;
    result3(i,:) = [iter31 - 1, time31, iter32 - 1, time32];

    tic
    [U,e,iter41] = pcgaggVcycle_fpd(A,C,R,D1,D2,n,m,b,tol,omega1,omega11,10000000000);
    toc
    time41 = toc;
     tic
    [u,e,iter42] = ToepPCGVcycle_fpd(A,C,R,D1,D2,n,m,b,tol,omega1,10000000000);
    toc
    time42 = toc;
    result4(i,:) = [iter41 - 1, time41, iter42 - 1, time42];

    end
    
    result(case_num,1) = {result1};
    result(case_num,2) = {result2};
    result(case_num,3) = {result3};
    result(case_num,4) = {result4};

end


%%
function omega = compute_omega(f)

    g = @(theta) compute_eig(f, theta);
    theta0 = [pi/2, pi/2];
    [omega,fval,eflag,output] = patternsearch(g,theta0);

end

function mineig = compute_eig(f,theta)

     v = eig(f(theta));
     mineig = real(v(1));

end

function mineig = compute_eig_2(f,theta)

     v = eig(f(theta));
     mineig = real(v(2));

end

