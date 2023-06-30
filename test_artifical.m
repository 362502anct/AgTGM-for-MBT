%% Well condition case
clc;
l = 16;
rng(200);
A = rand(l,l);
A = A + A';
for i = 1:l
    for j = 1:l
        if i == j
            A(i,j) = abs(A(i,j));
        else
            A(i,j) = -abs(A(i,j));
        end
    end
end
[U,V] = eig(A);
A0 = A - 2 * V(1,1) * eye(l);
[U,V] = eig(A0);
p0 = U(:,1);
A1 = zeros(l,l);
A1(:,l) = (1/4) * A0(:,l);
A1(l,l) = -A1(l,l);
A0 = A0 - A1 - A1';
D = diag(diag(A0));

[U,V] = eig(A0);
invA0_5 = U * diag(1./sqrt(diag(V))) * U';
L = chol(A0,'lower');
Linv = inv(L);
Linvt = inv(L');
R = invA0_5 * A1 * invA0_5;
Rt = invA0_5 * A1' * invA0_5;
f= @(x)  eye(l) + R * exp(-1i * 1 * x) + Rt * exp(1i * 1 * x);
g = @(x,pt) f(x) - (f(x) * pt) * pt' * f(x) / (pt' * f(x) * pt);


% Find a good p
[lambdaf,vf,xf,p] = findp_d1(A0,A1);
p = real(p);
p_opt = p/sqrt(p'*p);

ft = f;
fun_p = @(ps) eig2(ps,ft,l);

%[p_opt,fval,eflag,output] = patternsearch(fun_p,p_opt);
%p_opt = p_opt/norm(p_opt);

for i = 1:5

n = 2^(5+i);
b = rand(n*l,1);
tol = 1e-13;
inix = rand(n*l,1);

[U,V] = eig(A0);
A0_5 = U * diag(sqrt(diag(V))) * U';
invA0_5 = inv(A0_5);

[lambdagmax,vgmax,xgmax] = grad_mf_d1(-eye(l),-invA0_5 * A1 * invA0_5,-invA0_5 * A1' * invA0_5,eye(l),0);
lambdagmax = real(lambdagmax);
[lambdagmin,vgmin,xgmin] = grad_mf_d1(eye(l),invA0_5 * A1 * invA0_5,invA0_5 * A1' * invA0_5,eye(l),1);
lambdagmin = real(lambdagmin);
omega = 2 / (lambdagmin - lambdagmax);
omega1 = -2 / lambdagmax;
omega = real(omega);
omega1 = real(omega1);

[u,e] = twogridagg_d1(A0,A1,n,b,inix,tol,1,p_opt,omega);

p1 = ones(l,1)/sqrt(l);
[u,e] = twogridagg_d1(A0,A1,n,b,inix,tol,1,p1,omega);

end

%% Construct Ill condition
clc
l = 4;
rng(100);
A = rand(l,l);
A = A + A';
for i = 1:l
    for j = 1:l
        if i == j
            A(i,j) = abs(A(i,j));
        else
            A(i,j) = -abs(A(i,j));
        end
    end
end
[U,V] = eig(A);
A0 = A - V(1,1) * eye(l);
[U,V] = eig(A0);
p = U(:,1);
A1 = zeros(l,l);
A1(:,l - 1) = (1/4) * A0(:,l - 1);
A1(l - 1,l - 1) = -A1(l - 1,l - 1);
A0 = A0 - A1 - A1';

%step = 0;
step = 0:1e-5:1e-4;
theta = zeros(size(step));
iter = zeros(length(step),5);

for i = 5:9

n = 2^(5+i);
b = rand(n*l,1);
inix = rand(n*l,1);
tol = 1e-8;

Zn1 = diag(ones(n - 1,1),1);
Zn1 = sparse(Zn1);
A1Zn1 = kron(Zn1,A1') + kron(Zn1',A1);
A = kron(speye(n),A0) + A1Zn1;

% Find a good p
B = A0 - (A0 * p * p' * A0) / (p' * A0 * p);
[lambdaf,vf,xf] = grad_mf_d1(A0,A1,A1',B,1);
[U,V] = eig(A0);
A0_5 = U * diag(sqrt(diag(V))) * U';
invA0_5 = inv(A0_5);

[lambdagmax,vgmax,xgmax] = grad_mf_d1(-eye(l),-invA0_5 * A1 * invA0_5,-invA0_5 * A1' * invA0_5,eye(l),0);
lambdagmax = real(lambdagmax);
[lambdagmin,vgmin,xgmin] = grad_mf_d1(eye(l),invA0_5 * A1 * invA0_5,invA0_5 * A1' * invA0_5,eye(l),1);
lambdagmin = real(lambdagmin);
alpha = 0.5;
omega = 2/( - lambdagmax);
omega = real(omega);

    perturb = rand(l,1);
    for j = 1:length(step)
        p1 = p + step(j) * perturb;
        p1 = p1/norm(p1);
        theta(j) = p1' * p;
        pg = ones(l,1)/sqrt(l);
        [u,e] = twogridagg_d1(A0,A1,n,b,inix,tol,1,p1,omega);
    %    [ug,eg] = twogridagg_d1(A0,A1,n,b,inix,tol,0,pg,omega); 
    % Notay selection
        iter(j,i - 4) = length(e) - 1;
    end
    
end
%%
function val = eig2(ps,ft,l)
    g = @(x,pt) ft(x) - (ft(x) * pt) * pt' * ft(x) / (pt' * ft(x) * pt);
    E = @(theta,p) sort(eig(g(theta,p)));
    e2 = zeros(l,1);
    e2(2) = 1;
    eig_2nd = @(theta,p) e2' * real(E(theta,p));
    p = ps;
    fun_eig = @(theta) eig_2nd(theta,p);
    theta0 = pi/2;
    [~,fval] = patternsearch(fun_eig,theta0);
    val = -fval;
end
