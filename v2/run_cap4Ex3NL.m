% It's provide a rigorous verification of numerical solution of IVP or BVP
% u'(t) - f(t,u(t)) = [0,...,0], u(t0) = x0 (u(L) = x1) where f: R^(n+1) \to R^n, f=f(t,u(t)) is
% C2 and u: [t0, L] \to Rn, u = u(t) = (u1(t),...,un(t)). That is, if there
% exists a true solution near to numerical solution obtained. For that,
% the first part, "compute_solution", compute a numerical solution "b" for the system 
% u'(t)- f((L - t0)t + t0, u(t) + x0)(L - t0) = [0,...,0] 
% subject to initial condition u(0) = [0,...,0] (u(1) = x1-x0) by fsolver and ode45. 
% Note that, "b" has the coefficients of numerical 
% solution "w" in (H^{1,0}_B)^n. After that, using
% Theorem 2 of cited paper above, "verify_solution"
% verifies if there exists a true solution of u'(t)- f((L - t0)t + t0, u(t) + x0)(L - t0) = [0,...,0] 
% u(0) = [0,...,0] (u(1) = x1-x0), in closed ball B[w,R] 
% of (H^{1,0}_B)^n. The "t_star \in [0,R]" provides existence radius and "t_2star"
% provides uniqueness radius in (H^{1,0}_B)^n.
% 
% Input: f, d_uf, m, n, R, parameters (if necessary), intvSol = [t0, L], initCond and finalCond (for BVP).
% Output: Message about verification state and plot of numerical solution. If the message is 'Proof was successful!' then 
% it shows eta, nu, K, t_star, t_2star, elapsed time in seconds.

set(0, 'DefaultAxesFontSize', 18)
set(0, 'DefaultAxesFontWeight', 'bold')
tic
clear; clc; close all; close all force;

%% Define initial guess and parameters
m = 42; % Projection dimension
n = 2; % system size
R = 1; % radius of  Newton-Kantorovik verification
intvSol = [0, 2]; % interval [t0, L], L > t0, of solution domain. 
initCond = {2, 0}; % u(t0) = initCond. 
finalCond = {}; % For IVP finalCond={}. For BVP finalCond = {u1(L), u2(L),...,un(L)}, put [] on emty spaces
c = 1; % Van der Pol parameter

%% Define the equation and its derivatives 
f = @(t,x) [x(2), c*(1 - x(1)^2)*x(2) - x(1)];  % Van der Pol equation
d_uf = @(t,x)[0, -2*c*x(1)*x(2)-1; 1, c-c*x(1)^2];

%% Compute numerical solution to u'(t) = f((L - t0)t + t0, u + x0)(L - t0), u(0) = 0
[b, t0, L, x0]   = compute_solution(f, d_uf, m, n, intvSol, initCond, finalCond);
% b is a matrix whose elements are the coefficients of numerical solution w in H^{1,0}_B(I)^n
% [t0,L], L > t0, solution domain
% x0 = initCond


% %% Calcula nu0 e F_Bnm para teste
% new_f = @(t,x) f((L - t0)*t + t0, x + x0).*(L - t0);
% new_d_uf = @(t,x) (L-t0).*d_uf((L - t0)*t + t0, x + x0);
% 
% % nu0 = compute_nu0T1(new_d_uf, b)
% F_Bnm = @(a) integral(@(t) ([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]')*([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]*a - new_f(t, compute_u(a, t))), 0, 1,'ArrayValued',true);
% norm(F_Bnm(b))
% tic
% [convergiu, solaprox, num_itr] = Newton(new_f, new_d_uf, b, 1e-6, 100);
% convergiu
% num_itr
% b = solaprox;
% toc
% norm(F_Bnm(b))
% x = reshape(b,[1,m*n]);
% Eta = sqrt(sum(x.^2) - 2*(x(1:m:end)*integral(@(t) new_f(t, compute_u(reshape(x, [m,n]), t)), 0, 1,'ArrayValued',true)') + sum(integral(@(t) new_f(t, compute_u(reshape(x, [m,n]), t)).^2, 0, 1,'ArrayValued',true)) - 2*sqrt(2)*sum(integral(@(t) (((reshape(x, [m,n])')*([1 cos((1:(m-1)).*(pi*t))]'))')*(new_f(t, compute_u(reshape(x, [m,n]), t))'), 0, 1,'ArrayValued',true)))
% 
% ETA = sqrt(sum(integral(@(t) ([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]*F_Bnm(b)).^2, 0, 1,'ArrayValued',true)))
% 
% 
% w = @(x) compute_u(b, x); 
% dw = @(x) compute_du(b, x);
% F_u_2 = @(x) (dw(x) - new_f(x,w(x))).^2;
% eta3 = sqrt(sum(integral(F_u_2, 0, 1,'ArrayValued',true)))%(L2)^n norm
% 
% 
% F_u = @(x) compute_du(b, x) - new_f(x, compute_u(b, x));
% F_a = zeros(m,n);
% F_a(1,:) = integral(F_u, 0, 1,'ArrayValued',true);
% for i = 2:m
%     F_a(i,:) = integral(@(x) F_u(x).*sqrt(2).*cos((i-1).*pi.*x), 0, 1, 'ArrayValued',true);
% end
% 
% ETA4 = sqrt(sum(integral(@(t) ([1 sqrt(2).*cos((1:(m-1)).*pi.*t)]*F_a).^2, 0, 1,'ArrayValued',true)))
% 
% 

%% Verify the numerical solution
[eta, nu, K, t_star, t_2star] = verify_solution(b, f, t0, L, x0, d_uf, R);
% eta is the norm of F(w)= w'(t) - f(t, w(t)) in L2(I)^n, where w is numerical solution
% nu is bijectivity modulus
% K is a D2_uf bound
% t_star \in [0,R] provides existence radius in (H^{1,0}_B)^n
% t_2star provides uniqueness radius in (H^{1,0}_B)^n.
ElapsedTime = toc;
disp(['Elapsed Time in seconds = ', num2str(ElapsedTime)])

%% Plot the solution
w = @(x) compute_u(b, (x-t0)/(L-t0)) + x0;
x = t0 + ((0:100)/100)*(L - t0);
A = zeros(size(x,2), n); 
for i=1:size(x,2) 
    A(i,:) = w(x(i));
end
% 2D - Plot I
figure(1);
p = plot(x, A(:,1), 'k', x, A(:,2), '--gs'); xlim([t0,L]);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
legend({'Função coordenada w_1', 'Função coordenada w_2'},'Location','northwest')
xlabel('Tempo t')
pos1 = get(gcf,'Position'); % get position of Figure(1) 
set(gcf,'Position', pos1 - [pos1(3)/2,0,0,0]) % Shift position of Figure(1) 
% 3D - Plot 
figure(3);
p1 = plot3(A(:,1),A(:,2), x, 'k');
p1.LineWidth = 3;
xlabel('Eixo X')
ylabel('Eixo Y')
zlabel('Tempo t')
grid on
axis square
legend({'Solução numérica w'},'Location','northwest')