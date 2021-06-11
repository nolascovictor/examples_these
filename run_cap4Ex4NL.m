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
m = 10; % Projection dimension
n = 3; % system size
R = 1; % radius of  Newton-Kantorovik verification
intvSol = [0,8]; % interval [t0, L], L > t0, of solution domain. Default case, that is, intvSol = [], is [0,1]
initCond = {6,6,0}; % u(t0) = initCond. 
finalCond = {}; % For IVP finalCond={}. For BVP finalCond = {u1(L), u2(L),...,un(L)}, put [] on emty spaces
chi = 10;              % Lorentz parameter;
epsilon = 8/3;         % Lorentz parameter;
gamma = 28;            % Lorentz parameter;
pmtrs = [chi epsilon gamma]; 

%% Define the equation and its derivatives 
f = @(t,u) [chi*(u(2) - u(1)), epsilon*u(1) - u(2) - u(1)*u(3), -gamma*u(3) + u(1)*u(2)]; % Lorentz's equation
d_uf = @(t,u) [-chi, chi, 0; epsilon - u(3), -1, -u(1); u(2), u(1), -gamma]; 

%% Compute numerical solution to u'(t) = f((L - t0)t + t0, u + x0)(L - t0), u(0) = 0
[b, t0, L, x0]   = compute_solution(f, m, n, intvSol, initCond, finalCond);
% b is a matrix whose elements are the coefficients of numerical solution w in H^{1,0}_B(I)^n
% [t0,L], L > t0, solution domain
% x0 = initCond

%% Verify the numerical solution
% [eta, nu, K, t_star, t_2star] = verify_solution(b, f, t0, L, x0, d_uf, R);
% eta is the norm of F(w)= w'(t) - f(t, w(t)) in L2(I)^n, where w is numerical solution
% nu is bijectivity modulus
% K is a D2_uf bound
% t_star \in [0,R] provides existence radius in (H^{1,0}_B)^n
% t_2star provides uniqueness radius in (H^{1,0}_B)^n.
ElapsedTime = toc;
disp(['Elapsed Time in seconds = ', num2str(ElapsedTime)])

%% Plot the solution
w = @(x) compute_u(b, (x-t0)/(L-t0)) + x0;
x = t0 + ((0:10)/10)*(L - t0);
A = zeros(size(x,2), n);
for i=1:size(x,2) 
    A(i,:) = w(x(i));
end
% 2D - Plot
p = plot(x, A(:,1), 'k', x, A(:,2), 'ro-', x, A(:,3), '--gs');xlim([t0,L]);
p(1).LineWidth = 3;
p(2).LineWidth = 2;
p(3).LineWidth = 1;
legend({'Função coordenada w_1', 'Função coordenada w_2', 'Função coordenada w_3'},'Location','northwest')
xlabel('Tempo t')
pos1 = get(gcf,'Position'); % get position of Figure(1) 
set(gcf,'Position', pos1 - [pos1(3)/2,0,0,0]) % Shift position of Figure(1) 
% 3D - Plot 
figure(3);
p1 = plot3(A(:,1),A(:,2), A(:,3), 'k');
p1.LineWidth = 3;
xlabel('Eixo X')
ylabel('Eixo Y')
zlabel('Eixo Z')
grid on
axis square
legend({'Solução numérica w'},'Location','northwest')