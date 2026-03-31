% Parameters
lambda = 0.003;
delta_RA = 11;
K = 120;

% Jacobian Matrix
syms S RA RP

F1 = lambda*S*(1 - (S+RA+RP)/K);
F2 = lambda*(1+delta_RA)*RA*(1 - (S+RA+RP)/K);
F3 = lambda*RP*(1 - (S+RA+RP)/K);

sol = solve([F1==0, F2==0, F3==0], [S RA RP], 'Real', true);

disp('Equilibria:')
disp([sol.S sol.RA sol.RP])

% Define Jacobian matrix
vars = [S RA RP];
J = jacobian([F1 F2 F3], vars);

% Define equilibrium point as where all of ODEs are equal to zero, one such
% point being (0,0,0)
eq = [0 0 0];

% Evaluate the Jacobian at the chosen equilibrium point
J_eval = double(subs(J, vars, eq));

eigvals = eig(J_eval);

disp('Eigenvalues at (0,0,0):')
disp(eigvals)

% Create grid of possible values for S and RA between 0 and carrying capacity K = 120
[S, RA] = meshgrid(linspace(0,K,500), linspace(0,K,500));

% Compute RP from S + RA + RP = K for each (S, RA) pair
RP = K - S - RA;

% Define invalid region (negative RP)
mask = RP >= 0;

% Keep only points that satisfy S, RA, RP ≥ 0
S = S(mask);
RA = RA(mask);
RP = RP(mask);

% Plot all valid equilibrium points on the plane S + RA + RP = K
figure;
scatter3(S, RA, RP, 20, RP, 'filled')
hold on

xlabel('S')
ylabel('R_A')
zlabel('R_P')

title('Equilibrium Manifold: S + R_A + R_P = 120')
grid on

