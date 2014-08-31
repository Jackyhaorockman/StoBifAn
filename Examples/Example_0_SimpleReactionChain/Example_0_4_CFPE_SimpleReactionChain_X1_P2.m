function Example_0_4_CFPE_SimpleReactionChain_X1_P2
% EXAMPLE_REACTIONCHAIN_X1_P2 Solves for the parametric steady state
% distribution for the simple death-birth process (single species) with
% varying birth and death rates.
%
%
% ------------------------------
% StoBifan 1.0, 2014
%
% This is Stochastic Bifurcation Analyser, written by Shuohao Liao
% Mathematical Institute, University of Oxford
% webpage: http://maths.ox.ac.uk/liao
%
% For all questions, bugs and suggestions please email
% liao@maths.ox.ac.uk
% -----------------------------



% ====== Define the model parameters =====

% stochiometric matrix
v = [1; -1]; 

% initial reaction rate constants
rate = [1, 1]; 

% number of reactant molecules in each reaction
React = [0; 1];

% ==========


% ===== Define the simulation parameters ====

% lower and upper bounds for the truncated state space
x_lim = [0, 200];

% total number of variable model parameter (birth and death rates)
p_var = 2;

% lower and uppper bounds for the birth and death rates
rate_interval = [50, 150; % birth
                0.8, 1.5]; % death

% number of grid points for each coordinate of the state space.
d = 6; % number of grid points: n = 2^d

% number of grid points for each coordinate of the parameter space.
pd = 4; % number of grid points: n = 2^pd

% specify the parameter to its reaction index:
ip = {{1},{2}}; % the 1st reaction involves the 1st parameter, and the 2nd reaction involves the 2nd parameter.

% choose the solution to be plotted. 
% Note: the number should not exceed 2^pd!
ind_plot = [10, 10]; % the 100th grid points for birth rate; and same for the death rate.

% error tolerance for tensor rank truncation
tol_rank = 1e-12;

% define the tolerance for inverse power iterations
tol_IVP = 1e-1;

% define the error tolerance for tensor linear solver
tol_solve = 1e-4;

% ==========


% ====== Construct the Fokker-Planck operator =====
CFPE = Operator_CFPE_para(x_lim, d, v, rate, React, tol_rank, p_var, pd, rate_interval, ip, @CFPE_alpha_mass_action, 0);

RHS = Initial_CFPE_para(x_lim, d, 1, 0, 0, p_var, pd); % uniform distribution as initial guess

% ======


% ====== Solve for the steady state distribution =====

sol = Adaptive_Shifted_IVP(CFPE, RHS, tol_IVP, tol_solve);

% ======

% ====== plot the solution for different parameter values


[sol, ~] = tt_extract_full(sol, [d, pd, pd], [-1, ind_plot]);

sol = full(sol); sol = abs(sol); sol = sol / sum(sol) / ((x_lim(2) - x_lim(1))/(2^d-1));

figure; set(gca, 'fontsize', 15);

plot(linspace(x_lim(1), x_lim(2), 2^d), sol);

xlabel('Number'); 

ylabel('Prob')

% =======