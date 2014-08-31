function Example_1_1_CFPE_Schlogl_X1_P0
% EXAMPLE_1_1_CFPE_SCHLOGL_X1_P0 Provides an template code for solving and
% plotting the steady state distribution of the Schlogl model with constant
% reaction rate.
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
v = [1; -1; 1; -1]; 

% reaction rate constants
rate = [0.18; 2.5e-4; 2250; 37.5];


% number of reactant molecules in each reaction
React = [2; 3; 0; 1];

% ==========


% ===== Define the simulation parameters ====

% lower and upper bounds for the truncated state space
x_lim = [0, 800];

% number of grid points for each coordinate of the state space.
d = 7; % number of grid points: n = 2^d

% error tolerance for tensor rank truncation
tol_rank = 1e-12;

% choose full matrix format [2] or tensor format [1]
task_Operator = 1;

% define the type of initial guess: Delta
task_RHS = 3;

% define the mean and variantion for the Gaussian distribution as the
% initial guess
mean = 500;
sigma = 0;

% define the tolerance for inverse power iterations
tol_IVP = 1e-3;

% define the error tolerance for tensor linear solver
tol_solve = 1e-5;
% ==========


% ====== Construct the Fokker-Planck operator =====

CFPE = Operator_CFPE(x_lim, d, v, rate, React, tol_rank, task_Operator);

RHS = Initial_CFPE(x_lim, d, task_RHS, mean, sigma);

% ======


% ====== Solve for the steady state distribution and plot it =====

Adaptive_Shifted_IVP_plot(CFPE, RHS, tol_IVP, tol_solve, x_lim, d, 1, 0);

% ======

end