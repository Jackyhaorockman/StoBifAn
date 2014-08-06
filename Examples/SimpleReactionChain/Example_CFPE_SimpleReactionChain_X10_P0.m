function Example_CFPE_SimpleReactionChain_X10_P0
% EXAMPLE_REACTIONCHAIN_X1_P0 Solves for the steady state distribution for
% the reaction chain with ten chemical species with constant reaction
% rates.
%
% 0 --> X1 --> X2 --> ... --> X10 --> 0
%
%
% ------------------------------
% TT-Toolbox 1.0, 2014
%
% This is StoBifAn Toolbox, written by Shuohao Liao
% Mathematical Institute, University of Oxford
% webpage: http://maths.ox.ac.uk/liao
%
% For all questions, bugs and suggestions please email
% liao@maths.ox.ac.uk
% -----------------------------


% ====== Define the model parameters =====

% stochiometric matrix
v = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    -1, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, -1, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, -1, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -1, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, -1, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, -1, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, -1, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, -1, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 0, -1, 1;
    0, 0, 0, 0, 0, 0, 0, 0, 0, -1]

% reaction rate constants
rate = [200, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

% number of reactant molecules in each reaction
React = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1];

% ==========


% ===== Define the simulation parameters ====

% lower and upper bounds for the truncated state space
x_lim = [100, 300;
        100, 300;
        100, 300;
        100, 300;
        100, 300;
        100, 300;
        100, 300;
        100, 300;
        100, 300;
        100, 300];

% number of grid points for each coordinate of the state space.
d = 6; % number of grid points: n = 2^d

% error tolerance for tensor rank truncation
tol_rank = 1e-12;

% choose full matrix format [2] or tensor format [1]
task_Operator = 1;

% define the type of initial guess: Uniform
task_RHS = 1;

% define the mean and variantion for the Gaussian distribution as the
% initial guess
mean = 0;
sigma = 0;

% define the tolerance for inverse power iterations
tol_IVP = 1e-2;

% define the error tolerance for tensor linear solver
tol_solve = 1e-3;
% ==========


% ====== Construct the Fokker-Planck operator =====

CFPE = Operator_CFPE(x_lim, d, v, rate, React, tol_rank, task_Operator);

RHS = Initial_CFPE(x_lim, d*ones(1,10), task_RHS, mean, sigma);

% ======


% ====== Solve for the steady state distribution and plot it =====

Adaptive_Shifted_IVP_plot(CFPE, RHS, tol_IVP, tol_solve, x_lim, d*ones(10,1), 3, 7);

% ======

end