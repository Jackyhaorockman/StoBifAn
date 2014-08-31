function Example_2_1_1_CFPE_FitzhughNagumo_X2_p4
% EXAMPLE_2_1_1_CFPE_FITZHUGHNAGUMO_X2_P4 Present an example code for solving the
% steady state distribution of the Fitzhugh-Nagumo neuron model, which
% contains 2 chemical species (activator and inhibitor) and 4 model
% parameters. 
%
%
% For more information about this example, please see the related
% publication:
% 
% Liao, S., Vejchodsky, T. & Erban, R. (2014). Parameter estimation and 
% bifurcation analysis of stochastic models of gene regulatory networks: 
% tensor-structured methods. arXiv preprint arXiv:1406.7825.
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

% Stoichiometric matrix
v = [1, 0;
    -1, 0;
    -1, 0;
    1, 0;
    -1, 0;
    1, 0;
    0, 1;
    0, -1];

% Reaction rates
para_alpha.CT = 5000;
para_alpha.aa = 1;
para_alpha.k4p = 0.018;
para_alpha.tP = .001;

rate_ini = [2000 1 1 1 1];

% computational domain in parameter space
rate = [2000; 0.2; 0.112; 2.5; 0.105];
pcv = 15/100;
rate_interval = diag(rate) * ones(length(rate),1) * [1  - pcv, 1 + pcv];
rate_interval (1,:) = [];

React = v; % this variable will NOT be used.

% ==========


% ===== Define the simulation parameters ====

% lower and upper bounds for the truncated state space
x_lim = [-200, 1800;
        0, 700];
    
% number of grid points for each coordinate of the state space.
d = 8; % number of grid points: n = 2^d

% error tolerance for tensor rank truncation
tol_rank = 1e-12;

% total number of variable model parameter
p_var = 4;

% number of grid points for each coordinate of the parameter space.
pd = 7; % number of grid points: n = 2^pd

% specify the parameter to its reaction index:
ip = {{0};
    {0};
    {1};
    {1};
    {0};
    {2};
    {4};
    {3,4}};

% define the type of initial guess: uniform
task_RHS = 1;

% define the mean and variantion for the Gaussian distribution as the
% initial guess
mean = 0;
sigma = 0;

% define the tolerance for inverse power iterations
tol_IVP = 1e-3;

% define the error tolerance for tensor linear solver
tol_solve = 1e-5;

% define the path to save the data
file_name_tt = 'parametric_steady_state1';

file_name_time = 'cpu_time1'; 

file_name_err = 'error_convergence1';

file_name_dt = 'adaptive_time_step1';

% ==========


% ====== Construct the Fokker-Planck operator =====

CFPE = Operator_CFPE_para(x_lim, d, v, rate_ini, React, tol_rank, p_var, pd, rate_interval, ip, @CFPE_alpha_Example_FitzhughNagumo, 0);

RHS = Initial_CFPE_para(x_lim, d, task_RHS, mean, sigma, p_var, pd);

% ======


% ====== Solve for the steady state distribution and plot it =====

Adaptive_Shifted_IVP_save(CFPE, RHS, tol_IVP, tol_solve, file_name_tt, file_name_time, file_name_err, file_name_dt)

% ======