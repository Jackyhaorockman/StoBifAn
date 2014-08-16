function Example_CFPE_Schlogl_X1_P4
% EXAMPLE_CFPE_SCHLOGL_X1_P4 Provides an template code for solving the
% steady state distribution of the parametric chemical Fokker-Planck 
% equation for the Schlogl model. All four biophysical parameters are
% treated as variables here, and the solution is saved in .sdv data
% structure.
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

% reaction rate constants (mean value)
rate = [0.18; 2.5e-4; 2250; 37.5];
rate_interval = [rate*0.97, rate*1.03]; % the interval is defined as the 3% variation around the mean value

% number of reactant molecules in each reaction
React = [2; 3; 0; 1];

% ==========


% ===== Define the simulation parameters ====

% lower and upper bounds for the truncated state space
x_lim = [0, 1000];

% number of grid points for each coordinate of the state space.
d = 10; % number of grid points: n = 2^d

% error tolerance for tensor rank truncation
tol_rank = 1e-12;

% total number of variable model parameter (here, we make all reaction rates as variables.)
p_var = 4;

% number of grid points for each coordinate of the parameter space.
pd = 7; % number of grid points: n = 2^pd

% specify the parameter to its reaction index:
ip = {{1},{2},{3},{4}};

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
file_name_time = 'cpu_time'; 
file_name_err = 'error_convergence';
file_name_dt = 'adaptive_time_step';
% ==========


% ====== Construct the Fokker-Planck operator =====

CFPE = Operator_CFPE_para(x_lim, d, v, [1 1 1 1], React, tol_rank, p_var, pd, rate_interval, ip, @CFPE_alpha_mass_action, 0);

RHS = Initial_CFPE_para(x_lim, d, task_RHS, mean, sigma, p_var, pd);

% ======


% ====== Solve for the steady state distribution and plot it =====

Adaptive_Shifted_IVP_save(CFPE, RHS, tol_IVP, tol_solve, file_name_tt, file_name_time, file_name_err, file_name_dt)

% ======

end