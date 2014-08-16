function Example_CFPE_CellCycle_X4_p1
% EXAMPLE_CFPE_CELLCYCLE_X4_P1
%
%



% ====== Define the model parameters =====

% Stoichiometric matrix
v = [0,0,0,1,0;    % 0  --> X4
    0,0,0,-1,0;   % 0  <-- X4
    -1,1,0,-1,0;    % X1 + X4 --> X2
    0,-1,1,0,0;   % X2 --> X3
    0,1,-1,0,0;   % X2 <-- X3
    0,0,-1,0,1;   % X3 --> X5
    0,0,0,0,-1;   % X5 --> 0
    1,0,0,0,0;   % 0 --> X1
    -1,0,0,0,0;   % X1 --> 0
    -1,0,0,0,0;   % X1 --> 0
    ];            % v_ji denotes the change in the number of Xi produced by one Rj reaction.

% Reaction rates
para_alpha.CT = 5000;
para_alpha.aa = 1;
para_alpha.k4p = 0.018;
para_alpha.tP = .001;

rate = [0.015 * para_alpha.CT; 
    0;
    200 / para_alpha.CT;
	80; % 10 - 1000;
    0; 
	1; % bifurcation parameter
	0.6;
	1000;
    100];

rate_interval = [0.25, 0.4]; % for k_6

React = v; % this variable will NOT be used.

% ==========


% ===== Define the simulation parameters ====

% lower and upper bounds for the truncated state space
x_lim = [10, 70; 
        0, 1500; 
        0, 1200; 
        20, 70; 
        0, 700];
    
% number of grid points for each coordinate of the state space.
d = 8; % number of grid points: n = 2^d

% error tolerance for tensor rank truncation
tol_rank = 1e-12;

% total number of variable model parameter
p_var = 1;

% number of grid points for each coordinate of the parameter space.
pd = 5; % number of grid points: n = 2^pd

% specify the parameter to its reaction index:
ip = {{0},{0},{0},{0},{0},{1},{0},{0},{0},{0}};

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

CFPE = Operator_CFPE_para(x_lim, d, v, rate, React, tol_rank, p_var, pd, rate_interval, ip, @CFPE_alpha_mass_action, 0);

RHS = Initial_CFPE_para(x_lim, d, task_RHS, mean, sigma, p_var, pd);

% ======


% ====== Solve for the steady state distribution and plot it =====

Adaptive_Shifted_IVP_save(CFPE, RHS, tol_IVP, tol_solve, file_name_tt, file_name_time, file_name_err, file_name_dt)

% ======