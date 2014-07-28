function Adaptive_Shifted_IVP_save(Operator, rhs, err_tol, tt_tol, file_name_tt, file_name_time, file_name_err, file_name_dt)
% Adaptive shifted inverse power method for computing the steady state
% distribution of stochastic chemical systems.
%
% 
% Inputs:
%
% OPERATOR, is a tensor-structured operator (TT-matrix) which should
% contain a unique eigenvalue of smallest magnitude.
%
% RHS, is an initial guess for the ground state vector in tensor train
% format.
%
% ERR_TOL, defines the error tolerance for tensor solver.
%
% FILE_NAME_TT, is the director path and file name for saving tensor
% solutions.
%
% FILE_NAME_TIME, is the file name for saving the computational time.
%
% FILE_NAME_ERR, is the file name for saving the error for each inverse
% iterations.
%
% FILE_NAME_DT, is the file name for saving the adaptive time step for each
% inverse iteration (reciprical of the shift value).
%
% NOTE: there are more parameters in the tensor solver (alternative minimum
% energy method), amen_solve2, that can be tuned for better convergence,
% but change the default values only when you have a good understanding
% about what are you doing.
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



% computation parameters
dt = 1e-5;  % adaptive parameter: time step
dt_scale = 1.1;

rhs = rhs / norm(rhs);
err = norm(CME_para * rhs);
if_up = 0;
if_down = 0;
if_repeat = 0;
num_iter = 0;


I = tt_eye(2, length(size(Operator)));

% computation

tic
while err >= err_tol
    
    clc
    
    num_iter
    err
    dt
    if_repeat
    
    A = I - dt*Operator;
    
    % --- using amen_solve ----
    %      opts.nswp = 30;
    %      opts.kickrank=6;
    %      [rhs, ~, swp]=amen_solve(A, rhs, tt_tol, opts, rhs);
    % --- using amen_solve2 ----
    [rhs,testdata,~] = amen_solve2(A, rhs, tt_tol, 'nswp', 30, 'x0', rhs, 'ismex', 'true');
    swp = find(testdata{1}(1,:)==0,1,'first');
    if isempty(swp); swp = 30; end;
    % --------------------------
    
    % ------- change the time step -----
    if swp <= 3
        
        dt = dt*dt_scale;
        rhs = rhs_pre;
        if_save = 0;
        if_up = 1;
        
    elseif swp >= 30
        
        dt = dt/dt_scale;
        rhs = rhs_pre;
        if_save = 0;
        if_down = 1;
        
    elseif swp >= 20
        
        dt = dt/dt_scale;
        rhs = rhs / norm(rhs);
        err = norm(CME_para * rhs);
        rhs_pre = rhs;
        if_save = 1;
        if_up = 0;
        if_down = 0;
        if_repeat = 0;
        
    else
        
        rhs = rhs / norm(rhs);
        err = norm(CME_para * rhs);
        rhs_pre = rhs;
        if_save = 1;
        if_up = 0;
        if_down = 0;
        if_repeat = 0;
        
    end
    % -----------------
    
    % --- change time-step scaling when go back & forth ----
    
    if (if_up == 1) && (if_down == 1)
        
        if_repeat = if_repeat + 1;
        if_up = 0;
        if_down = 0;
        
    end
    
    if if_repeat >= 3
        
        dt_scale = mod(dt_scale,1) + 1.1;
        if_repeat = 0;
        
    end
    
    % --- save the data ----------------
    
    if if_save == 1
        
        num_iter = num_iter + 1;
        
        % save the data
        
        tt_save(file_name_tt, rhs, 0);
        save_time(num_iter) = toc;
        save(file_name_time, 'save_time');
        save_err(num_iter) = err;
        save(file_name_err, 'save_err');
        save_dt(num_iter) = dt;
        save(file_name_dt, 'save_dt');
        
    end
    
    tic % re-start the time counting
end

end