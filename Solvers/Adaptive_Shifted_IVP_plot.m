function Adaptive_Shifted_IVP_plot(Operator, rhs, err_tol, tt_tol, x_lim, d, d1, d2)
% ADAPTIVE_SHIFTED_IVP_SAVE  Implement adaptive shifted inverse power
% method for computing the steady state distribuiton of stochastic chemical
% systems.
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
% X_LIM is a N x 2 matrix, where the N is the total number of the molecular
% species. The first colume refers to the lower boundary and the second row
% refers to the higher boundary.
%
% D, is a vector of length N, of which the entries correspond to the number
% of modes in each dimension.
%
% D1, the index of the first of the two dimensions to be plotted.
%
% D2, the index of the second of the two dimensions to be plotted, if only
% single dimension is to be plotted, set d2 = 0.
%
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
rhs_pre = rhs;
err = norm(Operator * rhs);
if_up = 0;
if_down = 0;
if_repeat = 0;
num_iter = 0;
swp = 0;
figure;

I = tt_eye(2, length(size(Operator)));

% computation

tic
while err >= err_tol
    
    clc
    
    num_iter
    err
    dt
    swp
    if_repeat
    
    A = I - dt*Operator;
    
    % --- using amen_solve ----
    %      opts.nswp = 30;
    %      opts.kickrank=6;
    %      [rhs, ~, swp]=amen_solve(A, rhs, tt_tol, opts, rhs);
    % --- using amen_solve2 ----
    [rhs,testdata,~] = amen_solve2(A, rhs, tt_tol, 'nswp', 30, 'x0', rhs);
    swp = find(testdata{1}(1,:)==0,1,'first');
    if isempty(swp); swp = 30; end;
    % --------------------------
    
    % ------- change the time step -----
    if swp <= 4
        
        dt = dt*dt_scale;
        rhs = rhs_pre;
        if_plot = 0;
        if_up = 1;
        
    elseif swp >= 30
        
        dt = dt/dt_scale;
        rhs = rhs_pre;
        if_plot = 0;
        if_down = 1;
        
    elseif swp >= 20
        
        dt = dt/dt_scale;
        rhs = rhs / norm(rhs);
        err = norm(Operator * rhs);
        rhs_pre = rhs;
        if_plot = 1;
        if_up = 0;
        if_down = 0;
        if_repeat = 0;
        
    else
        
        rhs = rhs / norm(rhs);
        err = norm(Operator * rhs);
        rhs_pre = rhs;
        if_plot = 1;
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
    
    if if_plot == 1
        
        num_iter = num_iter + 1;
        
        if d2 == 0 % plot one dimensional distribution
            
            h = (x_lim(2) - x_lim(1)) / (2^d - 1);
            
            x_coordinate = linspace(x_lim(1), x_lim(2), 2^d);
            
            if length(d) > 1 % solution contains higher dimensional distribution
                
                x_dis = tt_extract_2d(rhs, d, d1, d2);
                
            else % solution contains only 1D distribution
                
                x_dis = rhs;
                
            end
            
            x_dis = full(x_dis); x_dis = x_dis / sum(x_dis); x_dis = abs(x_dis);
            
            x_dis = x_dis / h;
            
            plot(x_coordinate, x_dis, '-k', 'linewidth', 2); hold off;
            
            xlabel(sprintf('X_%d', d1));
            ylabel('Prob');
            
        else % plot two dimensional distribution
            
            h = (x_lim(:,2) - x_lim(:,1)) ./ (2.^d - 1);
            
            x1_coordinate = linspace(x_lim(d1,1), x_lim(d1,2), 2^d(d1));
            x2_coordinate = linspace(x_lim(d2,1), x_lim(d2,2), 2^d(d2));
            
            x_dis = tt_extract_2d(rhs, d, d1, d2);
            
            x_dis = x_dis / sum(x_dis(:)); 
            
            x_dis = x_dis / (h(d1) * h(d2));
            
            surf(x2_coordinate, x1_coordinate, x_dis, 'edgecolor', 'none');
            
            axis tight; colorbar; 
            
            xlabel(sprintf('X_%d',d2));
            ylabel(sprintf('X_%d',d1));
            
            hold off;
            
        end
        
        title(sprintf('eig = %1.6f', err));
        
        drawnow;
        
    end
    
    tic % re-start the time counting
end

end