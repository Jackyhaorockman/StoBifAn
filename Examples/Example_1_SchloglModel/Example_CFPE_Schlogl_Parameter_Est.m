function Example_CFPE_Schlogl_Parameter_Est
% EXAMPLE_CFPE_SCHLOGL_PARAMETER_EST Provide a template file for parameter
% estimation using the parameter continuation in tensor-structured data
% format. 
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

% load the tensor solution (see Example_CFPE_Schlogl_X1_P4)
dis = tt_tensor('parametric_steady_state');

% ----------- simulation parameters --------------------

var.beta = [1, 100, 0.001];   % weight coefficients for different moment orders.
var.beta = var.beta / sum(var.beta); % normalise the coefficients

var.mu_hat = [261.3168, 2.0295e4, -2.0414e5]; % emperical moments

mu_order = [1; 2; 3]; % moment orders to be considered.

m_index = [1, 1, 1];  % the index corresponds to the state space.

d = [10, 7, 7, 7, 7]; % number of modes for state and parameter spaces.

var.dp = [7 7 7 7];       % number of modes for dimensions in the parameter space.

x_lim = [0, 1000];    % lower and upper bound of the state space. 

var.dx = (x_lim(2) - x_lim(1)) / (2^d(1)-1); % size of grid cell

rate = [0.18; 2.5e-4; 2250; 37.5];  % mean value of the rate parameters

rate_interval = [rate*0.97, rate*1.03]; % the interval is defined as the 3% variation around the mean value

x_lim = [x_lim; rate_interval]; % combine interval of state space with parameter space.

[var.tt_moment, var.tt_integral] = parameter_cost_function_prepare(dis, mu_order, m_index, d, x_lim, var.mu_hat(1));

var.tol = (2.5/1000)^2; % tolerance for the method of moment

par_bound = [1, 2^var.dp(1);
            1, 2^var.dp(2);
            1, 2^var.dp(3)];

% ---------------------------------

% ------------ parametric continuation -------------------

para_match = parameter_continuation([64 64 64], par_bound, @check_fcn, var);

% --------------------------------

% ------- plot the combinations -----

figure; set(gca, 'fontsize', 15);

plot(para_match');

xlabel('parameter index');

ylabel('grid index')

% ------------------------



% ------------- check function ----------

    function check = check_fcn(pos, var)
        % this is a check function for estimating parameters, based on the
        % method of moment.
        
        % NOTE: here we only search for 3d parameter space, thus, we fix k4
        % at its mean (or true) value. 
        p_index = [pos(1), pos(2), pos(3), 64]; 
        
        cost_fcn = parameter_cost_function(var.tt_moment, var.tt_integral, var.beta, var.mu_hat, p_index, var.dp, var.dx);
        
        if cost_fcn <= var.tol
            
            check = 1;
            
        else
            
            check = 0;
            
        end
        
    end

% ------------------------

end