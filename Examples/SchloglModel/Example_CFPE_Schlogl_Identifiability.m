function Example_CFPE_Schlogl_Identifiability
% EXAMPLE_CFPE_SCHLOGL_IDENTIFIABILITY Computes the cost function of the
% method of moments for the Schlogl model. The function compares the first
% three moment orders with the prescribed imperical data:
%
%       mean     = 261.3168;
%       variance = 2.0295e4;
%       skewness = -2.0414e5;
%
% The cost function is only computed for parameter pairs, i.e., k1-k3 or
% k2-k4, for the sake of plotting.
%
% For more information about this example, please see the related
% publication:
% 
% Liao, S., Vejchodsky, T. & Erban, R. (2014). Parameter estimation and 
% bifurcation analysis of stochastic models of gene regulatory networks: 
% tensor-structured methods. arXiv preprint arXiv:1406.7825.
%
%
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

% !!! NOTE: the simulation parameters here should be consistent with the
% parameters assigned to the function: Example_CFPE_Schlogl_X1_P4.

dis = tt_tensor('parametric_steady_state');

beta = [1 100 .001];   % weight coefficients for different moment orders.

mu_hat = [261.3168, 2.0295e4, -2.0414e5]; % emperical moments

mu_order = [1; 2; 3]; % moment orders to be considered.

m_index = [1, 1, 1];  % the index corresponds to the state space.

d = [10, 7, 7, 7, 7]; % number of modes for state and parameter spaces.

dp = [7 7 7 7];       % number of modes for dimensions in the parameter space.

x_lim = [0, 1000];    % lower and upper bound of the state space. 

dx = (x_lim(2) - x_lim(1)) / (2^d(1)-1); % size of grid cell

rate = [0.18; 2.5e-4; 2250; 37.5];  % mean value of the rate parameters

rate_interval = [rate*0.97, rate*1.03]; % the interval is defined as the 3% variation around the mean value

x_lim = [x_lim; rate_interval]; % combine interval of state space with parameter space.

[tt_moment, tt_integral] = parameter_cost_function_prepare(dis, mu_order, m_index, d, x_lim, mu_hat(1));


% extract the value of the cost function for each parameter combination in
% the 2D parameter space
cost = zeros(2^dp(1), 2^dp(3));

for i1 = 1 : 2^dp(1) % loop through the first parameter k1
    
    for i3 = 1 : 2^dp(3) % loop through the third parameter k3
        
        clc

        p_index = [i1, 2^dp(2)/2, i3, 2^dp(4)/2]
        
        cost(i1,i3) = parameter_cost_function(tt_moment, tt_integral, beta, mu_hat, p_index, dp, dx);

    end
    
end

% rescale the cost function for better visualisation
cost = log(cost + 1e-2); % rescale to log scale

cost = cost - min(cost(:)); % lower bound rescaled to 0

cost = cost / max(cost(:)); % upper bound rescaled to 1

% plot the rescaled function

figure; set(gca, 'fontsize', 20);

surf( linspace(rate_interval(3,1), rate_interval(3,2), 2^dp(3)), ...
    linspace(rate_interval(1,1), rate_interval(1,2), 2^dp(1)), ...
    cost,...
    'edgecolor', 'none'); 

axis tight; view(2);

ylabel('k_3'); xlabel('k_1');
        
end