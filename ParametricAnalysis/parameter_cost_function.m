function output = parameter_cost_function(tt_moment, tt_integral, beta, mu_hat, p_index, d, h)
% PARAMETER_COST_FUNCTION Compute the cost function for parameter
% estimation using the method of moments from the tensorised solution,
% input_tensor. The cost function will be a weighted sum of the difference
% between the model and empirical moments, i.e. (for 1D case),
%
%   J(k) = \sum_{l=1}^L beta_l * (\hat{mu}_l - mu_l)^2 / \hat{mu}^2.
%
% where L is the maximum order to be compared, \hat{mu}_l denotes the
% l-th order empirical moment and mu_l be the model moment. 
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


% compute the tt index for the location specified by p_index
ind = [];

for i = 1 : length(d)
    
    ind = [ind, tt_sub2ind(2*ones(1,d(i)), p_index(i))];
    
end

% compute the cost function
J = 0;

for i = 1 : length(beta)
    
    mu = tt_moment{i}(ind) / tt_integral(ind) / h;
    
    J = J + beta(i)*(mu - mu_hat(i))^2 / mu_hat(i)^2;
    
end

output = J;

