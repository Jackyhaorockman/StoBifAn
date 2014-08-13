function output = parameter_cost_function(tt_moment, tt_integral, beta, mu_hat, p_index, d, dx)
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
% tt_moment and tt_integral, are the solutions obtained from the function
% TPA/ParametricAnalysis/parameter_cost_function_prepare.m
%
% beta, is a vector of length L, which defines the weight coefficients for
% each moment order.
%
% mu_hat, is a vector of length L, which defines the empirical moments.
%
% d, is a vector of length M, where M is the dimensionality of the
% parameter space. Thus, d(i) denotes the number of tensor modes for each
% dimension.
%
% dx, is the length (area/volume) of a single grid cell. This will be used
% for normalisation in the moment orders.
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
    
    ind = [ind, tt_ind2sub(2*ones(1,d(i)), p_index(i))];
    
end

% compute the cost function
J = 0;

for i = 1 : length(beta)
    
    mu = tt_moment{i}(ind) / tt_integral(ind) / dx;
    
    J = J + beta(i)*(mu - mu_hat(i))^2 / mu_hat(i)^2;
    
end

output = J;

