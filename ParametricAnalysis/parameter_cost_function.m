function output = parameter_cost_function(input_tensor, mu_order, beta, mu_hat, m_index, d, x_lim)
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


[n_m, n_x] = size(mu_order)

% For normalisation purpose, we need to first compute the integral. 
output = tt_moment_full( input_t, zeros(1,n_x), m_index, d, x_lim, 0);

