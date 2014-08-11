function [tt_moment, tt_integral] = parameter_cost_function_prepare(input_t, mu_order, m_index, d, x_lim)
% PARAMETER_COST_FUNCTION_PREPARE Compute nessary tensor data for
% evaluating the cost function (see PARAMETER_COST_FUNCTION). 
%
% tt_moment, is a cell of length L, where L is the total number of moment
% orders to be considered. Each element, tt_moment{l}, represents a
% tensor-structured data for the l-th moment values for different parameter
% combinations.
%
% tt_integral, is a tensor-structured data which stores the sum over the
% state space for different parameter values. This will be used for
% normalisation purpose (see PARAMETER_COST_FUNCTION).
%
% The inputs required include:
%
% input_t, the original tensor-structured parametric solution.
%
% mu_order, an L x N matrix where L is the total number of moment orders
% and N is the total number of the chemical species. For the l-th row, the
% number on i-th colume corresponds to the moment order for the i-th
% chemical species.
%
% m_index, a vector of length N which specifies the index of the dimension
% which corresponds to the state space where moments are to be computed.
%
% d, is a vector of length N + M, specifies the number of tensor modes for 
% each dimension. (M is the total number of parameter values)
%
% x_lim, is a (N+M) x 2 matrix, where the N is the total number of the 
% molecular species. The first colume refers to the lower boundary and the 
% second row refers to the higher boundary. For the other M rows, the value
% will not be used.
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


[n_m, n_x] = size(mu_order);

tt_integral = tt_moment_full( input_t, zeros(n_x, 1), m_index, d, x_lim, 0);

mean = zeros(n_x, 1);

for i = 1 : n_x
    
    ei = zeros(1, n_x); ei(i) = 1;
    
    mean(i) = tt_moment_full(input_t, ei, m_index, d, x_lim, 0);
    
end

tt_moment = cell(n_m,1);

for i = 1 : n_m
    
    tt_moment{i} = tt_moment_full( input_t, mu_order(i,:), m_index, d, x_lim, mean);

end
    