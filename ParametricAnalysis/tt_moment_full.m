function output = tt_moment_full( input_t, mu, m_index, d, x_lim, M)
% TT_MOMENT_FULL Computes the \mu_{i_1,i_2,...,i_N}-th order moment of the
% multi-dimensional probability distribution specified by input_t. The
% total number of dimension of state space is N, and the number is M for
% the parameter space
%
% mu, is a vector of length N, specifies the moment order to each
% dimension.
%
% m_index, is a vector of length N, specifies the index (or ranking) of the
% N state space within the N+M dimensional solution.
%
% d, is a vector of length N + M, specifies the number of tensor modes for 
% each dimension.
%
% x_lim, is a (N+M) x 2 matrix, where the N is the total number of the 
% molecular species. The first colume refers to the lower boundary and the 
% second row refers to the higher boundary. For the other M rows, the value
% will not be used.
%
% M, is a vector of length N, specifies the mean value if computing the
% higher order moment around the mean is required. To compute the raw
% moments, set M = 0.
%
%
% Note: if the input tensor is the parametric steady state distribution,
% then the output will be a tensor which contains the moment values for all
% possible combinations of the parameter values.
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


n_dim = length(mu);

for j = 1 : n_dim
    
    [input_t, d, x_lim] = tt_moment( input_t, mu(1), d, m_index(1), x_lim, M(1));
    
    mu(1) = [];
    
    m_index(1) = [];
    
    M(1) = [];
    
end

output = input_t;

end