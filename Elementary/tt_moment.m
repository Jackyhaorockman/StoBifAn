function [output, d] = tt_moment( input_t, mu, d, km, x_lim, M)
% TT_MOMENT Computes the mu-th order moment in the km-th dimension from the
% tensorised data provided by input_t. Here, d is a vector containing the
% number of modes for each dimension. x_lim is an N x 2 matrix storing the
% lower and upper bound for each coordinate. For mu > 1, exact value for
% the mean, M, is required, and to compute the raw moment, let M = 0. For
% mu = 1, M value will not be used.
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
% compute the mean variance for int_D dimension


% reshape the quantised tt representation to tt represenation for the
% dimension to be integrated
output = reshape(input_t, [2*ones(1,sum(d(1:km-1))), 2^d(km),...
    2*ones(1,sum(d(km+1:end)))]);

if mu == 1
    
    M = 0;
    
end

% grid value
x_data = linspace(x_lim(km,1), x_lim(km,2), 2^d(km));

% mode multiplication to the mom_D-th dimension
output = ttm(input_t, sum(d(1:km-1)) + 1, (x_data - M).^mu);

d(km) = [];

% to squeeze out the dimension being integrated
output = reshape(output, 2*ones(1,sum(d)));

end