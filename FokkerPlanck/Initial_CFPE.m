function output = Initial_CFPE(x_lim, d, task, mean, sigma)
% INITIAL_CFPE Returns initial guess vector for solving CFPE equations in
% tensor train format. Four types of initial conditions are supported
% assigned by the value of [task].
%
% For [task = 1], initial guess is assumed as uniform distribution over the
% state space, i.e., all the values over the grid points equal to 1.
% 
% For [task = 2], initial guess is assumed randomly distribution with zero
% mean and unit variance.
% 
% For [task = 3], initial guess is assumed as delta distribution at the
% grid points closest to the point specified by [mean].
%
% For [task = 4], initial guess is assumed to be Gaussian distributed, with
% mean and variance(sigma) specified.
%
% Required parameters are:
%
% X_LIM is a N x 2 matrix, where the N is the total number of the molecular
% species. The first colume refers to the lower boundary and the second row
% refers to the higher boundary.
%
% D determines the number of grid nodes in each dimension, i.e., n = 2^D.
%
% TASK as discssed above.
%
% MEAN, is a N x 1 vector of the location of the mean population for each
% separate chemical species.
%
% SIGMA, is a N x 1 vector of the standard variance of the population
% around its mean.
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


switch task
    case 1 % return an initial guess of all ones over all tensor grid points
        
        output = tt_ones(2, sum(d));
        
    case 2 % return an initial guess of random numbers
        
        output = tt_rand(2, sum(d), 10); 
        % NOTE: here 10 is the default value for the rank bound, increase
        % it if you prefer more randomness or decrease it if you want to
        % save storage.
        
    case 3 % return an initial guess of delta distributio at the mean value
        
        num_spe = length(d);
        
        output = [];
        
        for i = 1 : num_spe
            
            % grid size on the i-th dimension
            h = (x_lim(i,2) - x_lim(i,1)) / (2^d(i) - 1);
            
            % index of the grid points which support the delta
            % distribution (rounding to the nearest grid nodes if
            % necessary).
            dlt = round((mean(i) - x_lim(i,1)) / h + 1);
            
            % TT representation for the unit vector with 1 at the dlt-th
            % points
            x = tt_unit(2, d(i), dlt);
            
            output = tkron(output, x);
            
        end
        
    case 4 % return an initial guess as the normal distribution.
        
        num_spe = length(d);
        
        output = [];
        
        for i = 1 : num_spe
            
            % vector containing the node points in the i-th dimension
            x_full = linspace(x_lim(i,1), x_lim(i,2), 2^d(i));
            
            % probability density function over the grid points
            x_prob = normpdf(x_full, mean(i), sigma(i));
            
            % reshape into TT format (default rounding value: 1e-16)
            x = tt_qformfull(x_prob, 2, d(i), 1e-16);
            x = tt_tensor(x);
            
            output = tkron(output, x);
            
        end
end

end