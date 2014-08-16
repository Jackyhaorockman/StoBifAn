function alpha = CFPE_alpha_mass_action(para)
% CFPE_ALPHA_MASS_ACTION Return a cell array of length M, which contains
% the propensity functions for each reaction of any mass-action reaction
% system.
%
% The input parameter includes:
%
% para.d_x, is the number of tt modes for each dimension. 
%
% para.x_lim, is an N x 2 matrix, storing the lower and upper bound of the
% state space.
%
% para.v, is the stochiometric matrix.
%
% para.React, is the matrix denotes the number of substrate molecules in
% each reaction.
% 
% para.rate, is a vector of length M, defining the reaction rate constants
% for all mass-action reactions.
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

n = 2^para.d_x;

h = (para.x_lim(:,2) - para.x_lim(:,1)) / n;

[num_react, num_spe]  = size(para.v);

I = tt_eye(2,para.d_x);


% Construct diagonal matrices for each dimension, where the diagonal
% elements corresponds to the grid coordinates.
X = cell(num_spe,1);

for i = 1 : num_spe
    
    X{i} = diag(para.x_lim(i,1)*tt_ones(2*ones(1,para.d_x)) + h(i)*tt_x(2*ones(1,para.d_x)));
    
end

% Compute the propensity function for each reaction.
alpha = cell(num_react,1);

for i = 1 : num_react
    
    alpha_i = cell(num_spe,1);
    
    for j = 1 : num_spe
        
        if para.React(i,j) == 0
            
            alpha_i{j} = I;
            
        else
            
            alpha_i{j} = X{j};
            
            if para.React(i,j) > 1
                
                for k = 1 : para.React(i,j)-1
                    
                    alpha_i{j} = alpha_i{j} * (X{j} - k*I);
                    
                end
                
                alpha_i{j} = round(alpha_i{j},para.tol_rank);
                
            end
            
        end
        
    end
    
    alpha{i} = tt_merge(alpha_i);
    
end


for i = 1 : num_react
    
    alpha{i} = para.rate(i) * alpha{i};
    
end

end