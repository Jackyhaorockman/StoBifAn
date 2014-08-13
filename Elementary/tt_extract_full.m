function [tensor, d_new] = tt_extract_full(tensor, d, ind)
% TT_EXTRACT_FULL Extract the required data from the high-dimensional
% tensor-structured solutions, tensor, specified by the vector, ind. 
%
% d, is a vector of length N, which corresponds to the number of tensor
% modes for different dimensions (N dimensions in total).
%
% ind, is also a vector of length, which specifies the type of solutions to
% be extract from each dimension.
% For example, if ind(i) = 4, (or any value larger than 0), implies the
% solution at the 4-th grid points along the i-th dimension will be
% extracted. 
% If ind(i) = 0, the i-th dimensional will be integrated out.
% If ind(i) = -1 (or any value smaller than 0), the i-th dimensional data
% will be kept as its origianl tensor format.
%
% d_new, is a vector which specifies the number of tt modes for each
% dimensions in the output tensor. 
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



n_dim = length(d);

dm = 0;

d_new = [];

j = 0;

for i = 1 : n_dim
    
    if ind(i) == 0 % just integrate out this dimension
        
        ei = ones(2^d(i),1);
        
        tensor = reshape(tensor, [2*ones(1,dm), 2^d(i), 2*ones(1, sum(d(i+1 : end)))]);
        
        tensor = ttm(tensor, dm+1, ei);
        
        tensor = reshape(tensor, [2*ones(1,dm), 2*ones(1, sum(d(i+1 : end)))]);
        
    elseif ind(i) <= 0 % nothing to be done to this dimension
        
        dm = dm + d(i);
        
        j = j + 1;
        
        d_new(j) = d(i);
        
    else % extract the part of the solution
        
        ei = zeros(2^d(i), 1);
        
        ei(ind(i)) = 1;
        
        tensor = reshape(tensor, [2*ones(1,dm), 2^d(i), 2*ones(1, sum(d(i+1 : end)))]);
        
        tensor = ttm(tensor, dm+1, ei);
        
        tensor = reshape(tensor, [2*ones(1,dm), 2*ones(1, sum(d(i+1 : end)))]);
        
    end
   
    
    
end


end