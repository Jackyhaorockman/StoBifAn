function output = tt_1D_operator(A, I, n, d)
% TT_1D_OPERATOR  Extend a single dimensional matrix operator to N
% dimensions, where the d-th dimension is A, and all the others are
% identity matrix, I. For example:
%
% Operator = tt_1D_operator(A, I, 4, 2) is of the form:
%
% Operator = I \otimes A \otimes I \otimes I.
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

output = [];

for ii = 1 : n
    
    if any(ii == d)
        
        output = tkron(output,A);
        
    else
        
        output = tkron(output,I);
        
    end
    
end

end