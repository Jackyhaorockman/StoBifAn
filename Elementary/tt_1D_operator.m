function output = tt_1D_operator(A, I, n, d)
% Extend a single dimensional matrix operator to n dimensions, where the
% d-th dimension is A, and all the others are the identity matrix, I.
%
%
% ------------------------------
% TT-Toolbox 1.0, 2014
%
% This is StoBifAn Toolbox, written by Shuohao Liao
% Mathematical institute, University of Oxford
% webpage: http://maths.ox.ac.uk/liao
%
% For all questions, bugs and suggestions please mail
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