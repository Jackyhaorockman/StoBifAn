function output = tt_merge(x)
% TT_MERGE   Construct the Kronecker product of N single dimensional matrix
% operators. The input, x, should be a cell of N tensor matrices.
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

output = [];

for ii = 1 : length(x)
    
    output = tkron(output,x{ii});
    
end

end
