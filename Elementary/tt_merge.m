function output = tt_merge(x)
% TT_MERGE   Construct the Kronecker product of N single dimensional matrix
% operators. The input, x, should be a cell of N tensor matrices.
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

for ii = 1 : length(x)
    
    output = tkron(output,x{ii});
    
end

end
