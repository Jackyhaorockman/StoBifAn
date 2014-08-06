function [output, d] = tt_integration( tensor, d, int_D )
% TT_INTEGRATION Integrate out the int_D-th dimension from high-dimensional
% data, where the number of modes for each dimension are given in the
% vector, d.
%
% output, is the resulting tensor with int_D-th dimension being integrated.
%
% d (output), is the same vector as the input vector, d, except the
% int_D-th element is removed.
%
%
% ------------------------------
% TT-Toolbox 1.0, 2014
%
% This is StoBifAn Toolbox, written by Shuohao Liao
% Mathematical Institute, University of Oxford
% webpage: http://maths.ox.ac.uk/liao
%
% For all questions, bugs and suggestions please email
% liao@maths.ox.ac.uk
% -----------------------------


output = reshape(tensor, [2*ones(1,sum(d(1:int_D-1))), 2^d(int_D),...
    2*ones(1,sum(d(int_D+1:end)))]);

output = ttm(output, sum(d(1:int_D-1)) + 1, ones(2^d(int_D) ,1));

output = reshape(output, [2*ones(1,sum(d(1:int_D-1))), 2*ones(1,sum(d(int_D+1:end)))]);

d(int_D) = [];

end
