function output = tt_extract_1D(tensor, d1, d2, n)
% TT_EXTRACT_1D Extract the solution corresponding to a specific parameter
% value (for a single parameter). It is assume that the 1D parametric
% dimension of size d2 is added to the end, i.e., tensor = physical
% dimension + parametric dimension.
% 
% TENSOR, is the tensor-structured parametric probability distribution.
% 
% D1, is the number of the (virtual) physical dimensions, and also
% determines the size of the output solution.
%
% D2, is the number of dimensions for a single parameter. It will be
% integrated out.
%
% N, is the index of the parameter value, of which the tensor solution is
% extracted.
%
%

N = 2^d2;

sz = size(tensor);
sz(d1+1:end) = [];

lb = 0;
ub = N;

x = cell(d2,1);

for D = d2 :-1: 1
    if n > (lb + ub)/2
        x{D} = [0 1]';
        lb = (lb + ub)/2;
    else
        x{D} = [1 0]';
        ub = (lb + ub)/2;
    end
end

for i = 1 : d2
    tensor = ttm(tensor, d1+i, x{i});
end

output = reshape(tensor, sz);
        
end