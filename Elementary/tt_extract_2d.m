function tensor = tt_extract_2d(tensor, d, d1, d2)
% TT_EXTRACT_2D Extract the 2D solution along d1-th and d2-th dimensions
% from the tensor-structured solution. d, is a vector whose entries record
% the number of modes for each dimension. To return single dimensional
% result, let d2 = 0.
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


D = length( d );

if d2 == 0
    n = d;
    check = zeros(1,D); check(d1) = 1;
    
    while any( check == 0 )
        dim = find(check == 0);
        [tensor, n] = tt_integration(tensor, n, dim(1));
        check(dim(1)) = [];
    end
    tensor = full(tensor);
    tensor = abs(tensor);
    
else
    
    n = d;
    check = zeros(1,D); check(d1) = 1; check(d2) = 1;
    
    while any( check == 0 )
        dim = find(check == 0);
        [tensor, n] = tt_integration(tensor, n, dim(1));
        check(dim(1)) = [];
        dim(1) = [];
    end
    tensor = full(tensor);
    tensor = abs(tensor);
    tensor = reshape(tensor,2^d(d1),2^d(d2));
end
end