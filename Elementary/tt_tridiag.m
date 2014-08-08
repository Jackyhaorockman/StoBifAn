function operator = tt_tridiag(d,alpha, beta, gamma)
% TT_TRIDIAG  Construct a tri-diagonal matrix of size n^d x n^d, where the
% diagonal entries are all ALPHA's, the super-diagonal are all BETA's, and
% lower-diagonal are all GAMMA's.
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


tt = cell(d,1);
II = eye(2);
J = zeros(2);
J(1,2) = 1;

for key = 1 : d
    
    if (key == 1)
        
        tt{key} = zeros(2,2,3);
        tt{key}(:,:,1) = alpha*II + beta*J + gamma*J';
        tt{key}(:,:,2) = gamma*J;
        tt{key}(:,:,3) = beta*J';
        
    elseif (key == d)
        
        tt{key} = zeros(2,2,3);
        tt{key}(:,:,1) = II;
        tt{key}(:,:,2) = J';
        tt{key}(:,:,3) = J;
        
    else
        
        tt{key}=zeros(2,2,3,3);
        tt{key}(:,:,1,1) = II;
        tt{key}(:,:,2,2) = J;
        tt{key}(:,:,3,3) = J';
        tt{key}(:,:,2,1) = J';
        tt{key}(:,:,3,1) = J;
        
    end
    
end

operator = tt_matrix(tt);

end