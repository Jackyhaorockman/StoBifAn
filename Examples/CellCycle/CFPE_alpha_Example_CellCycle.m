function alpha = CFPE_alpha_Example_CellCycle(para)
% CFPE_ALPHA_EXAMPLE_CellCycle Returns the propensity functions for
% all reactions involved in the cell cycle model in tensor format. 
%
% see Example_CFPE_CellCycle_X5_p1.m
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

alpha = cell(num_react,1);

% compute coefficients
P3 = X{3}*X{3}*para.rate(4)/para.CT^2 + para.k4p*I;

alpha{1} = para.rate(1) * para.aa * tt_merge({I, I, I, I, I});

alpha{2} = para.rate(2) * tt_merge({I, I, I, X{4}, I});

alpha{3} = para.rate(3) * tt_merge({X{1}, I, I, X{4}, I});

alpha{4} =      tt_merge({I, X{2}, P3, I, I});

alpha{5} = para.rate(5) * para.tP * tt_merge({I, I, X{3}, I, I});

alpha{6} = para.rate(6) * tt_merge({I, I, X{3}, I, I});

alpha{7} = para.rate(7) * tt_merge({I, I, I, I, X{5}});

alpha{8} = para.rate(8) * tP * CT * tt_merge({I, I, I, I, I});

alpha{9} = para.rate(8) * para.tP * (tt_merge({X{1}, I, I, I, I}) + tt_merge({I, X{2}, I, I, I}) + tt_merge({I, I, X{3}, I, I}));

alpha{10} = para.rate(9) * tt_merge({X{1}, I, I, I, I});


end