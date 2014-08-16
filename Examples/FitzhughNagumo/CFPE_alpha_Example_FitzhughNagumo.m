function alpha = CFPE_alpha_Example_FitzhughNagumo(para)

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

% compute coefficients
% Reaction rates
V = para.rate(1);
theta = para.rate(2);
omega = para.rate(3);
gamma = para.rate(4);
phi = para.rate(5);


alpha = cell(num_react,1);

alpha{1} = 1/V *tt_merge({X{1}*(X{1}-I), I});

alpha{2} = 1/V^2 * tt_merge({X{1}*(X{1}-I)*(X{1}-2*I), I});

alpha{3} = theta * tt_merge({X{1}, I});

alpha{4} = 1/V * theta * tt_merge({X{1}*(X{1}-I), I});

alpha{5} = tt_merge({I, X{2}});

alpha{6} = omega * V * tt_merge({I, I});

alpha{7} = phi * tt_merge({X{1}, I});

alpha{8} = phi * gamma * tt_merge({I, X{2}});


end