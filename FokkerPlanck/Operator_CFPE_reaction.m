function CFPE = Operator_CFPE_reaction(x_lim, d, v, rate, React, tol_rank)
% OPERATOR_CPFE_REACTION  Compute Reaction-Wise Finite Difference Operaotor
% for the Chemical Fokker-Planck equation in tensor train format
% (MASS-ACTION reaction only).
%
% CPFE = OPERATOR_CFPE_REACTION(X_LIM, D, V, RATE, REACT, TOL_RANK)
%
% CFPE, is a cell of length M, where M is the total number of reaction
% channels. The j-th element, CFPE{j}, express the part of the finite
% difference operator of the CFPE contributed by the j-th reaction.
%
%
% X_LIM is a N x 2 matrix, where the N is the total number of the molecular
% species. The first colume refers to the lower boundary and the second row
% refers to the higher boundary.
%
% D determines the number of grid nodes in each dimension, i.e., n = 2^D.
%
% V is the stoichiometric matrix, where the V_ji refers to the change of
% the i-th molecular population in the j-th reaction.
%
% RATE is a vector of reaction constants.
%
% REACT is a matrix of the same size as the stiochiometric matrix, but
% REACT_ji denotes the number of i-th substrate consumned in the j-th
% reaction.
%
% TOL_RANK is a user-defined coefficient for the level of tensor truncation
% errors.
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


n = 2^d;
h = (x_lim(:,2) - x_lim(:,1)) / n;

[num_react, num_spe]  = size(v);

Lap = - tt_qlaplace_dd(d);
Cdiff = tt_tridiag(d, 0, 1/2, -1/2);
I = tt_eye(2,d);

% Construct diagonal matrices for each dimension, where the diagonal
% elements corresponds to the grid coordinates.
X = cell(num_spe,1);

for i = 1 : num_spe
    
    X{i} = diag(x_lim(i,1)*tt_ones(2*ones(1,d)) + h(i)*tt_x(2*ones(1,d)));
    
end

% Construct diagonal matrices for each dimension, where the diagonal
% elements corresponds to the grid coordinates.
X = cell(num_spe,1);

for i = 1 : num_spe
    
    X{i} = diag(x_lim(i,1)*tt_ones(2*ones(1,d)) + h(i)*tt_x(2*ones(1,d)));
    
end

% Compute the propensity function for each reaction.
alpha = cell(num_react,1);

for i = 1 : num_react
    
    alpha_i = cell(num_spe,1);
    
    for j = 1 : num_spe
        
        if React(i,j) == 0
            
            alpha_i{j} = I;
            
        else
            
            alpha_i{j} = X{j};
            
            if React(i,j) > 1
                
                for k = 1 : React(i,j)-1
                    
                    alpha_i{j} = alpha_i{j} * (X{j} - k*I);
                    
                end
                
                alpha_i{j} = round(alpha_i{j},tol_rank);
                
            end
            
        end
        
    end
    
    alpha{i} = tt_merge(alpha_i);
    
end


for i = 1 : num_react
    
    alpha{i} = rate(i) * alpha{i};
    
end


% Compute the 1st, 2nd and mixed derivatives in FDM
Dx = cell(num_spe,1);
Dxx = cell(num_spe,1);
Dxy = cell(num_spe,num_spe+1);

for i = 1 : num_spe
    
    Dx{i} = (1/(h(i))) * tt_1D_operator(Cdiff, I, num_spe, i);
    Dxx{i} = (1/h(i)^2) * tt_1D_operator(Lap, I, num_spe, i);
    
    if i < num_spe
        
        for j = i+1:num_spe
            
            Dxy{i,j} = (1/(h(i))) * (1/(h(j))) * tt_1D_operator(Cdiff, I, num_spe, [i,j]);
            
        end
        
    end
    
end

% Compute the elementary part of the Fokker-Planck operator that
% contributed from each separate reaction.
D = cell(num_spe,num_react);
DD = cell(num_spe,num_react);

for j = 1 : num_react
    
    for i = 1 : num_spe
        
        if v(j,i) ~= 0
            
            D{i,j} = D{i,j} + v(j,i)*alpha{j};
            D{i,j} = round(D{i,j}, tol_rank);
            
            DD{i,j} = DD{i,j} + v(j,i)^2*alpha{j};
            DD{i,j} = round(DD{i,j}, tol_rank);
            
        end
        
    end
    
    Diffusion = [];
    Drift = [];
    
    for i = 1 : num_spe
        
        if ~isempty( DD{i,j} )
            
            Diffusion = Diffusion + 0.5*Dxx{i}*DD{i,j};
            Diffusion = round(Diffusion, tol_rank);
            
        end
        
        if ~isempty( D{i,j} )
            
            Drift = Drift - Dx{i}*D{i,j};
            Drift = round(Drift, tol_rank);
            
        end
        
    end
    
    CFPE{j} = round(Diffusion + Drift, tol_rank);
    
end


end