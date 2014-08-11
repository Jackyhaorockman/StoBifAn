function CME = Operator_CME(d_x, x_lim_low, k_rate, nu, react, tol)
% OPERATOR_CME_MASS_ACTION Attemps to construct the operator of the CME for
% the mass-action reaction system in tt format. The list of parameters
% are:
% 
% d_x, a row vector containing the number of grid points in each dimension.
%
% x_lim_low - a colomn vector containing the lower bound of the state
%       space for different species.
%
% k_rate - a row vector containing the reaction rate constants.
%
% nu - stochiometric matrix with row indicies for reactions and column
%       indicies for species.
%
% react - number of reactant molecules for each reaction;
%
% For example, the parameters for the simple system (O -> A -> B -> O) are:
%       d_x = [6, 6];
%       x_lim_low = [0; 0];
%       k_rate = [20, 1, .8];
%       nu = [ 1, 0;
%             -1, 1;
%             0, -1];
%       react = [0, 0;
%                1, 0;
%                0, 1];
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

%x_lim = [x_lim_low , [2^d_x(1)-1; 2^d_x(2)-1]];


% no. of chemical species
[num_react, num_spe]  = size(nu);

% vector of coordinates
x_coordinate = cell(num_spe,1);

for i = 1 : num_spe
    
    if x_lim_low(i) == 0
        
        x_coordinate{i} = tt_x(2*ones(1,d_x(i)));
        
    else
        
        x_coordinate{i} = x_lim_low(i)*tt_ones(2*ones(1,dx(i))) + tt_x(2*ones(1,dx(i)));
        
    end
    
end

% propensity function
alpha = cell(num_react,1);

for i = 1 : num_react
    
    alpha{i} = [];
    
    for j = 1 : num_spe
        
        if react(i,j) == 0
            
            temp_alph = tt_ones(2, d_x(j));
            
        else
            
            temp_alph = x_coordinate{j};
            
            for k = 2 : react(i,j)
                
                temp_alph = diag(temp_alph) * tt_matrix(tt_qshift(d_x(j),k-1)) * x_coordinate{j};
                
            end
            
        end
        
        alpha{i} = tkron(alpha{i}, temp_alph);
        
    end
    
    alpha{i} = alpha{i} * k_rate(i);
    
    alpha{i} = diag(alpha{i});
    
end

% Step operator
CME_step = cell(num_react,1);

I_full = tt_eye(2, sum(d_x));

for i = 1 : num_react
    
    CME_step{i} = [];
    
    temp_E = [];
    
    for j = 1 : num_spe
        
        temp_Ej = tt_matrix(tt_qshift(d_x(j),nu(i,j)));
        
        temp_E = tkron(temp_E, temp_Ej);
        
    end
    
    CME_step{i} = temp_E - I_full;
    
end

% CME operator
CME = [];

for i = 1 : num_react
    
    CME = CME + CME_step{i}*alpha{i};
    
end

CME = round(CME, tol);

end