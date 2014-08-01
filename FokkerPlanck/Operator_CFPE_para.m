function output = Operator_CFPE_para(x_lim, d, v, rate_ini, React, tol_rank, p_var, pd, plim, ip)
% OPEARTOR_CFPE_PARA  Construct the parametric finite difference operator
% for the chemical Fokker-Planck equation in tensor train format
% (MASS-ACTION REACTION ONLY), as a generalised form of the sylver matrix:
% 
% OUTPUT = A_0 * I * ... * I + A_1 * K_1 * ... * I + A_M * I * ... * K_M
% 
% where A_0 is the non-parametric part; A_j is the parametric part for the
% j-th (varying) parameter; K_j is the diagonal matrix whose diagonal
% entries are discrete sample value of the j-th bio-physical parameters;
% and I is the identity matrix.
%
%
% The following inputs are required:
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
% RATE_INI, is a vector of reaction constants, which will be used to
% contruct the elementary operators, A_0, A_1, ..., A_M, as described
% above.
%
% REACT, is a matrix of the same size as the stiochiometric matrix, but
% REACT_ji denotes the number of i-th substrate consumned in the j-th
% reaction.
%
% TOL_RANK is a user-defined coefficient for the level of tensor truncation
% errors.
%
% P_VAR, is the total number of (varying) parameters, i.e., the
% dimensionality of the parameter space.
%
% PD, determines the number of discrete sample points for each parameter is
% 2^PD.
%
% PLIM, is a M x 2 matrix, where the first colume refers to the lower
% boundary for the parameter space and the second colume to the higher one.
%
% IP, is a cell whose length equal to the total number of reactions. The
% j-th cell element should contain the index (or indices) of the
% parameter(s) that contribute to the j-th reaction (only in a
% multiplicative way). For example: ip = {{0}, {1}, {2,3}};
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



[num_react, num_spe]  = size(v);


A = Operator_CFPE_reaction(x_lim, d, v, rate_ini, React, tol_rank);

Ip = tt_eye(2,pd);

for i = 1 : p_var
    
    Z{i} = diag(tt_ones(2,pd)*plim(i,1) + tt_x(2,pd)*ph(i));
    
end

output = [];

for i = 1 : num_react
    
    for j = 1 : p_var % alternating between parameters
        
        if any(cell2mat(ip{i}) == j)
            
            A{i} = tkron(A{i},Z{j});
            
        else
            
            A{i} = tkron(A{i},Ip);
            
        end
        
    end
    
    output = output + A{i};
    
    output = round(output, tol_rank);
    
end


end