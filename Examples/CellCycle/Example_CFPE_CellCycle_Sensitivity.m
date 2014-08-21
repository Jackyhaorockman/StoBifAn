function Example_CFPE_CellCycle_Sensitivity
% EXAMPLE_CFPE_CELLCYCLE_SENSITIVITY Computes the sensitivity indicators
% defined for probability distribution of the cell cycle model.
%
%
% For more information about this example, please see the related
% publication:
% 
% Liao, S., Vejchodsky, T. & Erban, R. (2014). Parameter estimation and 
% bifurcation analysis of stochastic models of gene regulatory networks: 
% tensor-structured methods. arXiv preprint arXiv:1406.7825.
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


d = 8; % No. of grid points in each dimension

d_p = 5; % No. of points for bifurcation parameters

Dim = 5; % No. of dimensions

% the rang of the bifurcation values;
k_bif_ini = 0.25; k_bif_end = 0.4;

h_p = (k_bif_end - k_bif_ini) / (2^d_p-1);  % grid size for parameter.


x = tt_tensor('parametric_steady_state');  % load the tensor-structured probability

S = ones(2^d_p,1); % initiate the sensitivity indicator


% loop through all the grid points for the sensitivity parameter
for para = 1:2^d_p
    
    v_k = k_bif_ini + (para-1)*h_p;
    
    y = tt_extract_1D(x, Dim*d, d_p, para);
    
    y = y / norm(y);
    
    if para == 1;
        
        P_pre = y;
        
    else
        
        S(para) = norm(P_pre - y) / h_p * v_k / norm(y);
        
        P_pre = y;
        
    end
    
end

S(1) = []; % delete the first entry


% plot the solution
figure(1);

    set(gca, 'fontsize', 15);
    
    p_data = linspace(k_bif_ini, k_bif_end, 2^d_p);
    
    p_data(1) = [];
    
    plot(p_data, S, '-ok', 'linewidth', 2);