function Example_CFPE_FitzhughNagumo_RobustnessAnalysis
% EXAMPLE_CFPE_FITZHUGHNAGUMO_ROBUSTNESSANALYSIS Computes the averaged
% behaviour of the Fitzhugh-Nagumo model with respect to different
% parametric distributions (delta, normal, uniform and bimodal).
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

% !!! NOTE: the simulation parameters here should be consistent with the
% parameters assigned to the function: Example_CFPE_FitzhughNagumo_X2_p4.

d = 8;  % number of grid size in physical space

pd = 7;  % number of grid size in parameter space

x_lim = [-200, 1800; 0, 700];  % computational domain in physical space

h = (x_lim(:,2) - x_lim(:,1)) / (2^d - 1);  % grid size in the state space

pcv = 15/100; % variations around the mean, i.e., x \in [85% mean, 115% mean];

x = tt_tensor('parametric_steady_state'); % loading the parametric steady state distribution

index_para = [65 65 65 65]';  % the grid index of the parameter space to be plotted (for the delta distrbution)

d_tol = [d d pd pd pd pd];  % grid size allocation for state + parameter space.


x1_data = linspace(x_lim(1,1),x_lim(1,2),2^d);

x2_data = linspace(x_lim(2,1),x_lim(2,2),2^d);

% the delta distribution
x_delta = x;

for i = 1 : 4
    
    x_delta = tt_extract_1D(x_delta, sum(d_tol(1: end - i)), d_tol(end-i+1), index_para(5-i));
    
end

x_delta = full(x_delta); x_delta = abs(x_delta);

x_delta = x_delta / sum(x_delta) / h(1) / h(2); % normalisation

x_delta = reshape(x_delta, 2^d, 2^d);

% the normal distribution
x_normal = Robustness_normal(x, d_tol, [1,2], pcv, pcv/3);

x_normal = x_normal / sum(x_normal) / h(1) / h(2);

x_normal = reshape(x_normal, 2^d, 2^d);

% the uniform distribution
x_uniform = Robustness_uniform(x, d_tol, [1,2]);

x_uniform = x_uniform / sum(x_uniform) / h(1) / h(2);

x_uniform = reshape(x_uniform, 2^d, 2^d);

% the bimodal distribution
x_bimodal = Robustness_bimodal(x, d_tol, [1,2], pcv, pcv/3);

x_bimodal = x_bimodal / sum(x_bimodal) / h(1) / h(2);

x_bimodal = reshape(x_bimodal, 2^d, 2^d);


% plotting
figure(1);

    set(gca, 'fontsize', 15);
    
    surf(x2_data, x1_data, x_delta, 'edgecolor', 'none'); 
    
    axis tight; view(2); colorbar;
    
figure(2);

    set(gca, 'fontsize', 15);
    
    surf(x2_data, x1_data, x_normal, 'edgecolor', 'none'); 
    
    axis tight; view(2); colorbar;
    
figure(3);

    set(gca, 'fontsize', 15);
    
    surf(x2_data, x1_data, x_uniform, 'edgecolor', 'none'); 
    
    axis tight; view(2); colorbar;
    
figure(4);

    set(gca, 'fontsize', 15);
    
    surf(x2_data, x1_data, x_bimodal, 'edgecolor', 'none'); 
    
    axis tight; view(2); colorbar;


    
    
% ========== subfunctions =========

    function dis = Robustness_normal(dis, d_tol, d_plt, pcv, CV)
        % ROBUSTNESS_NORMAL Computes the averaged behaviour of the
        % Fitzhugh-Nagumo model interns of normally distributed parametric
        % variations. 
        %
        % The inputs:
        %
        % dis, is the tensor-structured parametric steady state
        % distribution
        %
        % d_tol, is a vector containing the number of grid points in each
        % dimensions
        %
        % d_plot, is a vector (or a number) specifying the index of the 
        % state space to be plotted.
        %
        % pcv, is the number defining the boundaries for the parameters.
        %
        % CV, is the variations.
        
        D = length(d_tol);
        
        n = d_tol;
        
        check = zeros(1,D); check(d_plt) = 1;
        
        while any( check == 0 )
            
            dim = find(check == 0);
            
            [dis, n] = tt_integration_normal(dis, n, dim(1), pcv, CV);
            
            check(dim(1)) = [];
            
        end
        
        dis = full(dis);
        
        dis = abs(dis);
        
    end


    function [output, d] = tt_integration_normal( tensor, d, int_D, pcv, CV )
        % TT_INTEGRATION_NORMAL Integrates the int_D-th parameter using
        % normal distribution.
        
        output = reshape(tensor, [2*ones(1,sum(d(1:int_D-1))), 2^d(int_D),...
            2*ones(1,sum(d(int_D+1:end)))]);
        
        pdf_normal = pdf('Normal', linspace(1-pcv,1+pcv,2^d(int_D)), 1, (CV));
        figure(2); plot(pdf_normal);
        
        output = ttm(output, sum(d(1:int_D-1)) + 1, pdf_normal');
        
        output = reshape(output, [2*ones(1,sum(d(1:int_D-1))), 2*ones(1,sum(d(int_D+1:end)))]);
        
        d(int_D) = [];
        
    end


    function tensor = Robustness_uniform(tensor, dtol, dplot)
        
        % return the matrix for the distribution along ith and jth dimension
        
        D = length(dtol);
        
        n = dtol;
        check = zeros(1,D); check(dplot) = 1;
        
        while any( check == 0 )
            dim = find(check == 0);
            [tensor, n] = tt_integration(tensor, n, dim(1));
            check(dim(1)) = [];
        end
        tensor = full(tensor);
        tensor = abs(tensor);
        
    end


    function tensor = Robustness_bimodal(tensor, dtol, dplot, pcv, CV)
        
        % return the matrix for the distribution along ith and jth dimension
        
        D = length(dtol);
        
        n = dtol;
        check = zeros(1,D); check(dplot) = 1;
        
        while any( check == 0 )
            dim = find(check == 0);
            [tensor, n] = tt_integration_bimodal(tensor, n, dim(1), pcv, CV);
            check(dim(1)) = [];
        end
        tensor = full(tensor);
        tensor = abs(tensor);
        
    end


    function [output, d] = tt_integration_bimodal( tensor, d, int_D, pcv, CV )
        
        % integrate the int_D th dimension.
        
        output = reshape(tensor, [2*ones(1,sum(d(1:int_D-1))), 2^d(int_D),...
            2*ones(1,sum(d(int_D+1:end)))]);
        
        pdf_normal1 = pdf('Normal', linspace(1-pcv,1+pcv,2^d(int_D)), 1-pcv/1, (CV));
        pdf_normal2 = pdf('Normal', linspace(1-pcv,1+pcv,2^d(int_D)), 1+pcv/1, (CV));
        pdf_normal = pdf_normal1 + pdf_normal2;
        pdf_normal = pdf_normal / sum(pdf_normal);
        %figure(2); plot(pdf_normal);
        
        output = ttm(output, sum(d(1:int_D-1)) + 1, pdf_normal');
        
        output = reshape(output, [2*ones(1,sum(d(1:int_D-1))), 2*ones(1,sum(d(int_D+1:end)))]);
        
        d(int_D) = [];
        
    end


end