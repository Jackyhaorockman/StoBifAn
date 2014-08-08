function CFPE = Operator_CFPE(x_lim, d, v, rate, React, tol_rank, task)
% OPERATOR_CFPE  Construct the finte difference operator of the chemical
% Fokker-Planck equation (in either sparse matrix or tensor train format)
% for mass-action reactions.
%
% CFPE = CFPE_operator(X_LIM, D, V, RATE, REACT)
%
% The Direchlet boundary condition is applied for dimensions involved.
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
% TAKS, choose 1, if you want tensor operator, and 2 for sparse matrix
% format. WARNING: do not choose 2 if you are solving for problems with
% more than 3 dimensions.
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


switch task
    
    case 1 % tensor format
        
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
        
        
        % compute the terms including the first derivatives and the second
        % derivatives
        D = cell(num_spe,1);
        DD = cell(num_spe,1);
        
        for i = 1 : num_spe
            
            for j = 1 : num_react
                
                if v(j,i) ~= 0
                    
                    D{i} = D{i} + v(j,i)*alpha{j};
                    D{i} = round(D{i}, tol_rank);
                    
                    DD{i} = DD{i} + v(j,i)^2*alpha{j};
                    DD{i} = round(DD{i}, tol_rank);
                    
                end
                
            end
            
        end
        
        
        % computing the terms including the mixed derivatives of second order.
        DDD = cell(num_spe,num_spe);
        
        for i1 = 1 : num_spe-1
            
            for i2 = i1+1 : num_spe
                
                for j = 1 : num_react
                    
                    if v(j,i1)*v(j,i2) ~= 0
                        
                        DDD{i1,i2} = DDD{i1,i2} + v(j,i1)*v(j,i2)*alpha{j};
                        DDD{i1,i2} = round(DDD{i1,i2}, tol_rank);
                        
                    end
                    
                end
                
            end
            
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
        
        % Construct the Fokker-Planck operator by summing up the diffusion and the
        % drift terms.
        Diffusion = [];
        Drift = [];
        
        for i = 1 : num_spe
            
            if ~isempty( DD{i} )
                
                Diffusion = Diffusion + 0.5*Dxx{i}*DD{i};
                Diffusion = round(Diffusion, tol_rank);
                
            end
            
            if ~isempty( D{i} )
                
                Drift = Drift - Dx{i}*D{i};
                Drift = round(Drift, tol_rank);
                
            end
            
            for i2 = i+1 : num_spe
                
                if ~isempty( DDD{i,i2} )
                    
                    Drift = Drift + Dxy{i,i2}*DDD{i,i2};
                    Drift = round(Drift, tol_rank);
                    
                end
                
            end
            
        end
        
        % here it is, with lower tensor rank.
        CFPE = round(Diffusion + Drift, tol_rank);
        
    case 2 % sparse matrix format
        
        n = 2^d;
        h = (x_lim(:,2) - x_lim(:,1)) / (n-1);
        
        [num_react, num_spe]  = size(v);
        
        Lap = spdiags([ ones(n, 1), -2*ones(n, 1), ones(n, 1)], [-1 0 1], n, n);
        Cdiff = spdiags([ -ones(n,1), ones(n,1) ], [-1, 1], n, n);
        I = speye(n);
        
        X = cell(num_spe,1);
        
        for i = 1 : num_spe
            
            X{i} = spdiags( x_lim(i,1) + h(i)*[0:(n-1)]', 0, n, n );
            
        end
        
        
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
                        
                    end
                    
                end
                
            end
            
            temp = 1;
            
            for j = 1 : num_spe
                
                temp = kron(temp, alpha_i{j});
                
            end
            
            alpha{i} = temp;
            
        end
        
        
        for i = 1 : num_react
            
            alpha{i} = rate(i) * alpha{i};
            
        end
        
        
        D = cell(num_spe,1);
        DD = cell(num_spe,1);
        
        for i = 1 : num_spe
            
            for j = 1 : num_react
                
                if v(j,i) ~= 0
                    
                    if isempty( D{i} )
                        
                        D{i} = sparse(n^2,n^2);
                        DD{i} = sparse(n^2,n^2);
                        
                    end
                    
                    D{i} = D{i} + v(j,i)*alpha{j};
                    DD{i} = DD{i} + v(j,i)^2*alpha{j};
                    
                end
                
            end
            
        end
        
        
        DDD = cell(num_spe,num_spe);
        
        for i1 = 1 : num_spe-1
            
            for i2 = i1+1 : num_spe
                
                for j = 1 : num_react
                    
                    if v(j,i1)*v(j,i2) ~= 0
                        
                        if isempty( DDD{i1,i2} )
                            
                            DDD{i1,i2} = sparse(n^2,n^2);
                            
                        end
                        
                        DDD{i1,i2} = DDD{i1,i2} + v(j,i1)*v(j,i2)*alpha{j};
                        
                    end
                    
                end
                
            end
            
        end
        
        
        Dx = cell(num_spe,1);
        Dxx = cell(num_spe,1);
        Dxy = cell(num_spe,num_spe+1);
        
        for i = 1 : num_spe
            
            Dx{i} = (1/(h(i)));
            Dxx{i} = (1/h(i)^2);
            
            for j = 1 : num_spe
                
                if j == i
                    
                    Dx{i} = kron(Dx{i}, Cdiff);
                    Dxx{i} = kron(Dxx{i}, Lap);
                    
                else
                    
                    Dx{i} = kron(Dx{i}, I);
                    Dxx{i} = kron(Dxx{i}, I);
                    
                end
                
            end
            
            
            if i < num_spe
                
                for j = i+1:num_spe
                    
                    Dxy{i,j} = (1/(h(i))) * (1/(h(j)));
                    
                    for j1 = 1 : num_spe
                        
                        if any(j1 == [i,j])
                            
                            Dxy{i,j} = kron(Dxy{i,j}, Cdiff);
                            
                        else
                            
                            Dxy{i,j} = kron(Dxy{i,j}, I);
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        
        % CFPE operator
        Diffusion = sparse(n^2,n^2);
        Drift = sparse(n^2,n^2);
        
        for i = 1 : num_spe
            
            if ~isempty( DD{i} )
                
                Diffusion = Diffusion + 0.5*Dxx{i}*DD{i};
                
            end
            
            if ~isempty( D{i} )
                
                Drift = Drift - Dx{i}*D{i};

            end
            
        end
        CFPE = Diffusion + Drift;
        
end

end