function output = Initial_CFPE_para(x_lim, d, task, mean, sigma, p_var, pd)
% INITIAL_CFPE_PARA  Construct the initial guess (or condition) for the
% tensor-structured parametric Fokker-Planck operator (see function
% Operator_CFPE_para). 
%
% Inputs:
%
% P_VAR, is the total number of (varying) parameters, i.e., the
% dimensionality of the parameter space.
%
% PD, determines the number of discrete sample points for each parameter is
% 2^PD.
%
% X_LIM, D, TASK, MEAN, SIGMA, (see function: Initial_CFPE).
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

A = Initial_CFPE(x_lim, d, task, mean, sigma);

I = tt_eye(2, p_var * pd);

output = tkron(A, I);

end