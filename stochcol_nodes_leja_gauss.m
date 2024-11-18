function [nodes] = stochcol_nodes_leja_gauss(levels, varargin)
%STOCHCOL_NODES_LEJA_GAUSS computes the Leja nodes in [-\infty,+\infty]
%
% input:
%     levels    vector of levels
%     varargin  indicator whether the points should be computed from
%               scratch or use precomputed points
%
% output:
%     nodes     cell of Leja collocation points
%
% The function returns the Leja interpolation sparse grid nodes for a given
% set of levels 'i' (i>0).
%
% At level i, the number of collocation points are:
%
%   mi := 2*i-1
%
% Given the first Leja points X(1) = 0, X(2) = 1, X(3) = -1 the Leja points are defined
% recursively in such a way that the function
%
% F(Y,X(1:N-1)) = ABS(EXP(PROD(Y-X(1:N-1)))
%
% is maximized over [0,+\infty] by Y = X(N). The produced Leja points are symmetric
%
%   SCFEM function: AS 28 June 2024
% Copyright (c) 2024 A. Savinov

%%%%%%% 
% TO LOAD NODES

if numel(varargin) == 0
    try
        load('precomputation_leja_9_gaussian_full10_1.mat', 'nodes')
    catch
        load('precomputation_leja9_gaussian_full_sig_1.mat', 'nodes')
    end

    precomputed_nodes = nodes{1};
    clear nodes;
    nodes = cell(length(levels),1);
    for i = 1:length(levels)
        nodes{i} = precomputed_nodes(1:(2*levels(i)-1))';
    end
else
    %%%%%%%%
    % TO COMPUTE NODES
    nodes = cell(length(levels),1);
    n = max(levels);
    xleja = zeros(2*n-1,1);
    % generate leja points sequence
    xleja(1) = 0;
    xleja(2) = 1;
    xleja(3) = -1;
    x_potential = 0:1e-6:10;
    optimized_func = @(x) exp(-x.^2/2).*abs(x).*abs(x - xleja(2)).*abs(x - xleja(3));
    format long;
    for k = 4:2:2*n-1
        [~, max_ind] = max(optimized_func(x_potential));

        x_tiny = (x_potential(max_ind)-(1e-6)):1e-14:(x_potential(max_ind)+(1e-6));
        [~, max_ind] = max(optimized_func(x_tiny));

        xleja(k) = x_tiny(max_ind);
        xleja(k+1) = -x_tiny(max_ind);

        optimized_func = @(x) optimized_func(x).*abs(x - xleja(k)).*abs(x - xleja(k+1));
    end
    xleja = xleja(1:2*n-1);
    % choose set of leja points by level
    for i = 1:length(levels)
        nodes{i} = xleja(1:(2*levels(i)-1))';
    end
end