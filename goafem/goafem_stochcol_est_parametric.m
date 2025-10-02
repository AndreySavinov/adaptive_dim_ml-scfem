function [perrests_u, perrests_z] = ...
    goafem_stochcol_est_parametric(X_diff, errest2s_u, errest2s_z, ...
    gridd_diff, list, rule_id, varargin)
% GOAFEM_STOCHCOL_EST_PARAMETRIC computes indexwise parametric error estimator
%
% input:
%           X_diff    reduced margin
%       errest2s_u    spatial primal error corresponding to each new collocation point
%       errest2s_z    spatial dual error corresponding to each new collocation point
%       gridd_diff    indices of each collocation point (sparse grid parameter)
%             list    list of integrals of Lagrange polynomials
%          rule_id    indicator of either Leja or CC points
%
% output:
%       perrests_u    parametric error estimation of primal solution per
%                     index of X_dfiff 
%       perrests_z    parametric error estimation of dual solution
%                     index of X_dfiff 
%
% SEE ALSO: goafem_singlelevelSC
%
% TR; 28 September 2022

N = size(X_diff,1);
perrests_u = zeros(1, N);
perrests_z = zeros(1, N);
for k = 1:size(X_diff,1)
    X_new = X_diff(k, :);
    [perrests_u(k), perrests_z(k)] = ...
        goafem_stochcol_interpost(X_new, errest2s_u, errest2s_z, ...
        gridd_diff, list, rule_id);
end