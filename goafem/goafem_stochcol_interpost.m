function [perrest_u, perrest_z] = ...
    goafem_stochcol_interpost(X_new, errest2s_u, errest2s_z, ...
    gridd_diff, list, rule_id)
% GOAFEM_STOCHCOL_INTERPOST computes parametric errors.
%
% input:
%            X_new    reduced margin
%       errest2s_u    spatial primal error corresponding to each new collocation point
%       errest2s_z    spatial dual error corresponding to each new collocation point
%       gridd_diff    indices of each collocation point (sparse grid parameter)
%             list    list of integrals of Lagrange polynomials
%          rule_id    indicator of either Leja or CC points
%
% output:
%        perrest_u   parametric error estimation of primal solution
%        perrest_z   parametric error estimation of dual solution
%
% SEE ALSO: goafem_singlelevelSC
%
% TR; 28 September 2022


paras_sg_new = stochcol_sg(X_new, rule_id);
L_two_norm_new = stochcol_multilag_Ltwonorm(paras_sg_new, list);
gridd_new = paras_sg_new{4};
[~, IA, IB] = intersect(gridd_diff, gridd_new, 'rows');
perrest_u = errest2s_u(IA)*L_two_norm_new(IB);
perrest_z = errest2s_z(IA)*L_two_norm_new(IB);
