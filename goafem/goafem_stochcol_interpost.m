function [perrest_u, perrest_z] = ...
    goafem_stochcol_interpost(X_new, errest2s_u, errest2s_z, ...
    gridd_diff, list, rule_id)
% COMMENTS NEED FINISHING
% GOAFEM_STOCHCOL_INTERPOST computes parametric errors


paras_sg_new = stochcol_sg(X_new, rule_id);
L_two_norm_new = stochcol_multilag_Ltwonorm(paras_sg_new, list);
gridd_new = paras_sg_new{4};
[~, IA, IB] = intersect(gridd_diff, gridd_new, 'rows');
perrest_u = errest2s_u(IA)*L_two_norm_new(IB);
perrest_z = errest2s_z(IA)*L_two_norm_new(IB);
