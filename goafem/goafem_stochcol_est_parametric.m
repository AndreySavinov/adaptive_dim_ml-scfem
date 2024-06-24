function [perrests_u, perrests_z] = ...
    goafem_stochcol_est_parametric(X_diff, errest2s_u, errest2s_z, ...
    gridd_diff, list, rule_id)
% GOAFEM_STOCHCOL_EST_PARAMETRIC indexwise parametric error estimator
%
% COMMENTS NEED FINISHING

N = size(X_diff,1);
perrests_u = zeros(1, N);
perrests_z = zeros(1, N);
for k = 1:size(X_diff,1)
    X_new = X_diff(k, :);
    [perrests_u(k), perrests_z(k)] = ...
        goafem_stochcol_interpost(X_new, errest2s_u, errest2s_z, ...
        gridd_diff, list, rule_id);
end