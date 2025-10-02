% driver for generating reference SC solution and effectivity indices
%
% Latest update: AB; 22 November 2022
% Copyright (c) 2022 A. Bespalov, D. Silvester, F. Xu

fprintf('\ncomputing reference solution... ')

% Reference index set
w_ref = max(sum(X,2)) - M; % + 1;
X_ref = stochcol_getindexset(w_ref, M);
paras_sg_ref = stochcol_sg(X_ref, rule_id);
gridd_ref = paras_sg_ref{4};
clincombiset_ref = paras_sg_ref{5};
indlincombiset_ref = paras_sg_ref{6};
coords_ref = paras_sg_ref{9};

% Reference nodes
figure(999)
scatter(coords_ref(:,1),coords_ref(:,2),100,'o','filled');
axis square
title('reference nodes')

G_ref = stochcol_gmatrices(gridd_ref, clincombiset_ref, ...
    indlincombiset_ref, list);
% Reference solution is calulated on the uniformly refined mesh based on 
% the final mesh and using P2 finite elements
% Refine all elements of the final mesh
MMele_ref = (1:size(paras_fem{2}, 1))';
MMedge_ref = (1:size(paras_detail{7}, 1))';
paras_fem_ref = stochcol_mesh_refine_p2(MMele_ref, MMedge_ref, ...
    paras_fem, paras_detail, pmethod);
% Reference solutions
sols_ref = zeros(length(paras_fem_ref{1}), size(coords_ref, 1));
parfor k = 1:size(coords_ref, 1)
    % FE solution for a given collocation node
    [x_gal, ~] = stochcol_fem_solver(coords_ref(k, :), paras_fem_ref, ...
        aa, F_rhs{1});
    sols_ref(:, k) = x_gal;
end
fprintf('done.\n')
%%
fprintf('\ncomputing and plotting effectivity indices... ')
% Reference error is measured in Soblev norm, so we need the stiffness
% matrix corresponding to the unity diffusion coefficient
unit_coeff = @(x1, x2) ones(size(x1));
[A_unit_coeff,~] = stochcol_fem_setup(paras_fem_ref{1}, paras_fem_ref{2}, ...
    unit_coeff, F_rhs{1});
parfor k = 1:iter
    ref_err(k) = stochcol_ref_err(sols_u_iter{k}, paras_sg_iter{k}, ...
        paras_fem_iter{k}, sols_u_ref, paras_sg_ref, paras_fem_ref, ...
        polys, G_ref, A_unit_coeff);
    eff_ind_iter(k) = error_iter_u(k)/ref_err(k);
    eff_ind_d_iter(k) = error_dir_iter_u(k)/ref_err(k);
end
%%
% Effectivity indices
dof=dof(1,1:iter);
figure(97);
semilogx(dof, eff_ind_iter,'o-k', dof, eff_ind_d_iter,'s-b')
hold on
grid on
xlabel('degree of freedom')
ylabel('\Theta')
legend('\Theta', '\Theta_{direct}', 'Location', 'Best')
axis tight

fprintf('done.\n')
