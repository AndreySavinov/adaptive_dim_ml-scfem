% driver for generating effectivity indices for the multilevel SC-FEM solution 
% requires preloaded reference solution that can be precomputed as follows:
% - run the driver singlelevelSC.m for the given test problem
% - generate reference solution by running the driver referenceSC.m
% - save the output
%   save filename.mat sols_ref paras_sg_ref paras_fem_ref polys G_ref A_unit_coeff
%
% Latest update: AB; 15 December 2022
% Copyright (c) 2022 A. Bespalov, D. Silvester, F. Xu

if ~exist('sols_ref'),
error('Oops .. need to preload reference solution'), end

%------ loop over all the adaptive steps
for k = 1:iter
meshesP = paras_fem_iter{k};
sols_ml = sols_ml_iter{k};
paras_sg = paras_sg_iter{k};
coords=paras_sg{9};
nocp=length(coords(:,1));
% generate union mesh
[xy_union, evt_union] = stochcol_meshes_union(meshesP(1:nocp,1)', nocp);
xy_union1 = xy_union(:,1); xy_union2 = xy_union(:,2);
nvtx=length(xy_union1);
fprintf('\niteration step %g\n',k)
fprintf('%g collocation points\n',nocp)
fprintf('union mesh has %g vertices\n',nvtx)
paras_fem_union = {xy_union,evt_union};
                                              
sols_union=nan(nvtx,nocp);
for kk = 1:nocp
% interpolate solutions to the union mesh
xy_k1 = meshesP{kk,1}(:,1);
xy_k2 = meshesP{kk,1}(:,2);
sols_union(:,kk) = griddata(xy_k1, xy_k2,sols_ml{kk}, xy_union1, xy_union2);
end

% call single grid reference code
    ref_err(k) = stochcol_ref_err(sols_union, paras_sg, ...
        paras_fem_union, sols_ref, paras_sg_ref, paras_fem_ref, ...
        polys, G_ref, A_unit_coeff);
    eff_ind_iter(k) = error_iter(k)/ref_err(k);
    eff_ind_d_iter(k) = error_d_iter(k)/ref_err(k);
fprintf('effectivity is %g\n',eff_ind_d_iter(k))
end  % adaptive step loop
                                              
% Effectivity indices
dof=dof(1,1:iter);
figure(111);
semilogx(dof, eff_ind_iter,'o-k', dof, eff_ind_d_iter,'s-b')
hold on
grid on
xlabel('degrees of freedom')
ylabel('\Theta')
legend('\Theta', '\Theta_{direct}','Location', 'Best')
axis tight
