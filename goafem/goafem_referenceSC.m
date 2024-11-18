% Driver for generating reference solution, goal functional related 
% quantities and effectivity indices
%
% SEE ALSO: referenceSC
%           goafem_referenceSC_known
%
% AS, TR; modified 28 June 2024

if nargin(F_rhs{1}) < 3
fprintf('\nRefining parameters...')
startrefTime = tic;

% Reference index set
X_ref = stochcol_getgrid([X; stochcol_margin_set(X)]);
M = size(X_ref, 2);

if Q~=0
    X_ref = [X_ref, ones(size(X_ref, 1), Q)];
    for q=1:Q
        additional_dim = ones(1, M+Q);
        additional_dim(end-q+1) = 2;
        X_ref = [X_ref; additional_dim];
    end
    X_ref = stochcol_getgrid(X_ref);
end

paras_sg_ref = stochcol_sg(X_ref, rule_id);
gridd_ref = paras_sg_ref{4};
clincombiset_ref = paras_sg_ref{5};
indlincombiset_ref = paras_sg_ref{6};
coords_ref = paras_sg_ref{9};

% Reference nodes
figure(11)
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
%paras_fem_ref = stochcol_mesh_refine(MMele_ref, MMedge_ref, ...
%     paras_fem, paras_detail, pmethod);
paras_fem_ref = stochcol_mesh_refine_p2(MMele_ref, MMedge_ref, ...
    paras_fem, paras_detail, pmethod);

fprintf('\nComputing reference solutions...')
% Reference solution calculation
sols_u_ref = zeros(length(paras_fem_ref{1}), size(coords_ref, 1));
sols_z_ref = zeros(length(paras_fem_ref{1}), size(coords_ref, 1));
Gu = zeros(size(coords_ref, 1), 1);
parfor k = 1:size(coords_ref, 1) %parfor
    % FE solution for a given collocation node
    F_rhs_tmp = F_rhs;
    [u_gal, z_gal, A_ref] = goafem_stochcol_fem_solver(coords_ref(k, :), ...
        paras_fem_ref, aa, F_rhs_tmp, G_rhs);
    sols_u_ref(:, k) = u_gal;
    sols_z_ref(:, k) = z_gal;
    Gu(k) = u_gal' * A_ref * z_gal;
end
end
%%
if nargin(F_rhs{1}) < 3
fprintf('\nComputing reference goal functional...')
% Compute reference goal functional
[~,~,~,~,wp_ref] = stochcol_multilag(paras_sg_ref{4}, paras_sg_ref{1}, paras_sg_ref{2}, paras_sg_ref{3}, rule_id, fun_p);
tic;
input_tmp = input;
input_tmp(1) = size(paras_sg_ref{9},2);
KL_DATA = stoch_kl(input_tmp,expansion_type);
if G_rhs{5} >= 0
    G_ref = stochcol_gmatrices(paras_sg_ref{4}, paras_sg_ref{5}, paras_sg_ref{6}, list);
    if G_rhs{5} == 0 
        Mass_ref = stochcol_mass_matrix(paras_fem_ref{1}, paras_fem_ref{2}, G_rhs{1});
        Goal_unc_ref = sum(dot(G_ref, sols_u_ref' * Mass_ref * sols_u_ref));
    elseif G_rhs{5} == 1 
        [~, Mass_ref] = stochcol_mass_matrix(paras_fem_ref{1}, paras_fem_ref{2}, G_rhs{1});
        Goal_unc_ref = sum(dot(G_ref, sols_u_ref' * Mass_ref * sols_u_ref));
    elseif G_rhs{5} == 2
        [~, wz_ref] = goafem_stochcol_fem_setup(paras_fem_ref{1}, paras_fem_ref{2}, a_unit, G_rhs);
        u_col_point_ref = wz_ref.'*sols_u_ref;
        Goal_unc_ref = 100*u_col_point_ref*G_ref*u_col_point_ref.';
    elseif G_rhs{5} == 3
        [~, wz_ref] = goafem_stochcol_fem_setup(paras_fem_ref{1}, paras_fem_ref{2}, a_unit, G_rhs);
        u_col_point_ref = wz_ref.'*sols_u_ref;
        Goal_unc_ref = 100*(u_col_point_ref*G_ref*u_col_point_ref.' - (u_col_point_ref*wp_ref)^2);
    end
    B_ref = goafem_doBilinearForm(KL_DATA, paras_fem_ref, paras_sg_ref, list, ...
        listy, listy2, sols_u_ref, sols_z_ref, F_rhs, G_rhs, rf_type, ...
        X_ref, fun_p, rule_id, aa, polys, 0);

    [~, fz_ref] = goafem_stochcol_fem_setup(paras_fem_ref{1}, paras_fem_ref{2}, a_unit, F_rhs);
    F = fz_ref.'*sols_z_ref*wp_ref;
    Goal_ref = Goal_unc_ref + F - B_ref;   
else
    toc
    tic;
    Goal_ref = 0;
    for k = 1:size(coords_ref,1)
        Goal_ref = Goal_ref + wp_ref(k) * Gu(k);
    end
    Goal_unc_ref = Goal_ref;
    toc
    tic;
    B_ref = goafem_doBilinearForm(KL_DATA, paras_fem_ref, paras_sg_ref, list, ...
            listy, listy2, sols_u_ref, sols_z_ref, F_rhs, G_rhs, rf_type, ...
            X_ref, fun_p, rule_id, aa, polys, 0);
     toc;   
    Goal_ref = 2*Goal_ref - B_ref;
end
else
    Goal_unc_ref = pi*(sqrt(10) -1)/450;
end
%%
%iter
%Goal_unc_ref = 0.124120599;
% Approximate error in goal functional for uncorrected/corrected variants
GoalErr  = abs(Goal_unc_iter - Goal_unc_ref * ones(1,size(Goal_iter,2)));
GoalErr2 = abs(Goal_iter - Goal_unc_ref * ones(1,size(Goal_iter,2)));

% This quantity attempts to give a quantitive description of the 'effect'
% of the bilinear form, but is typically not necessary for plots.
% GoalErr3 = abs(B_iter - Goal_unc_ref * ones(1,size(Goal_iter,2))); 

%if nargin(F_rhs{1}) < 3
%endrefTime = toc(startrefTime);
%fprintf('\n<strong>Elapsed time:</strong> %.2f sec\n',endrefTime);
%end

iterno = default('\nHow many iterations should be plotted? (Default = ALL of them)',size(dof,2));
if iterno > size(dof,2)
    fprintf('\nNo. of iterations too high. Reducing to maximum no.')
    iterno = size(dof,2);
end
if iterno == 0
    return
end
c0 = 0.75*GoalErr2(1)*dof(1)^(3/4);
c4 = 1.00*(dof(iterno))^(2/3)*min([error_iter(1:iterno), error_dir_iter(1:iterno)]);
l7 = min(l5, min(GoalErr(1:iterno)));
l8 = max(l6, max(GoalErr(1:iterno)));
l9 = max(3*dof(1:iterno));
if min(c0/c4,c4/c0) > 0.1
    c0 = (c4 + c0)/2;
    c4 = c0;
end

% Plot of Error in the Goal Functional
dof_iterno = dof(1,1:iterno);
figure(15);
p2 = loglog(...
       [10^-1 10^8], ((c3+c4)/4)*[10^-1 10^8].^(-2/3), 'k--');
hold on
p1 = loglog(dof_iterno, abs(GoalErr(1:iterno)),'+-b', ...
            dof_iterno, abs(GoalErr2(1:iterno)),'s--r', ...
            dof_iterno, error_dir_iter(1:iterno),'d--m');%, ...
           % dof_iterno, error_dir_iter0(1:iterno),'o--g');
grid on
xlabel('degrees of freedom')
ylabel('Goal Functional Error')
legend([p1;p2],{'$|Q(u_{ref}) - Q(u_{iter})|$', ...
                '$|Q(u_{ref}) - \bar{Q}(u_{iter},z_{iter})|$' , ...
                '$(\eta_{iter} + \mu_{iter})(\sigma_{iter} + \tau_{iter})$ (total error)',...
                '$\mathcal{O}(N^{-2/3})$'}, ...
       'Location', 'Best', 'interpreter', 'latex')
axis ([l1 l9 l7 l8])

effs_1 = error_dir_iter(1:iterno)./GoalErr(1:iterno);
effs_2 = error_dir_iter(1:iterno)./GoalErr2(1:iterno);


% Plot of Effectivity indices
figure(16);
loglog(dof_iterno, effs_1, 'o-b', ...
       dof_iterno, effs_2, 'd-r')
hold on
grid on
xlabel('degrees of freedom')
ylabel('Effectivity Indices')
legend('$\bar\eta / |Q(u_{ref}) - \bar{Q}(u_{iter},z_{iter})|$' , ...
       '$\eta / |Q(u_{ref}) - \bar{Q}(u_{iter},z_{iter})|$' , ...
       'Location', 'Best', 'interpreter', 'latex')
axis([min(dof_iterno)/2 max(dof_iterno)*2 min([effs_1,effs_2])/1.1 max([effs_1,effs_2])*1.1])
%% to save .dat files for latex
dof = dof_iterno;
error_d = error_dir_iter(1:iterno);
error_unc = abs(GoalErr(1:iterno));
error_corr = abs(GoalErr2(1:iterno));
straight_line = ((c3+c4)/4)*dof_iterno.^(-2/3);
