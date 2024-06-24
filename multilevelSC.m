%MULTILEVELSC solves stochastic diffusion problem using adaptive multilevel SC-FEM
%
% Main driver for running multilevel adaptive stochastic collocation algorithm
%
% Latest update: AB; 13 December 2022
% Copyright (c) 2022 A. Bespalov, D. Silvester, F. Xu

clear all; close all;
stochcol_testproblems;
stochcol_adaptive_global_settingsP1;
if sn ~= 3 
   delta = default('set the error tolerance (default is 6e-3)',6e-3);
else %sn == 3
   delta = default('set the error tolerance (default is 1e-1)',1e-1);
end
adaptmax = default('set the number of adaptive steps (default is 40)',40);
cpointsmax = 7;
tot_err_est = Inf;
tot_err_direct = Inf;
nocpmax = 40; % max number of collocation points (for initialisation purposes)


% preallocation of arrays
    dof = nan(1,adaptmax);
    err_p_iter = nan(1,adaptmax);
    err_s_iter = nan(1,adaptmax);
    error_iter = nan(1,adaptmax);
    
    err_p_d_iter = nan(1,adaptmax);
    err_s_d_iter = nan(1,adaptmax);
    error_d_iter = nan(1,adaptmax);
    
% preallocation of cell data
    sols_ml_iter = cell(1,adaptmax);
    paras_sg_iter = cell(1,adaptmax);
    paras_fem_iter = cell(1,adaptmax);
    
iter = 0; glevel = 0;
startLoopTime = tic;
while tot_err_direct >= delta && iter <= adaptmax
    if iter == 0 % First iteration step
        % Initial index set is for a single collocation point
        X = stochcol_getindexset(0, M);
        % Several attributes of the general sparse grid interpolation given
        % the index set and the sequence of 1D collocation points, such as
        % the underlying grid, the coordinates of the underlying grid, the
        % multivariate Lagrange polynomials (represented by the coefficient
        % list of the single-term multivariate Lagrange polynomials and the
        % set of indexing the degree of the single-term polynomials), etc.
        paras_sg = stochcol_sg(X, rule_id);
        gridd = paras_sg{4};
        clincombiset = paras_sg{5};
        indlincombiset = paras_sg{6};
        coords = paras_sg{9};
        nocp = size(coords,1); % number of collocation points

        % plotting initial collocation point(s)
        figure(901)
        scatter(coords(:,1),coords(:,2),100,'o','filled');
        axis square
        box on
        grid on
        title('initial collocation points in first two directions')
        
        % Setting up the coarsest FE mesh/space
        % Several attributes of the finite element space: vertex coordinate
        % vector, element mapping matrix, boundary vertex vector, boundary
        % elment mapping matrix
        paras_fem0 = stochcol_fem_grid_generator(pmethod, dom_paras);
        % Collection of edge information for flux jump computation and
        % outputs of linear detail space Y grid-generator
        paras_detail0 = stochcol_mesh_detail(paras_fem0);
        
        % preallocation of FE-mesh data arrays
        meshesP = cell(nocpmax,length(paras_fem0));
        meshesPdetail = cell(nocpmax,length(paras_detail0));
        % assigning the coarsest FE mesh to collocation point(s)
        for i=1:nocp
            meshesP(i,:) = paras_fem0;
            meshesPdetail(i,:) = paras_detail0;
        end
        
        % plotting the coarsest FE mesh
        if pmethod == 1 % P1
            plot_mesh(paras_fem0{2}, paras_fem0{1}, 'initial FE mesh');
        elseif pmethod == 2 % P2
            plot_mesh(paras_fem0{6}, paras_fem0{5}, 'initial FE mesh');
        end
        
        % We represent the general sparse grid interpolation of the finite
        % element approximation as a cell of vertex values associated with
        % collocation points in the order they are introduced
        %
        % preallocation of cells for sampled FE solutions
        sols_ml = cell(1,nocp);
        sols_0 = cell(1,nocp);
        
        % For each collocation point, we calculate an estimate of the
        % energy error of the FE solution and a list of
        % marked elements and edges
        %
        % preallocation of these data arrays
        errests = nan(1,nocp);
        elerrs = cell(1,nocp);
        ederrs = cell(1,nocp);
        
        % use parfor to enable parallel computing
        parfor k = 1:nocp 
            % FE solution for each collocation point
            if nargin(rhs) == 3 
                rhs_fun = @(x1, x2) rhs(x1, x2, coords(k, :));
            else
                rhs_fun = rhs;
            end
            [x_gal, ~] = stochcol_fem_solver(coords(k,:), meshesP(k,:), ...
                aa, rhs_fun);
            % energy error of the FE solution and element/edge-based error indicators
            [errest, elerr, ederr, ~, ~] = stochcol_fem_estimator_combmark(pmethod,...
                paras_fem_errest, meshesP(k,:), meshesPdetail(k,:), ...
                x_gal, coords(k, :), aa, aax1, aax2, rhs_fun);
            % updating relevant cells/vectors
            sols_ml{k} = x_gal;
            sols_0{k} = x_gal;
            errests(k) = errest;
            elerrs{k} = elerr;
            ederrs{k} = ederr;
        end
        
        % L2-norm of multivariate Lagrange polynomials
        L_two_norm = stochcol_multilag_Ltwonorm(paras_sg, list);
        % construct the detail index set --- reduced margin
        X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X));
        % Attributes of the general sparse grid interpolation corresponding
        % to the augmented index set (the union of the original index 
        % set and the detail index set)
        paras_sg_full = stochcol_sg([X; X_diff], rule_id);
        % set the sparse grid corresponding to the detail index set
        % the grid is used to label ColPoints
        [gridd_diff, gridd_diff_ind] = setdiff(paras_sg_full{4}, ...
            paras_sg{4}, 'rows');
        coords_diff = paras_sg_full{9}(gridd_diff_ind,:);
        
    else
        fprintf('\n\nIteration %i \n',iter)
        fprintf(['   spatial error indicator', ...
            ' is %10.4e \n'],serrest)
        fprintf(['parametric error indicator', ...
            ' is %10.4e \n'],perrest)
        if serrest < perrest % use e.g. serrest < 4*perrest for more aggressive param. enrichment
%        if serrest_d < perrest_d % use e.g. serrest < 4*perrest for more aggressive param. enrichment
            fprintf('Parametric refinement ... new indices added \n');
            glevel = glevel+1;
            % Refined index set
            disp(Mset)
            X = stochcol_getgrid([X; Mset]);
            % Attributes of the general sparse grid interpolation based on
            % the refined index set
            paras_sg_new = stochcol_sg(X, rule_id);
            % Grid of collocation pointss based on previous index set
            gridd = paras_sg{4}; 
            % Grid and coordinates of collocation points based on the
            % refined index set
            gridd_new = paras_sg_new{4};
            coords_new = paras_sg_new{9};
            
            % FE solutions, energy error estimates, marked elements and
            % edges for collocation points based on the refined index set.
            % All FE solutions and some of the estimates are from previous iteration
            %
            % preallocation of cells for FE mesh and solution data for refined sparse grid
            meshesP_new = cell(size(coords_new, 1),length(paras_fem0));
            meshesPdetail_new = cell(size(coords_new, 1),length(paras_detail0));
            sols_ml_new = cell(1,size(coords_new, 1));
            sols_0_new = cell(1,size(coords_new, 1));
            elerrs_new = cell(1,size(coords_new, 1));
            ederrs_new = cell(1,size(coords_new, 1));
            errests_new = nan(1,size(coords_new, 1));
            
            [~, IA, IB] = intersect(gridd_new, gridd, 'rows');
            % IA - positions of the inherited ColPoints within the new grid
            % IB - positions of the inherited ColPoints within the previously used grid
            [~, IC, ID] = intersect(gridd_new, gridd_diff, 'rows');
            % IC - positions of the new selected 'probed' ColPoints within the new grid
            % ID - positions of the new selected 'probed' ColPoints within the set of all 'probed' ColPoints
            
            % Reuse the FE data from previous iteration
            %
            % the following is simplified from the original code by FX
            % to avoid the loop over k = 1:length(IA) (any pitfalls?)
            meshesP_new(IA,:) = meshesP(IB,:);
            meshesPdetail_new(IA,:) = meshesPdetail(IB,:);
            sols_ml_new(IA) = sols_ml(IB);
            sols_0_new(IA) = sols_0(IB);
            elerrs_new(IA) = elerrs(IB);
            ederrs_new(IA) = ederrs(IB);
            errests_new(IA) = errests(IB);
            
            % Assign the coarsest mesh to newly added collocation points
            for i=1:length(IC)
                meshesP_new(IC(i),:) = paras_fem0; % meshesP_diff(ID,:);
                meshesPdetail_new(IC(i),:) = paras_detail0; % meshesPdetail_diff(ID,:);
            end
            % reuse the FE solutions on the coarsest mesh for newly added collocation points
            sols_0_new(IC) = sols_0_diff(ID);
            % initialise the corresponding error estimates by infinity
            errests_new(IC) = deal(Inf);

            % Updating the sparse grid
            paras_sg = paras_sg_new;
            gridd = paras_sg{4};
            clincombiset = paras_sg{5};
            indlincombiset = paras_sg{6};
            coords = paras_sg{9};          
            % updating the number of collocation ponts
            nocp = size(coords,1); % number of collocation points
            
            % Updating all FE data
            meshesP(1:nocp,:) = meshesP_new;
            meshesPdetail(1:nocp,:) = meshesPdetail_new;
            sols_ml = sols_ml_new;
            sols_0 = sols_0_new;
            elerrs = elerrs_new;
            ederrs = ederrs_new;
            errests = errests_new;
            
            % 2-norm of multivariate Lagrange polynomials
            L_two_norm = stochcol_multilag_Ltwonorm(paras_sg, list);
            % construct the detail index set --- reduced margin
            X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X));
%            X_diff = stochcol_margin_set(X); % uncomment for the (full) margin
            % Attributes of the general sparse grid interpolation corresponding
            % to the new index set, which is the union of the original index
            % set and the detail index set
            paras_sg_full = stochcol_sg([X; X_diff], rule_id);
            
            % Computing FE approximations for newly added collocation points
            % (via SOLVE -> ESTIMATE -> MARK -> REFINE loop)
            %
            % setting a tolerance based on the error estimates for other
            % (i.e., older) collocation points
            tol_diff = sum(errests(IA)'.*L_two_norm(IA))/length(IA); % alternatively, use min(errests(IA)'.*L_two_norm(IA))
            %
            for kk = 1:length(IC)
                k = IC(kk);
                mesh_k = meshesP(k,:);
                mesh_k_detail = meshesPdetail(k,:);
                %refin_count=0; % uncomment for counting the number of refinements
                if nargin(rhs) == 3 
                    rhs_fun = @(x1, x2) rhs(x1, x2, coords(k, :));
                else
                    rhs_fun = rhs;
                end
                while errests(k)*L_two_norm(k) >= tol_diff
                    [x_gal, ~] = stochcol_fem_solver(coords(k, :), mesh_k, aa, rhs_fun);
                    sols_ml{k} = x_gal;
                    meshesP(k,:) = mesh_k;
                    meshesPdetail(k,:) = mesh_k_detail;
                    [errests(k), elerrs{k}, ederrs{k}, MMele, MMedge] = ...
                        stochcol_fem_estimator_combmark(pmethod, paras_fem_errest, mesh_k, mesh_k_detail, ...
                        x_gal, coords(k,:), aa, aax1, aax2, rhs_fun);
                    mesh_k = stochcol_mesh_refine(MMele, MMedge, mesh_k, ...
                        mesh_k_detail, pmethod);
                    mesh_k_detail = stochcol_mesh_detail(mesh_k);
                    %refin_count = refin_count +1; % uncomment for counting the number of refinements
                end
                %refin_count-1 % uncomment for printing out the number of refinements
            end
            
            % set the sparse grid corresponding to the detail index set
            % the grid is used to label 'probed' ColPoints
            [gridd_diff, gridd_diff_ind] = setdiff(paras_sg_full{4}, paras_sg{4}, 'rows');
            coords_diff = paras_sg_full{9}(gridd_diff_ind,:);

        else
            fprintf('Spatial refinement...');
            % Mesh refinement
            for k = 1:nocp
                if ~isempty(MMeles{k})
                    %fprintf('\nRefinement of mesh %i...',k);
                    meshesP(k,:) = stochcol_mesh_refine(MMeles{k}, MMedges{k}, meshesP(k,:), ...
                        meshesPdetail(k,:), pmethod);
                    % Collecting edge information for flux jump computation and
                    % outputs of linear detail space Y grid-generator based on
                    % refined mesh(es)
                    meshesPdetail(k,:) = stochcol_mesh_detail(meshesP(k,:));
                    %fprintf(' done');
                end
            end
            fprintf(' done\n\n');
            
            % Preallocating FE data arrays
            sols_ml = cell(1,nocp);
            errests = nan(1,nocp);
            elerrs = cell(1,nocp);
            ederrs = cell(1,nocp);
            % Computing new FE solution, error estimates and indicators
            parfor k = 1:nocp
                if nargin(rhs) == 3 
                    rhs_fun = @(x1, x2) rhs(x1, x2, coords(k, :));
                else
                    rhs_fun = rhs;
                end
                [x_gal, ~] = stochcol_fem_solver(coords(k,:), meshesP(k,:), aa, rhs_fun);
                % energy error of the FE solution and element/edge-based error indicators
                [errest, elerr, ederr, ~, ~] = stochcol_fem_estimator_combmark(pmethod,...
                    paras_fem_errest, meshesP(k,:), meshesPdetail(k,:), ...
                    x_gal, coords(k, :), aa, aax1, aax2, rhs_fun);
                sols_ml{k} = x_gal;
                errests(k) = errest;
                elerrs{k} = elerr;
                ederrs{k} = ederr;
            end
        end
    end %------------- of adaptive loop

    output_matrix1 = [errests' L_two_norm errests'.*L_two_norm];

    % assembling local error indicators for spatial marking
    serrests = cell(1,nocp);
    nEstP = nan(1,nocp);
    if markedgelem == 1 % marking elements
        parfor k = 1:nocp
            nEstP(1,k) = length(elerrs{k});
            serrests{k} = elerrs{k} * L_two_norm(k,1);
        end
    else % marking edges
        parfor k = 1:nocp
            nEstP(1,k) = length(ederrs{k});
            serrests{k} = ederrs{k} * L_two_norm(k,1);
        end
    end
    
    % spatial marking
    MsetP = cell(1,nocp);
    MMeles = MsetP;
    MMedges = MsetP;
    cumsumnEstP = cumsum(nEstP);
    %
    est_vec = vertcat(serrests{:}); % vector with all spatial indicators
    subvec_glob2loc = [0, cumsumnEstP];
    subvec_start = subvec_glob2loc+1;
    subvec_end = cumsumnEstP;
%    marked = marking_strategy_fa(est_vec,markstrat,pmthreshold);
    marked = dorfler_marking((1:cumsumnEstP(end))', est_vec, smthreshold);
    for k = 1:nocp
        % set of marked spatial indicators
        MsetP{k} = marked((marked>=subvec_start(k))&(marked<=subvec_end(k)))-subvec_glob2loc(k);
        if size(MsetP{k},1)==0
            MsetP{k} = [];
        end
        % overall set of marked elements and edges
        [MMeles{k},MMedges{k}] = get_all_marked_elem(meshesPdetail{k,7},meshesPdetail{k,4},MsetP{k},markedgelem);
    end
    % end of spatial marking

    % spatial error indicator (assembled from components)
    serrest = errests*L_two_norm;

    % parametric error indicator
    %
    % preallocating data arrays for coarsest-mesh FE computation
    sols_0_diff = cell(1, size(coords_diff,1));
    errest2s = zeros(1, size(coords_diff,1));
    % setting the unit coefficient for computing spatial H1 seminorm
    a_unit = @(x1, x2) ones(size(x1));
    [A_unit_0,~] = stochcol_fem_setup(paras_fem0{1}, paras_fem0{2}, a_unit, a_unit);
    % the matrix of the coarsest-mesh FE solution vectors for current ColPoints
    sols_0_matrix = cell2mat(sols_0);
    
    % compute the coarsest-mesh FE approximations for prospective ColPoints
    for k = 1:size(coords_diff,1)
        if nargin(rhs) == 3 
            rhs_fun = @(x1, x2) rhs(x1, x2, coords_diff(k, :));
        else
            rhs_fun = rhs;
        end
        [x_gal_0, ~] = stochcol_fem_solver(coords_diff(k, :), paras_fem0, aa, rhs_fun);
        sols_0_diff{k} = x_gal_0;
        % vector of multivariate Lagrange polynomials for the collocation point
        LL_diff = stochcol_getinterpolant_2(gridd, ...
            clincombiset, indlincombiset, coords_diff(k,:), polys);
        % evaluating the current-sparse-grid-interpolated coarsest-mesh FE solutions
        % at the prospective ColPoint
        uyy = sols_0_matrix * LL_diff;
        % H1 seminorm of the interpolation error at the prospective ColPoint
        errest2s(k) = sqrt((x_gal_0 - uyy)' * A_unit_0 * (x_gal_0 - uyy));
    end
    
    % indexwise parametric error estimates based on the detail collocation points
    perrests = stochcol_est_parametric(X_diff, errest2s, gridd_diff, list, rule_id);
    % parametric error indicator
    perrest = sum(perrests);
    % marked index set from the detail index set
    [Mset, ~] = dorfler_marking(X_diff, perrests, pmthreshold);
 
    % total error estimate
    tot_err_est = serrest + perrest;
    
    % Compute direct estimates of error
    %
    % spatial
    serrest_d = stochcol_direct_estimator_ml(paras_sg, ...
                meshesP, meshesPdetail, sols_ml, list, pmethod, rhs, aa);
    % parametric
    paras_sg_diff = stochcol_sg(X_diff, rule_id);
    G_diff = stochcol_gmatrices(paras_sg_diff{4}, ...
            paras_sg_diff{5},paras_sg_diff{6}, list);
    ncpts = length(G_diff(:,1));
    ncptsf = nocp + size(coords_diff, 1);
    gridd_old_ind = setdiff([1:ncptsf]',gridd_diff_ind);
    % construct the array containing coarsest-mesh FE solution vectors
    % for all ColPoints (current and prospective)
    sols_0_all = nan(size(paras_fem0{1},1),ncptsf);
    sols_0_all(:,gridd_old_ind) = sols_0_matrix;
    sols_0_all(:,gridd_diff_ind) = cell2mat(sols_0_diff);
    if ncptsf > ncpts
        fprintf('redundant sparse grid solutions ..\n')
        [gdiff,indx]=setdiff(paras_sg_full{4}, paras_sg_diff{4}, 'rows');
        %gdiff
        iactive=setdiff([1:ncptsf]',indx);
        sols_0_all=sols_0_all(:,iactive);
        elseif ncptsf < ncpts, error('Oops ..  fatal logic issue'),
    end
    % compute the Bochner norm via hierarchical surplus
    perrest_d = sqrt(sum(dot(G_diff, sols_0_all' * A_unit_0 * sols_0_all)));

    % total direct error estimate
    tot_err_direct = serrest_d + perrest_d;

    % output error estimates
    if iter==0, fprintf('\n\nIteration %i \n',iter), end
    fprintf(['   spatial error estimate', ...
             ' is %10.4e  vs  %10.4e (spatial indicator)\n'],serrest_d,serrest)
    fprintf(['parametric error estimate', ...
             ' is %10.4e  vs  %10.4e (parametric indicator)\n'],perrest_d,perrest)
    fprintf(['overall estimate from indicators', ...
             ' is %10.4e'],tot_err_est)
    fprintf(['\n   overall direct error estimate', ...
    ' is <strong>%10.4e</strong>\n'],tot_err_direct)
    
    iter = iter + 1;
    
    % saving data corresponding to this iteration
    err_p_iter(iter) = perrest;
    err_s_iter(iter) = serrest;
    error_iter(iter) = tot_err_est;
    %
    err_p_d_iter(iter) = perrest_d;
    err_s_d_iter(iter) = serrest_d;
    error_d_iter(iter) = tot_err_direct;
    %
    paras_sg_iter{iter} = paras_sg;
    %
    paras_fem_iter{iter} = meshesP;
    sols_ml_iter{iter} = sols_ml;
    % computing and saving the number of dof
    dof(iter) = 0;
    for k=1:nocp
        dof(iter) = dof(iter) + length(sols_ml{k});
    end
end

endLoopTime = toc(startLoopTime);

% plot error estimates
figure(98);
loglog(dof, error_iter, 'o-k', dof, error_d_iter, 's-b')
hold on, grid on
xlabel('degrees of freedom')
ylabel('estimated error')
legend('$\bar\eta$ (total)', '$\eta$ (total)', ...
       'Location', 'Best','interpreter','latex');
axis tight
figure(99);
loglog(dof, error_d_iter, 'o-k', dof, err_s_d_iter, 's-b', dof, err_p_d_iter, 'd-r')
hold on
grid on
xlabel('degrees of freedom')
ylabel('estimated error')
legend('$\eta$ (total)', '$\mu$ (spatial)', '$\tau$ (parametric)', ...
       'Location', 'Best','interpreter','latex')
axis tight
%
figure(100);
loglog(dof, error_iter, 'o-k', dof, err_s_iter, 's-b', dof, err_p_iter, 'd-r')
hold on
grid on
xlabel('degrees of freedom')
ylabel('estimated error')
legend('$\bar\eta$ (total)', '$\bar\mu$ (spatial)', '$\bar\tau$ (parametric)', ...
       'Location', 'Best','interpreter','latex');
axis tight

% plotting the final mesh for the first collocation point
paras_fem = meshesP(1,:);
if pmethod == 1 % P1
    plot_mesh(paras_fem{2}, paras_fem{1}, 'Final FE mesh1');
    pause(1)
elseif pmethod == 2
    plot_mesh(paras_fem{6}, paras_fem{5}, 'Final FE mesh1');
end

% final collocation points
coords = paras_sg{9};
figure(988)
scatter(coords(:,1),coords(:,2),100,'o','filled');
axis square
box on
grid on
title('final collocation points in first two directions')

%% postprocess to generate statistics
[xy_union, evt_union] = stochcol_meshes_union(meshesP(1:nocp,1)', nocp);
xy_union1 = xy_union(:,1);
xy_union2 = xy_union(:,2);
%
G = stochcol_gmatrices(paras_sg{4}, paras_sg{5}, paras_sg{6}, list);

[~,~,~,~,multiweights_pdf]...
        = stochcol_multilag(paras_sg{4}, paras_sg{1}, paras_sg{2}, paras_sg{3}, rule_id, fun_p);

MEAN = zeros(size(xy_union,1), 1);
Moment_U_squared  = zeros(size(xy_union,1), 1);

sols_ml_union = nan(size(xy_union,1), nocp);
parfor k = 1:nocp
    % total pdf
    xy_k1 = meshesP{k,1}(:,1);
    xy_k2 = meshesP{k,1}(:,2);
    sols_ml_union(:, k)= griddata(xy_k1, xy_k2, sols_ml{k}, xy_union1, xy_union2);
    MEAN = MEAN + multiweights_pdf(k) * sols_ml_union(:, k); 
end

for k=1:nocp
    for l=1:nocp
        Moment_U_squared = Moment_U_squared + sols_ml_union(:,k).*sols_ml_union(:,l)*G(k,l); % second order moment
    end
end
VAR = Moment_U_squared - MEAN.^2;

STD = VAR.^(1/2);

fprintf('\nFinal sparse grid\n')
disp([paras_sg{4}])

if tot_err_direct <= delta
    fprintf('Tolerance was reached in %g iterations',iter-1)
    fprintf('\n    after %g parametric refinements',glevel)
else
    fprintf('<strong>Tolerance was not reached after %g iterations!</strong>',iter-1)
    fprintf('\nNumber of parametric refinements: %g',glevel)
end

fprintf('\n              Mean maximum %9.6f\n',max(MEAN))
fprintf('          Variance maximum %9.6f\n',max(VAR))
fprintf('Standard Deviation maximum %9.6f\n',max(STD))

if nargin(rhs) == 3
    % QoI computation
    QQ = stochcol_mass_matrix(xy_union, evt_union);
    QoI = sum(dot(G, sols_ml_union' * QQ * sols_ml_union));
    [A_unit_union,~] = stochcol_fem_setup(xy_union, evt_union, a_unit, a_unit);
    energy = sqrt(sum(dot(G, sols_ml_union' * A_unit_union * sols_ml_union)));
    fprintf('                      QoI is %12.9f\n',QoI/16)
    fprintf('          discrete energy is %12.9f\n',energy)
end
fprintf('\n<strong>Total elapsed time:</strong> %.2f sec\n\n',endLoopTime);
%fprintf('To compute a reference solution run the script <strong>referenceSCml</strong>\n\n');

% plot final solution
if pmethod == 1% P1
    sc_plotsol_p1_ml(dom_type, MEAN, VAR, evt_union, xy_union);
end
