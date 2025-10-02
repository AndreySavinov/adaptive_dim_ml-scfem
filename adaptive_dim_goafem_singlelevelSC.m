%DIM_ADAPT_GOAFEM_SINGLELEVELSC solves goal-oriented stochastic diffusion problem using adaptive single-level SC-FEM
%
% Main driver for running goal-oriented single-level adaptive stochastic collocation algorithm
%
% Latest update: AS; 12 June 2025
% Copyright (c) 2024 A. Bespalov, T. Round, A. Savinov

addpath('goafem')
adaptive_dim_goafem_stochcol_adaptive_global_settings;
if sn==1
    delta = default('\nSet the error tolerance (default is 3e-4)',3e-4);
elseif sn==2
    delta = default('\nSet the error tolerance (default is 4e-4)',4e-4);
elseif sn==3
    delta = default('\nSet the error tolerance (default is 2e-4)',2e-4);
elseif sn==4
    delta = default('\nSet the error tolerance (default is 8e-5)',8e-5);
else
    delta = default('\nSet the error tolerance (default is 1e-4)',1e-4);
end
adaptmax = default('\nSet the number of adaptive steps (default is 50)', 50);
cpointsmax = 7;
tot_err_est_u = Inf;
tot_err_direct_u = Inf;
tot_err_est_z = Inf;
tot_err_direct_z = Inf;
tot_err_direct = Inf;
%tot_err_direct0 = Inf;

% preallocation of arrays
    dof = nan(1,adaptmax);
    err_p_iter_u = nan(1,adaptmax);
    err_s_iter_u = nan(1,adaptmax);
    error_iter_u = nan(1,adaptmax);
    err_p_iter_z = nan(1,adaptmax);
    err_s_iter_z = nan(1,adaptmax);
    error_iter_z = nan(1,adaptmax);
    error_iter   = nan(1,adaptmax);
    Goal_unc_iter = nan(1,adaptmax);
    B_iter        = nan(1,adaptmax);
    Goal_iter     = nan(1,adaptmax);
    
    err_p_dir_iter_u = nan(1,adaptmax);
    err_s_dir_iter_u = nan(1,adaptmax);
    error_dir_iter_u = nan(1,adaptmax);
    err_p_dir_iter_z = nan(1,adaptmax);
    err_s_dir_iter_z = nan(1,adaptmax);
    error_dir_iter_z = nan(1,adaptmax);
    error_dir_iter   = nan(1,adaptmax);
    error_dir_iter0  = nan(1,adaptmax);

% preallocation of cell data
    sols_u_iter = cell(1,adaptmax);
    sols_z_iter = cell(1,adaptmax);
    paras_sg_iter = cell(1,adaptmax);
    paras_fem_iter = cell(1,adaptmax);
    
iter = 0; glevel = 0;
startLoopTime = tic;
while tot_err_direct >= delta && iter <= adaptmax
    if iter == 0 % First iteration step
        % Initial index set is for a single collocation point
        X = stochcol_getindexset(0, M);
        % Specimen Index Set
        % X = [1 1 1 1; 2 1 1 1];
        
        % Several attributes of the general sparse grid interpolation given
        % the index set and the sequence of 1D collocation nodes, such as
        % the underlying grid, the coordinates of the underlying grid, the
        % multivariate Lagrange polynomials (represented by the coefficient
        % list of the single-term multivariate Lagrange polynomials and the
        % set of indexing the degree of the single-term polynomials), etc.
        paras_sg = stochcol_sg([X, ones(size(X, 1), Q)], rule_id);
        gridd = paras_sg{4};
        clincombiset = paras_sg{5};
        indlincombiset = paras_sg{6};
        coords = paras_sg{9};
        input(1) = input(1) + Q;
        if ~exist('tout','var') 
            tout = 0; 
        end
        if tout == 1
            figure(1)
            scatter(coords(:,1),coords(:,2),100,'o','filled');
            axis square
            box on
            grid on
            title('initial collocation nodes in first two directions')
        end
        % Several attributes of the finite element: vertex coordinate
        % vector, element mapping matrix, boundary vertex vector, boundary
        % elment mapping matrix
        paras_fem0 = stochcol_fem_grid_generator(pmethod, dom_paras);
        paras_fem = paras_fem0;
        % Collection of edge information for flux jump computation and 
        % outputs of linear detail space Y grid-generator
        paras_detail0 = stochcol_mesh_detail(paras_fem0);
        paras_detail = paras_detail0;
        % plot the first mesh
        if tout == 1
            if pmethod == 1 % P1
                plot_mesh(paras_fem0{2}, paras_fem0{1}, 'initial FE mesh');
                pause(1)
            elseif pmethod == 2
                plot_mesh(paras_fem0{6}, paras_fem0{5}, 'initial FE mesh');
            end
        end
        % We represent the general sparse grid interpolation of the finite
        % element approximation as a cell of vertex values associated with
        % collocation nodes in the order they are introduced
        A = cell(1,size(coords, 1));
        sols_u = zeros(length(paras_fem0{1}), size(coords, 1));
        sols_z = zeros(length(paras_fem0{1}), size(coords, 1));
        elerrs_u = zeros(length(paras_fem0{2}), size(coords, 1));
        elerrs_z = zeros(length(paras_fem0{2}), size(coords, 1));
        ederrs_u = zeros(length(paras_detail0{7}), size(coords, 1));
        ederrs_z = zeros(length(paras_detail0{7}), size(coords, 1));
        % For each collocation node, we calculate an estimate of the
        % energy error of the FE solution and a list of
        % marked elements and edges
        errests_u = nan(1,size(coords, 1));
        errests_z = nan(1,size(coords, 1));
        Msetsu  = cell(1,size(coords, 1));
        Msetsz  = cell(1,size(coords, 1));
        % use parfor to enable parallel computing
        
        parfor k = 1:size(coords, 1)%parfor
            % FE solution for a current collocation node
            if nargin(F_rhs{1}) == 3 
                F_rhs_tmp = F_rhs;
                F_rhs_tmp{1} = @(x1, x2) F_rhs{1}(x1, x2, coords(k, :));
            else
                F_rhs_tmp = F_rhs;
            end
            [u_gal, z_gal, Anbc] = goafem_stochcol_fem_solver(coords(k, :), paras_fem0, ...
                aa, F_rhs_tmp, G_rhs);
            % Energy error of the FE solution and marked elements and edges
            [errest_u, errest_z, elerr_u, elerr_z, ederr_u, ...
                    ederr_z, Msetu, Msetz] = ...
                    goafem_stochcol_fem_estimator(pmethod, paras_fem_errest, ...
                        paras_fem0, paras_detail0, u_gal, z_gal, ...
                        coords(k, :), aa, aax1, aax2, F_rhs_tmp, G_rhs);
            sols_u(:, k) = u_gal;            
            sols_z(:, k) = z_gal;
            sols_u_ml{k} = u_gal;
            sols_z_ml{k} = z_gal;
            if rv_id == 3
                mesh_x = paras_fem{1}(:,1);
                mesh_y = paras_fem{1}(:,2);
                a_coord = aa(mesh_x, mesh_y, coords(k, :));
                a_min_z = min(a_coord);
                errests_u(k) = errest_u/a_min_z;
                errests_z(k) = errest_z/a_min_z;
            else
                errests_u(k) = errest_u;
                errests_z(k) = errest_z;
            end
            elerrs_u(:, k) = elerr_u;
            elerrs_z(:, k) = elerr_z;
            ederrs_u(:, k) = ederr_u;
            ederrs_z(:, k) = ederr_z;
            Msetsu{k} = Msetu;
            Msetsz{k} = Msetz;
            A{k} = Anbc;
            Gu(k) = u_gal' * Anbc * z_gal;
        end %end parfor k
        
        % 2-norm of multivariate Lagrange polynomials
        L_two_norm = stochcol_multilag_Ltwonorm(paras_sg, list);
        % construct the detail index set ---- reduced margin
        X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X), Q);
        % X_diff = stochcol_margin_set(X);
        % Attributes of the general sparse grid interpolation corresponding
        % to the new index set, which is the union of the original index 
        % set and the detail index set
        paras_sg_full = stochcol_sg([X, ones(size(X, 1), Q); X_diff], rule_id);
        % set grid due to the inclusion of the detail index set
        % the grid is used to label the collocation nodes
        [gridd_diff, grid_diff_ind] = setdiff(paras_sg_full{4}, ...
            paras_sg{4}, 'rows');
        coords_diff = paras_sg_full{9}(grid_diff_ind,:);
        sols_u0 = sols_u;
        sols_z0 = sols_z;
        close all
    else %iter > 0
        fprintf('\n\nIteration %i \n',iter)
        fprintf(['Primal spatial error indicator    ', ...
         'is %10.4e \n'],serrest_u)
        fprintf(['Primal parametric error indicator ', ...
         'is %10.4e \n'],perrest_u)
        fprintf(['Dual spatial error indicator      ', ...
         'is %10.4e \n'],serrest_z)
        fprintf(['Dual parametric error indicator   ', ...
         'is %10.4e \n'],perrest_z)
        if spind == 1
            fprintf('Parametric enrichment step ... new indices added \n');
            % Determine marked index set
            
            if new_dim_ind == 0 || Q == 0
                [Mset] = ...
                    goafem_mark_combine(Mset_u(:, 1:M), Mset_z(:, 1:M), markedgelem, 3, ...
                                        perrests_u(IE)', perrests_u(IF)'+perrests_z(IF)');
               
                disp(Mset)
                X = stochcol_getgrid([X; Mset]);
            else
                new_dim = find(Mset - 1);
                %M = M + 1;
                input(1) = new_dim;
                disp(Mset(:, 1:new_dim))
                fprintf('New dimensionality is %d \n', new_dim)
                paras_sg = stochcol_sg([X, ones(size(X,1), new_dim-M+Q)], rule_id);
                gridd_diff = [gridd_diff, ones(size(gridd_diff, 1), new_dim-M)];
                X_tmp = [X, ones(size(X, 1), new_dim-M)];
                for q=1:(new_dim-M)
                    new_row = ones(1, size(X_tmp, 2));
                    new_row(end-q+1) = 2;
                    X_tmp = [X_tmp; new_row];
                end
                X = stochcol_getgrid(X_tmp);
                M = new_dim;
            end
            
            
            glevel = glevel+1;
            % Attributes of the general sparse grid interpolation based on
            % the refined index set
            paras_sg_new = stochcol_sg([X, ones(size(X,1), Q)], rule_id);
            % Grid of collocation nodes based on previous index set
            gridd = paras_sg{4}; 
            % Grid and coordinates of collocation nodes based on the
            % refined index set
            gridd_new = paras_sg_new{4};
            coords_new = paras_sg_new{9};
            % FE solutions, energy error estimates, marked elements and
            % edges for collocation nodes  based on the refined index set
            % All FE solutions and part of estimates are from last iteration
            A_new = cell(1,size(coords_new, 1));
            sols_u_new = zeros(length(paras_fem{1}), size(coords_new, 1));
            sols_z_new = zeros(length(paras_fem{1}), size(coords_new, 1));
%             sols_u_new0 = zeros(length(paras_fem0{1}), size(coords_new, 1));
%             sols_z_new0 = zeros(length(paras_fem0{1}), size(coords_new, 1));
            elerrs_u_new = zeros(length(paras_fem{2}), size(coords_new, 1));
            elerrs_z_new = zeros(length(paras_fem{2}), size(coords_new, 1));
            ederrs_u_new = zeros(length(paras_detail{7}), size(coords_new, 1));
            ederrs_z_new = zeros(length(paras_detail{7}), size(coords_new, 1));
            errests_u_new = zeros(1,size(coords_new, 1));
            errests_z_new = zeros(1,size(coords_new, 1));
            Msetsu_new  = cell(1,size(coords_new, 1));
            Msetsz_new  = cell(1,size(coords_new, 1));
            [~, IA, IB] = intersect(gridd_new, gridd, 'rows');
            [~, IC, ID] = intersect(gridd_new, gridd_diff, 'rows');
            for k = 1:length(IA)
                % reuse the FE solutions and estimates from the last iteration
                A_new{IA(k)} = A{IB(k)};
                sols_u_new(:, IA(k)) = sols_u(:, IB(k));
                sols_z_new(:, IA(k)) = sols_z(:, IB(k));
%                 sols_u_new0(:, IA(k)) = sols_u0(:, IB(k));
%                 sols_z_new0(:, IA(k)) = sols_z0(:, IB(k));
                elerrs_u_new(:, IA(k)) = elerrs_u(:, IB(k));
                elerrs_z_new(:, IA(k)) = elerrs_z(:, IB(k));
                ederrs_u_new(:, IA(k)) = ederrs_u(:, IB(k));
                ederrs_z_new(:, IA(k)) = ederrs_z(:, IB(k));
                errests_u_new(IA(k)) = errests_u(IB(k));
                errests_z_new(IA(k)) = errests_z(IB(k));
                Msetsu_new{IA(k)} = Msetsu{IB(k)};
                Msetsz_new{IA(k)} = Msetsz{IB(k)};
            end
            for k = 1:length(IC)
                % reuse the FE solutions from the last iteration
                A_new{IC(k)} = A_diff{ID(k)};
                sols_u_new(:, IC(k)) = sols_u_diff(:, ID(k));
                sols_z_new(:, IC(k)) = sols_z_diff(:, ID(k));
%                 sols_u_new0(:, IC(k)) = sols_u_diff0(:, ID(k));
%                 sols_z_new0(:, IC(k)) = sols_z_diff0(:, ID(k));
            end
            % Compute new FE error estimates
            coords_temp = coords_diff(ID, :);
            sols_u_diff_temp = sols_u_diff(:, ID);
            sols_z_diff_temp = sols_z_diff(:, ID);
            elerrs_u_temp = zeros(length(paras_fem{2}), length(ID));
            elerrs_z_temp = zeros(length(paras_fem{2}), length(ID));
            ederrs_u_temp = zeros(length(paras_detail{7}), length(ID));
            ederrs_z_temp = zeros(length(paras_detail{7}), length(ID));
            errests_u_temp = zeros(1,length(ID));
            errests_z_temp = zeros(1,length(ID));
            Msetsu_temp = cell(1,length(ID));
            Msetsz_temp = cell(1,length(ID));
            
            parfor k = 1:length(ID) %parfor
                u_gal = sols_u_diff_temp(:,k);
                z_gal = sols_z_diff_temp(:,k);
                if nargin(F_rhs{1}) == 3 
                    F_rhs_tmp = F_rhs;
                    F_rhs_tmp{1} = @(x1, x2) F_rhs{1}(x1, x2, coords_temp(k, :));
                else
                    F_rhs_tmp = F_rhs;
                end
            [errest_u, errest_z, elerr_u, elerr_z, ederr_u, ...
                ederr_z, Msetu, Msetz] = ...
             goafem_stochcol_fem_estimator(pmethod, paras_fem_errest, ...
                    paras_fem, paras_detail, u_gal, z_gal, ...
                    coords_temp(k, :), aa, aax1, aax2, F_rhs_tmp, G_rhs);
                errests_u_temp(k) = errest_u;
                errests_z_temp(k) = errest_z;
                elerrs_u_temp(:, k) = elerr_u;
                elerrs_z_temp(:, k) = elerr_z;
                ederrs_u_temp(:, k) = ederr_u;
                ederrs_z_temp(:, k) = ederr_z;
                Msetsu_temp{k} = Msetu;
                Msetsz_temp{k} = Msetz;
            end % end parfor
            for k = 1:length(IC)
                elerrs_u_new(:,IC(k)) = elerrs_u_temp(:, k);
                elerrs_z_new(:,IC(k)) = elerrs_z_temp(:, k);
                ederrs_u_new(:,IC(k)) = ederrs_u_temp(:, k);
                ederrs_z_new(:,IC(k)) = ederrs_z_temp(:, k);
                errests_u_new(IC(k)) = errests_u_temp(k);
                errests_z_new(IC(k)) = errests_z_temp(k);
                Msetsu_new{IC(k)} = Msetsu_temp{k};
                Msetsz_new{IC(k)} = Msetsz_temp{k};
            end % end for
            paras_sg = paras_sg_new;
            gridd = paras_sg{4};
            clincombiset = paras_sg{5};
            indlincombiset = paras_sg{6};
            coords = paras_sg{9};
                        
            A = A_new;
            sols_u = sols_u_new;
            sols_z = sols_z_new;
%             sols_u0 = sols_u_new0;
%             sols_z0 = sols_z_new0;            
            elerrs_u = elerrs_u_new;
            elerrs_z = elerrs_z_new;
            ederrs_u = ederrs_u_new;
            ederrs_z = ederrs_z_new;
            errests_u = errests_u_new;
            errests_z = errests_z_new;
            Msetsu = Msetsu_new;
            Msetsz = Msetsz_new;
            for k = 1:size(A,2)
                Gu(k) = sols_u(:,k)' * A{k} * sols_z(:,k);
            end
            % 2-norm of multivariate Lagrange polynomials
            L_two_norm = stochcol_multilag_Ltwonorm(paras_sg, list);
            % construct the detail index set ---- reduced margin
            X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X), Q);
%            X_diff = stochcol_margin_set(X);
            % Attributes of the general sparse grid interpolation corresponding
            % to the new index set, which is the union of the original index
            % set and the detail index set
            paras_sg_full = stochcol_sg([X, ones(size(X,1), Q); X_diff], rule_id);
            % set grid due to the inclusion of the detail index set
            % the grid is used to label the collocation nodes
            [gridd_diff, grid_diff_ind] = setdiff(paras_sg_full{4}, ...
                paras_sg{4}, 'rows');
            coords_diff = paras_sg_full{9}(grid_diff_ind,:);
        else % spind == 0
            fprintf('Spatial refinement step...\n');
            % Mesh refinement
            [MMele, MMedge] = ...
                   goafem_mark_combine(Msetu, Msetz, markedgelem, 3, ... 
                                       elerrs_u, elerrs_u + elerrs_z, ... 
                                       ederrs_u, ederrs_u + ederrs_z, paras_detail); 
       
            fprintf('original number of elements is %g\n',length(paras_fem{2}(:,1)));
            paras_fem = stochcol_mesh_refine(MMele, MMedge, paras_fem, ...
                                             paras_detail, pmethod);
            fprintf('new number of elements      is %g\n\n',length(paras_fem{2}(:,1)));
            % Collection of edge information for flux jump computation and
            % outputs of linear detail space Y grid-generator based on
            % refined mesh
            paras_detail = stochcol_mesh_detail(paras_fem);
            A = cell(1, size(coords,1));
            sols_u = zeros(length(paras_fem{1}), size(coords, 1));
            sols_z = zeros(length(paras_fem{1}), size(coords, 1));
%             sols_u0 = zeros(length(paras_fem0{1}), size(coords, 1));
%             sols_z0 = zeros(length(paras_fem0{1}), size(coords, 1));
            elerrs_u = zeros(length(paras_fem{2}), size(coords, 1));
            elerrs_z = zeros(length(paras_fem{2}), size(coords, 1));
            ederrs_u = zeros(length(paras_detail{7}), size(coords, 1));
            ederrs_z = zeros(length(paras_detail{7}), size(coords, 1));
            errests_u = zeros(1,size(coords, 1));
            errests_z = zeros(1,size(coords, 1));
            Msetsu  = cell(1,size(coords, 1));
            Msetsz  = cell(1,size(coords, 1));
            for k = 1:size(coords, 1)
                if nargin(F_rhs{1}) == 3 
                    F_rhs_tmp = F_rhs;
                    F_rhs_tmp{1} = @(x1, x2) F_rhs{1}(x1, x2, coords(k, :));
                else
                    F_rhs_tmp = F_rhs;
                end
                [u_gal, z_gal, Anbc] = goafem_stochcol_fem_solver( ...
                    coords(k, :), paras_fem, aa, F_rhs_tmp, G_rhs);
%                 [u_gal0, z_gal0, ~] = goafem_stochcol_fem_solver( ...
%                     coords(k, :), paras_fem0, aa, F_rhs, G_rhs);
                [errest_u, errest_z, elerr_u, elerr_z, ederr_u, ...
                    ederr_z, Msetu, Msetz] = ...
                goafem_stochcol_fem_estimator(pmethod, paras_fem_errest, ...
                    paras_fem, paras_detail, u_gal, z_gal, ...
                    coords(k, :), aa, aax1, aax2, F_rhs_tmp, G_rhs);
                sols_u(:, k) = u_gal;
                sols_z(:, k) = z_gal;
%                 sols_u0(:, k) = u_gal0;
%                 sols_z0(:, k) = z_gal0;
                elerrs_u(:, k) = elerr_u;
                elerrs_z(:, k) = elerr_z;
                ederrs_u(:, k) = ederr_u;
                ederrs_z(:, k) = ederr_z;
                if rv_id == 3
                    mesh_x = paras_fem{1}(:,1);
                    mesh_y = paras_fem{1}(:,2);
                    a_coord = aa(mesh_x, mesh_y, coords(k, :));
                    a_min_z = min(a_coord);
                    errests_u(k) = errest_u/a_min_z;
                    errests_z(k) = errest_z/a_min_z;
                else
                    errests_u(k) = errest_u;
                    errests_z(k) = errest_z;
                end
                Msetsu{k} = Msetu;
                Msetsz{k} = Msetz;
                A{k} = Anbc;
                Gu(k) = u_gal' * Anbc * z_gal;
            end % end for
        end
    end % if iteration no.
    
    % spatial error indicator (assembled from components)
    serrest_u = errests_u * L_two_norm;
    serrest_z = errests_z * L_two_norm;
    
    % generate union of marked elements or edges (single grid code)
    Msetu  = munion(Msetsu);
    Msetz  = munion(Msetsz);

    % parametric error indicator
    % compute interpolated FE approximations corresponding to new nodes
    A_diff = cell(1,size(coords_diff,1));
    sols_u_diff = zeros(length(paras_fem{1}), size(coords_diff, 1));
    sols_z_diff = zeros(length(paras_fem{1}), size(coords_diff, 1));
%     sols_u_diff0 = zeros(length(paras_fem0{1}), size(coords_diff, 1));
%     sols_z_diff0 = zeros(length(paras_fem0{1}), size(coords_diff, 1));
    errest2s_u = zeros(1, size(coords_diff, 1));
    errest2s_z = zeros(1, size(coords_diff, 1));
    a_unit = @(x1, x2) ones(size(x1)); 
    if nargin(F_rhs{1}) == 3 
            F_rhs_tmp = F_rhs;
            F_rhs_tmp{1} = @(x1, x2) F_rhs{1}(x1, x2, coords(1, :));
    else
            F_rhs_tmp = F_rhs;
    end
    [A_unit,~,~] = goafem_stochcol_fem_setup(paras_fem{1}, ...
                   paras_fem{2}, a_unit, F_rhs_tmp);
    parfor k = 1:size(coords_diff, 1)%parfor
        if nargin(F_rhs{1}) == 3 
            F_rhs_tmp = F_rhs;
            F_rhs_tmp{1} = @(x1, x2) F_rhs{1}(x1, x2, coords_diff(k, :));
        else
            F_rhs_tmp = F_rhs;
        end
        % FE solution for a propective collocation node
        [u_gal, z_gal, Anbc] = goafem_stochcol_fem_solver(coords_diff(k, :), ...
                               paras_fem, aa, F_rhs_tmp, G_rhs);
        A_diff{k} = Anbc;                   
        sols_u_diff(:, k) = u_gal;
        sols_z_diff(:, k) = z_gal;
        
%         [u_gal0, z_gal0, ~] = goafem_stochcol_fem_solver(coords_diff(k, :), ...
%                                paras_fem0, aa, F_rhs, G_rhs);
%         sols_u_diff0(:, k) = u_gal0;
%         sols_z_diff0(:, k) = z_gal0;
        
        % vector of multivariate Lagrange polynomials for the collocation node
        LL_diff = stochcol_getinterpolant_2(gridd, ...
                  clincombiset, indlincombiset, coords_diff(k,:), polys);
        % interpolated FE solution at the collocation node
        uyy = sols_u*LL_diff;
        zyy = sols_z*LL_diff;
        % H1 seminorm of the interpolation error at the collocation node
        errest2s_u(k) = sqrt((u_gal - uyy)' * A_unit * (u_gal - uyy));
        errest2s_z(k) = sqrt((z_gal - zyy)' * A_unit * (z_gal - zyy));
    end
    
    % LL_n = stochcol_getinterpolant_2(gridd, clincombiset, indlincombiset, coords_diff, polys)
    
    % indexwise parametric error estimates based on the detail collocation nodes
    [perrests_u, perrests_z] = ...
        goafem_stochcol_est_parametric(X_diff, errest2s_u, errest2s_z, ...
        gridd_diff, list, rule_id);
    % parametric error indicator
    perrest_u = sum(perrests_u);
    perrest_z = sum(perrests_z);    
    
    if Q ~= 0 %dimensionaly adaptive
        n = size(X_diff,1);
        costs = ones(1, n);
       for i = 1:n
            row = X_diff(i, :);
            for j=1:M+1
                if row(j)>1
                    if rule_id == 2
                        costs(i) = costs(i)*(2^(row(j)-1)+1);
                    elseif rule_id == 3
                        costs(i) = costs(i)*(2*row(j)-1);
                    else
                        costs(i) = costs(i)*row(j);
                    end
                end
            end
        end 
    
        % parametric error indicator for new dimension
        perrests_new_dim_u = perrests_u(1:Q)./costs(1:Q);
        perrests_new_dim_z = perrests_z(1:Q)./costs(1:Q);
        [new_dim_err_max, ind_max] = max(perrests_new_dim_u + perrests_new_dim_z);
        % marked index set from the detail index set
        %disp(X_diff)
        %disp(perrests)
            
        if  sum(perrests_new_dim_u + perrests_new_dim_z) > sum(perrests_u(Q+1:n)./costs(Q+1:n)) + sum(perrests_z(Q+1:n)./costs(Q+1:n))
            Mset = X_diff(ind_max, :);
            new_dim_ind = 1;
        else
            [Mset_u, ~] = dorfler_marking(X_diff(Q+1:n, :), perrests_u(Q+1:n), pmthreshold);
            [Mset_z, ~] = dorfler_marking(X_diff(Q+1:n, :), perrests_u(Q+1:n) +  perrests_z(Q+1:n), pmthreshold); 
            [~, ~, IE] = intersect(Mset_u, X_diff, 'rows');
            [~, ~, IF] = intersect(Mset_z, X_diff, 'rows');
            new_dim_ind = 0;
        end
    else
        [Mset_u, ~] = dorfler_marking(X_diff, perrests_u, pmthreshold);
        [Mset_z, ~] = dorfler_marking(X_diff, perrests_u + perrests_z, pmthreshold); 
        [~, ~, IE] = intersect(Mset_u, X_diff, 'rows');
        [~, ~, IF] = intersect(Mset_z, X_diff, 'rows');
        new_dim_ind = 0;
    end
    
    
    spind = goafem_indicator(serrest_u, serrest_z, perrest_u, perrest_z, 3);
     
    % total error indicator
    tot_err_est_u = serrest_u + perrest_u;
    tot_err_est_z = serrest_z + perrest_z;
    if G_rhs{5} >= 0
        tot_err_est = tot_err_est_u *(tot_err_est_u + tot_err_est_z);
    else
        tot_err_est = tot_err_est_u * tot_err_est_z;
    end
                           
    % compute direct estimates of error
    % spatial
    [serrestdir_u, serrestdir_z] = goafem_stochcol_direct_estimatorX(...
        paras_sg, paras_fem, list, pmethod, F_rhs, G_rhs, aa);
    % parametric
    paras_sg_diff = stochcol_sg(X_diff, rule_id);
    G_diff = stochcol_gmatrices(paras_sg_diff{4}, ...
            paras_sg_diff{5},paras_sg_diff{6}, list);
    ncpts = length(G_diff(:,1));
    ncptsf = size(coords, 1)+ size(coords_diff, 1);
    grid_old_ind = setdiff([1:ncptsf]',grid_diff_ind);
    % construct array containing all the cpoint solution vectors
    sols_u_all = nan(length(paras_fem{1}), ncptsf);
    sols_z_all = nan(length(paras_fem{1}), ncptsf);
%     sols_u_all0 = nan(length(paras_fem0{1}), ncptsf);
%     sols_z_all0 = nan(length(paras_fem0{1}), ncptsf);
    sols_u_all(:,grid_old_ind) = sols_u; sols_u_all(:,grid_diff_ind) = sols_u_diff;
    sols_z_all(:,grid_old_ind) = sols_z; sols_z_all(:,grid_diff_ind) = sols_z_diff;
%     sols_u_all0(:,grid_old_ind) = sols_u0; sols_u_all0(:,grid_diff_ind) = sols_u_diff0;
%     sols_z_all0(:,grid_old_ind) = sols_z0; sols_z_all0(:,grid_diff_ind) = sols_z_diff0;
        if ncptsf > ncpts
        fprintf('\nredundant sparse grid solutions ..\n')
        [gdiff,indx]=setdiff(paras_sg_full{4}, paras_sg_diff{4}, 'rows');
%        gdiff
        iactive=setdiff([1:ncptsf]',indx);
        sols_u_all = sols_u_all(:,iactive);
        sols_z_all = sols_z_all(:,iactive);
%         sols_u_all0 = sols_u_all0(:,iactive);
%         sols_z_all0 = sols_z_all0(:,iactive);
        elseif ncptsf < ncpts, error('Oops ..  fatal logic issue'),
        end
    % compute the energy norm of the SG solution via hierarchical surplus
    %[A_unit,~,~] = goafem_stochcol_fem_setup(paras_fem{1}, ...
    %               paras_fem{2}, a_unit, F_rhs);
    perrestdir_u = sqrt(sum(dot(G_diff, sols_u_all' * A_unit * sols_u_all)));
    perrestdir_z = sqrt(sum(dot(G_diff, sols_z_all' * A_unit * sols_z_all)));
    
    %perrestdir_u0 = sqrt(sum(dot(G_diff, sols_u_all0' * A_unit0 * sols_u_all0)));
    %perrestdir_z0 = sqrt(sum(dot(G_diff, sols_z_all0' * A_unit0 * sols_z_all0)));
    
    tot_err_direct_u = serrestdir_u + perrestdir_u;
    tot_err_direct_z = serrestdir_z + perrestdir_z;
    
    %tot_err_direct_u0 = serrestdir_u + perrestdir_u0;
    %tot_err_direct_z0 = serrestdir_z + perrestdir_z0;
    
    if G_rhs{5} >= 0
        tot_err_direct = tot_err_direct_u *( tot_err_direct_z +  tot_err_direct_u);
        %tot_err_direct0 = tot_err_direct_u0 *( tot_err_direct_z0 +  tot_err_direct_u0);
    else
        tot_err_direct = tot_err_direct_u * tot_err_direct_z;
        %tot_err_direct0 = tot_err_direct_u0 * tot_err_direct_z0;
    end
                             
    % output error estimates
     if iter==0, fprintf('\n\nIteration %i \n',iter), end
    fprintf(['primal spatial error estimate     ', ...
    'is %10.4e  vs  %10.4e (spatial indicator)   \n'],serrestdir_u,serrest_u)
    fprintf(['dual spatial error estimate       ', ...
    'is %10.4e  vs  %10.4e (spatial indicator)   \n'],serrestdir_z,serrest_z)
    fprintf(['primal parametric error estimate  ', ...
    'is %10.4e  vs  %10.4e (parametric indicator)\n'],perrestdir_u,perrest_u)
    fprintf(['dual parametric error estimate    ', ...
    'is %10.4e  vs  %10.4e (parametric indicator)\n'],perrestdir_z,perrest_z)
    fprintf(['overall estimate from indicators  ', ...
        'is %10.4e '],tot_err_est)
    fprintf(['\noverall direct error estimate     ', ...
    'is<strong> %10.4e </strong>\n'],tot_err_direct)

    [~,~,~,~,wp] = stochcol_multilag(paras_sg{4}, paras_sg{1}, paras_sg{2}, paras_sg{3}, rule_id, fun_p);
    % Goal Functional Calculation
    input_tmp = input;
    input_tmp(1) = size(paras_sg{9},2);
    KL_DATA = stoch_kl(input_tmp,expansion_type);
    if G_rhs{5} >= 0
        G = stochcol_gmatrices(paras_sg{4}, paras_sg{5}, paras_sg{6}, list);
        if G_rhs{5} == 0
            Mass = stochcol_mass_matrix(paras_fem{1}, paras_fem{2}, G_rhs{1});
            Goal_unc = sum(dot(G, sols_u' * Mass * sols_u));
        elseif G_rhs{5} == 1
            [~, Mass] = stochcol_mass_matrix(paras_fem{1}, paras_fem{2}, G_rhs{1});
            Goal_unc = sum(dot(G, sols_u' * Mass * sols_u));            
        elseif G_rhs{5} == 2
            [~, wz] = goafem_stochcol_fem_setup(paras_fem{1}, paras_fem{2}, a_unit, G_rhs);
            u_col_point = wz.'*sols_u;
            Goal_unc = u_col_point*G*u_col_point.';
        elseif G_rhs{5} == 3
            [~, wz] = goafem_stochcol_fem_setup(paras_fem{1}, paras_fem{2}, a_unit, G_rhs);
            u_col_point = wz.'*sols_u;
            %E = G_rhs{6};
            E = u_col_point*wp;
            %G_rhs{6} = E;
            Goal_unc = 100*(u_col_point*G*u_col_point.' - E^2);
        end
        B = goafem_doBilinearForm(KL_DATA, paras_fem, paras_sg, list, ...
                    listy, listy2, sols_u, sols_z, F_rhs, G_rhs, rf_type, ...
                    X, fun_p, rule_id, aa, polys, Q); 
        if nargin(F_rhs{1}) == 3
            F = goafem_compute_F(paras_fem, paras_sg, sols_z, F_rhs, G_rhs, ...
                    X, fun_p, rule_id, aa, polys);
        else
            [~, fz] = goafem_stochcol_fem_setup(paras_fem{1}, paras_fem{2}, a_unit, F_rhs);
            F = fz.'*sols_z*wp;
        end
        Goal = Goal_unc + F - B;
    else
        Goal = 0;
        for k = 1:size(coords,1)
            Goal = Goal + wp(k) * Gu(k);
        end
        Goal_unc = Goal;
        B = goafem_doBilinearForm(KL_DATA, paras_fem, paras_sg, list, ...
                listy, listy2, sols_u, sols_z, F_rhs, G_rhs, rf_type, ...
                X, fun_p, rule_id, aa, polys, Q);
        Goal = 2*Goal - B;
        
        %[~, fz] = goafem_stochcol_fem_setup(paras_fem{1}, paras_fem{2}, a_unit, F_rhs);
        %    F = fz.'*sols_z*wp;
    end
    
    % exponential code 
    iter = iter + 1;
    
    if pmethod == 1 && mod(iter, 5) == 0% P1
        TITLE = ['FE mesh at ', num2str(iter), ' iteration'];
        plot_mesh(paras_fem{2}, paras_fem{1}, TITLE);
        pause(1)
    end
    
    
    err_p_iter_u(iter) = perrest_u;
    err_s_iter_u(iter) = serrest_u;
    error_iter_u(iter) = tot_err_est_u;
    err_p_iter_z(iter) = perrest_z;
    err_s_iter_z(iter) = serrest_z;
    error_iter_z(iter) = tot_err_est_z; 
    error_iter(iter)   = tot_err_est;
        
    err_p_dir_iter_u(iter) = perrestdir_u;
    err_s_dir_iter_u(iter) = serrestdir_u;
    error_dir_iter_u(iter) = tot_err_direct_u;
    err_p_dir_iter_z(iter) = perrestdir_z;
    err_s_dir_iter_z(iter) = serrestdir_z;
    error_dir_iter_z(iter) = tot_err_direct_z;
    error_dir_iter(iter)   = tot_err_direct;
%    error_dir_iter0(iter)   = tot_err_direct0;

    sols_u_iter{iter}    = sols_u;
    sols_z_iter{iter}    = sols_z;
    B_iter(iter)         = B;
    Goal_unc_iter(iter)  = Goal_unc;
    fprintf('\nValue of uncorrected GF %.9f',Goal_unc);
    fprintf('\nValue of corrected GF %.9f\n',Goal);
    %fprintf('\n RHS %.9f\n',F);
    %fprintf('\n Bilinear form %.9f\n',B);
    Goal_iter(iter)      = Goal;
    paras_sg_iter{iter}  = paras_sg;
    paras_fem_iter{iter} = paras_fem;
    dof(iter) = size(sols_u,1)*size(sols_u,2);
    if nargin(F_rhs{1}) == 3
        fprintf('\nValue of uncorrected GF %.9f',Goal_unc/16);
        fprintf('\nValue of corrected GF %.9f\n',Goal/16);
        Goal_unc_iter(iter)  = Goal_unc/16;
        Goal_iter(iter)      = Goal/16;
    end

end % end while (adaptive loop)

% Resize quantities
goafem_postquantity;

endLoopTime = toc(startLoopTime);
%%
% plot error estimate
if ~exist('iplot','var')
    iplot = 1;
end
if iplot == 1
    goafem_convergenceplot
end

% final mesh
if pmethod == 1 % P1
    plot_mesh(paras_fem{2}, paras_fem{1}, 'Final FE mesh');
    pause(1)
% elseif pmethod == 2
%    plot_mesh(paras_fem{6}, paras_fem{5}, 'Final FE mesh');
end
% final nodes
coords = paras_sg{9};
figure(7)
scatter(coords(:,1),coords(:,2),100,'o','filled');
axis square
box on
grid on
title('final collocation nodes in first two directions')

% postprocess to generate statistics
G = stochcol_gmatrices(paras_sg{4}, paras_sg{5}, paras_sg{6}, list);
[~,~,~,~,wp] = stochcol_multilag(paras_sg{4}, paras_sg{1}, paras_sg{2}, paras_sg{3}, rule_id, fun_p);

MEAN_u = zeros(size(sols_u,1), 1);
MEAN_z = MEAN_u;
VAR_u  = MEAN_u;
VAR_z  = MEAN_u;
for k = 1:size(coords,1)
   MEAN_u = MEAN_u + wp(k) * sols_u(:,k);
   MEAN_z = MEAN_z + wp(k) * sols_z(:,k); 
end

for k = 1:size(coords,1)
    for l = 1:size(coords, 1)
        VAR_u = VAR_u + sols_u(:,k).*sols_u(:, l)*G(k,l); % second order moment
        VAR_z = VAR_z + sols_z(:,k).*sols_z(:, l)*G(k,l); 
    end
end
VAR_u = VAR_u - MEAN_u.^2;
VAR_z = VAR_z - MEAN_z.^2;

STD_u = VAR_u.^(1/2);
STD_z = VAR_z.^(1/2);

% plot solutions
if pmethod == 1
    goafem_sc_plotsol_p1(dom_type,MEAN_u,MEAN_z,VAR_u,VAR_z,paras_fem{2},paras_fem{1})
end

disp([paras_sg{4}(:, 1:M)])
fprintf('\n Final sparse grid\n')
fprintf('\n Tolerance reached in %g iterations',iter-1)
fprintf('\n    after %g parametric refinements',glevel)
fprintf('\n              Primal Mean maximum %9.6f\n',max(MEAN_u))
fprintf('          Primal Variance maximum %9.6f\n',max(VAR_u))
fprintf('Primal Standard Deviation maximum %9.6f\n',max(STD_u))
fprintf('\n              Dual Mean maximum %9.6f\n',max(MEAN_z))
fprintf('          Dual Variance maximum %9.6f\n',max(VAR_z))
fprintf('Dual Standard Deviation maximum %9.6f\n',max(STD_z))
fprintf('\n<strong>Elapsed time:</strong> %.2f sec\n',endLoopTime);

%goafem_referenceSC