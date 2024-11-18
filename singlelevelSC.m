%SINGLELEVELSC solves stochastic diffusion problem using adaptive single-level SC-FEM
%
% Main driver for running single-level adaptive stochastic collocation algorithm
%
% Latest update: AS; 28 June 2024
% Copyright (c) 2022 A. Bespalov, A. Savinov, D. Silvester, F. Xu

clear all; close all;
%%%%%%%%%%
addpath('goafem');
x1_0 = default('Domain x1-coord lower limit (default is 0.25)', 0.25);
        x1_1 = default('Domain x1-coord upper limit (default is 0.75)', 0.75);
        x2_0 = default('Domain x2-coord lower limit (default is 0.25)', 0.25);
        x2_1 = default('Domain x2-coord upper limit (default is 0.75)', 0.75);
        supp = @(x1,x2) 1.*(x1_0 < x1 & x1 < x1_1 & x2_0 < x2 & x2 < x2_1);
        area = (x1_1 - x1_0)*(x2_1 - x2_0);
        G_rhs = {@(x1,x2) (1/area).*supp(x1,x2), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), @(x1,x2) zeros(size(x1)), ...
                 -1};
%%%%%%%%%
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
% preallocation of arrays
    dof = nan(1,adaptmax);
    err_p_iter = nan(1,adaptmax);
    err_s_iter = nan(1,adaptmax);
    error_iter = nan(1,adaptmax);
    
    err_p_d_iter = nan(1,adaptmax);
    err_s_d_iter = nan(1,adaptmax);
    error_d_iter = nan(1,adaptmax);

% preallocation of cell data
    sols_iter = cell(1,adaptmax);
    paras_sg_iter = cell(1,adaptmax);
    paras_fem_iter = cell(1,adaptmax);
    Goal_unc_iter = nan(1,adaptmax);
    
iter = 0; glevel = 0;
startLoopTime = tic;
while tot_err_direct >= delta && iter <= adaptmax
    if iter == 0 % First iteration step
        % Initial index set is for a single collocation point
        X = stochcol_getindexset(0, M);
        % Specimen Index Set
        % X = [1 1 1 1; 2 1 1 1; 1 2 1 1];
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
        simulationParams(1) = simulationParams(1) + Q;
        %figure(901)
        %scatter(coords(:,1),coords(:,2),100,'o','filled');
        %axis square
        %box on
        %grid on
        %title('initial collocation nodes in first two directions')

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
        if pmethod == 1 % P1
            plot_mesh(paras_fem0{2}, paras_fem0{1}, 'initial FE mesh');
            pause(1)
        elseif pmethod == 2 % P2
            plot_mesh(paras_fem0{6}, paras_fem0{5}, 'initial FE mesh');
        end
        % We represent the general sparse grid interpolation of the finite
        % element approximation as a cell of vertex values associated with
        % collocation nodes in the order they are introduced
        sols = zeros(length(paras_fem0{1}), size(coords, 1));
        % For each collocation node, we calculate an estimate of the
        % energy error of the FE solution and a list of
        % marked elements and edges
        errests = nan(1,size(coords, 1));
        MMeles  = cell(1,size(coords, 1));
        MMedges  = cell(1,size(coords, 1));
        % use parfor to enable parallel computing
        for k = 1:size(coords, 1) 
            if nargin(rhs) == 3 
                rhs_fun = @(x1, x2) rhs(x1, x2, coords(k, :));
            else
                rhs_fun = rhs;
            end
            % FE solution for a current collocation node
            [x_gal, ~] = stochcol_fem_solver(coords(k, :), paras_fem0, ...
                aa, rhs_fun);
            % Energy error of the FE solution and marked elements and edges
            [errest, ~, MMele, MMedge] = stochcol_fem_estimator(...
                pmethod, paras_fem_errest, paras_fem0, paras_detail0, ...
                x_gal, coords(k, :), aa, aax1, aax2, rhs_fun);
            sols(:, k) = x_gal;
            if rv_id == 3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
                mesh_x = paras_fem{1}(:,1);
                mesh_y = paras_fem{1}(:,2);
                a_coord = aa(mesh_x, mesh_y, coords(k, :));
                a_min_z = min(a_coord);
                errests(k) = errest/a_min_z;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                errests(k) = errest;
            end
            MMeles{k} = MMele;
            MMedges{k} = MMedge;
        end
        sols0 = sols;
        % 2-norm of multivariate Lagrange polynomials
        L_two_norm = stochcol_multilag_Ltwonorm(paras_sg, list);
        % construct the detail index set ---- reduced margin
        X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X), Q);
        % Attributes of the general sparse grid interpolation corresponding
        % to the new index set, which is the union of the original index 
        % set and the detail index set
        paras_sg_full = stochcol_sg([X, ones(size(X, 1), Q); X_diff], rule_id);
        % set grid due to the inclusion of the detail index set
        % the grid is used to label the collocation nodes
        [gridd_diff, grid_diff_ind] = setdiff(paras_sg_full{4}, paras_sg{4}, 'rows');
        coords_diff = paras_sg_full{9}(grid_diff_ind,:);
    else
        fprintf('\n\nIteration %i \n',iter)
        fprintf(['   spatial error indicator', ...
         ' is %10.4e \n'],serrest)
        fprintf(['parametric error indicator', ...
         ' is %10.4e \n'],perrest)
        if serrest < perrest % serrest_d < perrest_d
            fprintf('Parametric refinement ... new indices added \n');
            glevel = glevel+1;
            % Refined index set
            %disp(Mset)
            if all(Mset(:, M+1:M+Q) == 1, 'all') || Q == 0
                X = X(:, 1:M);
                Mset = Mset(:, 1:M);
                disp(Mset)
                X = stochcol_getgrid([X; Mset]);
            else
                new_dim = find(Mset - 1);
                %M = M + 1;
                simulationParams(1) = new_dim;
                disp(Mset(1:new_dim))
                fprintf('New dimensionality is %d \n', new_dim)
                paras_sg = stochcol_sg([X, ones(size(X,1), new_dim-M+Q)], rule_id);
                gridd_diff = [gridd_diff, ones(size(gridd_diff, 1), new_dim-M)];
                X = stochcol_getgrid([X, ones(size(X, 1), new_dim-M); Mset(1:new_dim)]);
                M = new_dim;
            end
            
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
            sols_new = zeros(length(paras_fem{1}), size(coords_new, 1));
            sols_new0 = zeros(length(paras_fem0{1}), size(coords_new, 1));
            errests_new = zeros(1,size(coords_new, 1));
            MMeles_new  = cell(1,size(coords_new, 1));
            MMedges_new  = cell(1,size(coords_new, 1));
            [~, IA, IB] = intersect(gridd_new, gridd, 'rows');
            [~, IC, ID] = intersect(gridd_new, gridd_diff, 'rows');
            for k = 1:length(IA)
                % reuse the FE solutions and estimates from the last iteration
                sols_new(:, IA(k)) = sols(:, IB(k));
                sols_new0(:, IA(k)) = sols0(:, IB(k));
                errests_new(IA(k)) = errests(IB(k));
                MMeles_new{IA(k)} = MMeles{IB(k)};
                MMedges_new{IA(k)} = MMedges{IB(k)};
            end
            for k = 1:length(IC)
                % reuse the FE solutions from the last iteration
                sols_new(:, IC(k)) = sols_diff(:, ID(k));
                sols_new0(:, IC(k)) = sols_diff0(:, ID(k));
            end
            % Compute new FE error estimates
            coords_temp = coords_diff(ID, :);
            sols_diff_temp = sols_diff(:, ID);
            errests_temp = zeros(1,length(ID));
            MMeles_temp = cell(1,length(ID));
            MMedges_temp = cell(1,length(ID));
            for k = 1:length(ID)
                x_gal = sols_diff_temp(:,k);
                if nargin(rhs) == 3 
                    rhs_fun = @(x1, x2) rhs(x1, x2, coords_temp(k, :));
                else
                    rhs_fun = rhs;
                end

                [errest, ~, MMele, MMedge] = stochcol_fem_estimator(pmethod, ...
                    paras_fem_errest, paras_fem, paras_detail, x_gal, coords_temp(k, :), ...
                    aa, aax1, aax2, rhs_fun);
                errests_temp(k) = errest;
                MMeles_temp{k} = MMele;
                MMedges_temp{k} = MMedge;
            end
            for k = 1:length(IC)
                errests_new(IC(k)) = errests_temp(k);
                MMeles_new{IC(k)} = MMeles_temp{k};
                MMedges_new{IC(k)} = MMedges_temp{k};
            end
            paras_sg = paras_sg_new;
            gridd = paras_sg{4};
            clincombiset = paras_sg{5};
            indlincombiset = paras_sg{6};
            coords = paras_sg{9};
            sols = sols_new;
            sols0 = sols_new0;
            errests = errests_new;
            MMeles = MMeles_new;
            MMedges = MMedges_new;
            % 2-norm of multivariate Lagrange polynomials
            L_two_norm = stochcol_multilag_Ltwonorm(paras_sg, list);
            % construct the detail index set ---- reduced margin
            X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X), Q);
            % Attributes of the general sparse grid interpolation corresponding
            % to the new index set, which is the union of the original index
            % set and the detail index set
            paras_sg_full = stochcol_sg([X, ones(size(X,1), Q); X_diff], rule_id);
            % set grid due to the inclusion of the detail index set
            % the grid is used to label the collocation nodes
            [gridd_diff, grid_diff_ind] = setdiff(paras_sg_full{4}, paras_sg{4}, 'rows');
            coords_diff = paras_sg_full{9}(grid_diff_ind,:);
        else
            fprintf('Spatial refinement...\n');
            % Mesh refinement
            fprintf('original number of elements is %g\n',length(paras_fem{2}(:,1)));
            paras_fem = stochcol_mesh_refine(MMele, MMedge, paras_fem, paras_detail, pmethod);
            fprintf('     new number of elements is %g\n\n',length(paras_fem{2}(:,1)));
            % Collection of edge information for flux jump computation and
            % outputs of linear detail space Y grid-generator based on
            % refined mesh
            paras_detail = stochcol_mesh_detail(paras_fem);
            sols = zeros(length(paras_fem{1}), size(coords, 1));
            sols0 = zeros(length(paras_fem0{1}), size(coords, 1));
            errests = zeros(1,size(coords, 1));
            MMeles  = cell(1,size(coords, 1));
            MMedges  = cell(1,size(coords, 1));
            for k = 1:size(coords, 1)
                if nargin(rhs) == 3 
                    rhs_fun = @(x1, x2) rhs(x1, x2, coords(k, :));
                else
                    rhs_fun = rhs;
                end
                [x_gal, ~] = stochcol_fem_solver(coords(k, :), ...
                    paras_fem, aa, rhs_fun);
                [x_gal0, ~] = stochcol_fem_solver(coords(k, :), ...
                    paras_fem0, aa, rhs_fun);
                [errest, ~, MMele, MMedge] = stochcol_fem_estimator(pmethod, ...
                    paras_fem_errest, paras_fem, paras_detail, x_gal, coords(k, :), ...
                    aa, aax1, aax2, rhs_fun);
                sols(:, k) = x_gal;
                sols0(:, k) = x_gal0;
                if rv_id == 3
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    mesh_x = paras_fem{1}(:,1);
                    mesh_y = paras_fem{1}(:,2);
                    a_coord = aa(mesh_x, mesh_y, coords(k, :));
                    a_min_z = min(a_coord);
                    errests(k) = errest/a_min_z;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    errests(k) = errest;
                end
                MMeles{k} = MMele;
                MMedges{k} = MMedge;
            end
            %Plotting statistics at selected iterations
            %if nargin(rhs) == 3 && mod(iter, 5) == 0
            %    plotstats(M, paras_sg,coords,fun_p,paras_fem,sols, dom_type), pause(.1);
            %end
        end
    end %------------- of adaptive loop
    

    % spatial error indicator (assembled from components)
    serrest = errests*L_two_norm;
    % generate union of marked elements and edges (single grid code)
    MMele  = munion(MMeles);
    MMedge = munion(MMedges);

    % parametric error indicator
    % compute interpolated FE approximations corresponding to the new
    % collocation nodes
    sols_diff = zeros(length(paras_fem{1}), size(coords_diff, 1));
    sols_diff0 = zeros(length(paras_fem0{1}), size(coords_diff, 1));
    errest2s = zeros(1, size(coords_diff, 1));
    a_unit = @(x1, x2) ones(size(x1));
    rhs_unit = @(x1, x2) ones(size(x1));
    [A_unit,~] = stochcol_fem_setup(paras_fem{1}, paras_fem{2}, a_unit, rhs_unit);
    for k = 1:size(coords_diff, 1)
        % FE solution for a propective collocation node
        if nargin(rhs) == 3 
            rhs_fun = @(x1, x2) rhs(x1, x2, coords_diff(k, :));
        else
            rhs_fun = rhs;
        end
        [x_gal, ~] = stochcol_fem_solver(coords_diff(k, :), ...
            paras_fem, aa, rhs_fun);
        [x_gal0, ~] = stochcol_fem_solver(coords_diff(k, :), ...
            paras_fem0, aa, rhs_fun);
        sols_diff(:, k) = x_gal;
        sols_diff0(:, k) = x_gal0;
        pause(.1);
        % vector of multivariate Lagrange polynomials for the collocation node
        LL_diff = stochcol_getinterpolant_2(gridd, ...
            clincombiset, indlincombiset, coords_diff(k,:), polys);
        % interpolated FE solution at the collocation node
        uyy = sols*LL_diff;
        % H1 seminorm of the interpolation error at the collocation node
        errest2s(k) = sqrt((x_gal - uyy)' * A_unit * (x_gal - uyy));
    end

    % indexwise parametric error estimates based on the detail collocation nodes
    perrests = stochcol_est_parametric(X_diff, errest2s, ...
            gridd_diff, list, rule_id);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Q ~= 0 %dimensionaly adaptive
        n = size(X_diff,1);
        costs = ones(1, n);
        for i = 1:n
            row = X_diff(i, :);
            max_dir = find(row > 1, M+Q, 'last');
            costs(i) = costs(i)*max_dir(end);
            for j=1:M+1
                if row(j)>1
                    costs(i) = costs(i)*(2^(row(j)-1) + 1);
                end
            end
        end 
    
        perrest = sum(perrests);
        % parametric error indicator for new dimension
        perrests_new_dim = perrests(1:Q)./costs(1:Q);
        [new_dim_err_max, ind_max] = max(perrests_new_dim);
        % marked index set from the detail index set
        %disp(X_diff)
        %disp(perrests)
            
        if new_dim_err_max > sum(perrests(Q+1:n)./costs(Q+1:n))
            %[~, ind_max] = max(perrests_new_dim);
            Mset = X_diff(ind_max, :);
        else
            [Mset, ~] = dorfler_marking(X_diff(Q+1:n, :), perrests(Q+1:n)./costs(Q+1:n), pmthreshold);
        end
    else
        [Mset, ~] = dorfler_marking(X_diff, perrests, pmthreshold);
        perrest = sum(perrests);
        
    end
 
    % total error indicator
    tot_err_est = serrest + perrest;
                           
    % compute direct estimates of error
    % spatial
    serrest_d = stochcol_direct_estimator(X, ...
        paras_sg, paras_fem, list, pmethod, rhs, aa);
    % parametric
    paras_sg_diff = stochcol_sg(X_diff, rule_id);
    G_diff = stochcol_gmatrices(paras_sg_diff{4}, ...
            paras_sg_diff{5},paras_sg_diff{6}, list);
    ncpts = length(G_diff(:,1));
    ncptsf = size(coords, 1)+ size(coords_diff, 1);
    grid_old_ind = setdiff([1:ncptsf]',grid_diff_ind);
    % construct array containing all the cpoint solution vectors
    sols_all = nan(length(paras_fem{1}), ncptsf);
    sols_all(:,grid_old_ind) = sols; sols_all(:,grid_diff_ind) = sols_diff;
        if ncptsf > ncpts
        fprintf('redundant sparse grid solutions ..\n')
        [gdiff,indx]=setdiff(paras_sg_full{4}, paras_sg_diff{4}, 'rows');
%        gdiff
        iactive=setdiff([1:ncptsf]',indx);
        sols_all=sols_all(:,iactive);
        elseif ncptsf < ncpts, error('Oops ..  fatal logic issue'),
        end
    % compute the energy norm of the SG solution via hierarchical surplus
    perrest_d = sqrt(sum(dot(G_diff, sols_all' * A_unit * sols_all)));
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
    
    err_p_d_iter(iter) = perrest_d;
    err_s_d_iter(iter) = serrest_d;
    error_d_iter(iter) = tot_err_direct;

    sols_iter{iter} = sols;
    paras_sg_iter{iter} = paras_sg;
    paras_fem_iter{iter} = paras_fem;
    dof(iter) = size(sols,1)*size(sols,2);
    [~, qz] = goafem_stochcol_fem_setup(paras_fem{1}, paras_fem{2}, a_unit, G_rhs);
    [~,~,~,~,wp] = stochcol_multilag(paras_sg{4}, paras_sg{1}, paras_sg{2}, paras_sg{3}, rule_id, fun_p);
    Goal_unc = qz.'*sols*wp;
    Goal_unc_iter(iter) = Goal_unc; 
    fprintf('\n GF error %.9f\n',abs(Goal_unc - 0.124120599));
end

endLoopTime = toc(startLoopTime);
%%
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

% final mesh
if pmethod == 1 % P1
    plot_mesh(paras_fem{2}, paras_fem{1}, 'Final FE mesh');
    pause(1)
elseif pmethod == 2
    plot_mesh(paras_fem{6}, paras_fem{5}, 'Final FE mesh');
end
% final nodes
coords = paras_sg{9};
figure(988)
scatter(coords(:,1),coords(:,2),100,'o','filled');
axis square
box on
grid on
title('final collocation nodes in first two directions')

%% postprocess to generate statistics
G = stochcol_gmatrices(paras_sg{4}, paras_sg{5}, paras_sg{6}, list);

[~,~,~,~,multiweights_pdf]...
        = stochcol_multilag(paras_sg{4}, paras_sg{1}, paras_sg{2}, paras_sg{3}, rule_id, fun_p);

MEAN = zeros(size(sols, 1), 1);
Moment_U_squared  = zeros(size(sols, 1), 1);
for k = 1:size(coords,1)
   MEAN = MEAN + multiweights_pdf(k) * sols(:,k); 
end

for k = 1:size(coords,1)
    for l = 1:size(coords, 1)
        Moment_U_squared = Moment_U_squared + sols(:, k).*sols(:, l)*G(k,l); % second order moment
    end
end
VAR = Moment_U_squared - MEAN.^2;

STD = VAR.^(1/2);


final_grid = [paras_sg{4}];
M = size(final_grid, 2);
fprintf('\nFinal dimension of parameter space is %d\n', M-Q)
fprintf('\nFinal sparse grid\n')
disp(final_grid(:, 1:M-Q))

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
    QQ = stochcol_mass_matrix(paras_fem{1}, paras_fem{2});
    QoI = sum(dot(G, sols' * QQ * sols));
    energy = sqrt(sum(dot(G, sols' * A_unit * sols)));
    fprintf('                      QoI is %12.9f\n',QoI/16)
    fprintf('          discrete energy is %12.9f\n',energy)
end
fprintf('\n<strong>Total elapsed time:</strong> %.2f sec\n\n',endLoopTime);
fprintf('To compute a reference solution run the script <strong>referenceSC</strong>\n\n');

% plot final solution
if pmethod == 1% P1
    sc_plotsol_p1(dom_type, MEAN, VAR, paras_fem{2}, paras_fem{1});  
end
