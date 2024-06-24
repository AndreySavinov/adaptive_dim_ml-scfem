function [serrest, xy_new_union] = stochcol_direct_estimator_ml(paras_sg, ...
     meshesP, meshesPdetail, sols_ml, list, pmethod, rhs, aa)
% computes direct spatial estimate of the error

G = stochcol_gmatrices(paras_sg{4}, paras_sg{5}, paras_sg{6}, list);
coords = paras_sg{9};

% preallocation of cell data
x_gal_diff = cell(1,size(coords,1));
xy_new = cell(1,size(coords,1));

for k = 1:size(coords, 1) % loop over all ColPoints
    % save parameter values as global variables
    if nargin(rhs) == 3
        rhs_fun = @(x1, x2) rhs(x1, x2, coords(k, :));
    else
        rhs_fun = rhs;
    end
    % uniform refinement of the (current) FE mesh
    % (can be improved using a dedicated T-IFISS function?)
    paras_fem = meshesP(k,:);
    paras_detail = meshesPdetail(k,:);
    MMele_full = (1:size(paras_fem{2}, 1))';      % all elements
    MMedge_full = (1:size(paras_detail{7}, 1))';  % all edges
    paras_fem_new = stochcol_mesh_refine(MMele_full, MMedge_full, ...
                    paras_fem, paras_detail, pmethod);
    % extracting vertex coordinates of the current and refined meshes
    xy1 = paras_fem{1}(:,1);
    xy2 = paras_fem{1}(:,2);
    xy_new{k} = paras_fem_new{1};
    xy_new1 = xy_new{k}(:,1);
    xy_new2 = xy_new{k}(:,2);
    % computing FE approximations on the refined mesh
    [x_gal_new, ~] = stochcol_fem_solver(coords(k, :), paras_fem_new, ...
        aa, rhs_fun);
    % interpolating the current-mesh FE solution onto the refined mesh
    x_gal_interp = griddata(xy1, xy2, sols_ml{k}, xy_new1, xy_new2);
    % computing the difference of two FE solutions
    x_gal_diff{k} = x_gal_new - x_gal_interp;
end

% finding the union (i.e., coarsest common refinement) of refined FE meshes
% over all ColPoints
[xy_new_union, evt_new_union] = stochcol_meshes_union(xy_new, size(coords, 1));
% extracting vertex coordinates of the union mesh
xy_new_union1 = xy_new_union(:,1);
xy_new_union2 = xy_new_union(:,2);
% setting the unit coefficient for computing spatial H1 seminorm
a_unit = @(x1, x2) ones(size(x1));
[A_unit_union,~] = stochcol_fem_setup(xy_new_union, evt_new_union, a_unit, a_unit);

% preallocation of cell data
diff_sols = nan(size(xy_new_union,1), size(coords, 1));

for k = 1:size(coords, 1)
    % interpolating the difference of two FE solutions onto the union mesh
    xy_new1 = xy_new{k}(:,1);
    xy_new2 = xy_new{k}(:,2);
    diff_sols(:,k) = griddata(xy_new1, xy_new2, x_gal_diff{k}, xy_new_union1, xy_new_union2);
end

% computing the Bochner norm via hierarchical surplus
serrest = sqrt(sum(dot(G, diff_sols' * A_unit_union * diff_sols)));

return
