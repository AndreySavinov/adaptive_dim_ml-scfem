function [serrest_u, serrest_z] = goafem_stochcol_direct_estimatorX( ...
    paras_sg, paras_fem, list, pmethod, rhs_fun, qoi_fun, aa)
% GOAFEM_STOCHCOL_DIRECT_ESTIMATORX computes direct spatial error
% estimates.
%
% input:
%          paras_sg   sparse grid parameters
%         paras_fem   spatial mesh parameter
%              list   list of integrals of Lagrange polynomials
%           pmethod   indicator P1/P2 basis in spatial domain
%           rhs_fun   RHS function cell
%           qoi_fun   QOI function cell
%                aa   diffuison coefficient handle
%
% output:
%       serrest_u     direct spatial error estimation of primal solution 
%       serrest_z     parametric error estimation of dual solution

% SEE ALSO: goafem_singlelevelSC
%
% TR, AS; 28 June 2024


G = stochcol_gmatrices(paras_sg{4}, paras_sg{5}, paras_sg{6}, list);
paras_detail = stochcol_mesh_detail(paras_fem);
MMele_full = (1:size(paras_fem{2}, 1))';      % all elements
MMedge_full = (1:size(paras_detail{7}, 1))';  % all edges
paras_fem_new = stochcol_mesh_refine(MMele_full, MMedge_full, ...
    paras_fem, paras_detail, pmethod);
a_unit = @(x1, x2) ones(size(x1));
%[A_unit,~,~] = goafem_stochcol_fem_setup(paras_fem{1}, paras_fem{2}, ...
%    a_unit, rhs_fun, qoi_fun);
if nargin(rhs_fun{1}) == 3 
   rhs_fun_tmp = rhs_fun;
   rhs_fun_tmp{1} = @(x1, x2) rhs_fun{1}(x1, x2, [0, 0]);
else
   rhs_fun_tmp = rhs_fun;
end
[A_unit_new,~,~] = goafem_stochcol_fem_setup(paras_fem_new{1}, ...
    paras_fem_new{2}, a_unit, rhs_fun_tmp);
coords = paras_sg{9};
xy = paras_fem{1};
xy1 = xy(:,1);
xy2 = xy(:,2);
xy_new = paras_fem_new{1};
xy_new1 = xy_new(:,1);
xy_new2 = xy_new(:,2);
% sols_x = zeros(length(xy), size(coords, 1));
% sols_z = zeros(length(xy), size(coords, 1));
sols_x_new = zeros(length(xy_new), size(coords, 1));
sols_z_new = zeros(length(xy_new), size(coords, 1));
diff_sols_u = zeros(size(sols_x_new));
diff_sols_z = zeros(size(sols_z_new));
parfor k = 1:size(coords, 1)    %parfor
    if nargin(rhs_fun{1}) == 3 
       rhs_fun_tmp = rhs_fun;
       rhs_fun_tmp{1} = @(x1, x2) rhs_fun{1}(x1, x2, coords(k, :));
    else
       rhs_fun_tmp = rhs_fun;
    end
    [u_gal, z_gal, ~] = goafem_stochcol_fem_solver(coords(k, :), ...
        paras_fem, aa, rhs_fun_tmp, qoi_fun);
    
    [u_gal_new, z_gal_new, ~] = goafem_stochcol_fem_solver( ...
        coords(k, :), paras_fem_new, aa, rhs_fun_tmp, qoi_fun, xy, u_gal);
    
    u_gal_interp = griddata(xy1, xy2, u_gal, xy_new1, xy_new2);
    z_gal_interp = griddata(xy1, xy2, z_gal, xy_new1, xy_new2);
    % sols_u(:, k) = u_gal;
    % sols_z(:, k) = z_gal;
    % sols_u_new(:, k) = u_gal_new;
    % sols_z_new(:, k) = z_gal_new;
    u_gal_diff = u_gal_new - u_gal_interp;
    z_gal_diff = z_gal_new - z_gal_interp;
    diff_sols_u(:, k) = u_gal_diff;
    diff_sols_z(:, k) = z_gal_diff;
end
serrest_u = sqrt(sum(dot(G, diff_sols_u' * A_unit_new * diff_sols_u)));
serrest_z = sqrt(sum(dot(G, diff_sols_z' * A_unit_new * diff_sols_z)));
return
