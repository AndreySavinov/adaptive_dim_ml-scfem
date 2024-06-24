function [u_gal, z_gal, Anbc, Qnbc, qenbc, fnbc] = goafem_stochcol_fem_solver(coord, paras_fem, a_fun, rhs, qoi, varargin)
% COMMENTS NEED FINISHING
%
%  [u_gal, z_gal, Anbc] = goafem_stochcol_fem_solver(coord, paras_fem, a_fun, rhs, qoi)
%
%

xy     = paras_fem{1};
evt    = paras_fem{2};
bound  = paras_fem{3};
nvtx = length(xy);

% solvingTime = tic;
ayy = @(x1, x2) a_fun(x1, x2, coord);
% check the minimum of a
a_val = ayy(xy(:, 1), xy(:, 2));
% fprintf('\nmininum of diffusion coefficient %9.5f\n', min(a_val))
if min(a_val) < 0
    error('Negative diffusion coefficient!')
end

% FEM matrix/vector assembly
[Anbc,fnbc,~, ~, ~] = goafem_stochcol_fem_setup(xy, evt, ayy, rhs);
% impose boundary condition and solve
nodes = (1:nvtx)';
intern = nodes(~ismember(nodes,bound));
xbd = xy(bound,1); % x coordinate of boundary nodes
ybd = xy(bound,2); % y coordinate of boundary nodes
bc = specific_bc(xbd,ybd);

u_gal = zeros(nvtx,1);
u_gal(bound) = bc;
[A,f] = imposebc(Anbc,fnbc,intern,bound,bc);
u_gal(intern) = A\f;

if numel(varargin) > 0
    xy_old = varargin{1};
    u_gal_old = varargin{2};
    u_gal_interp = griddata(xy_old(:,1), xy_old(:, 2), u_gal_old, xy(:, 1), xy(:, 2));
    [Anbc,~,gnbc,Qnbc,qenbc] = goafem_stochcol_fem_setup(xy, evt, ayy, rhs, qoi, u_gal_interp);
else
    [Anbc,~,gnbc,Qnbc,qenbc] = goafem_stochcol_fem_setup(xy, evt, ayy, rhs, qoi, u_gal);
end
[B,g] = imposebc(Anbc,gnbc,intern,bound,bc);
z_gal = zeros(nvtx,1);
z_gal(bound) = bc;
z_gal(intern) = B\g;

%u_gal(intern) = minres(A,b,1e-8,9999);
%fprintf('solution maximum %9.5f\n',max(u_gal))
%energy = sqrt(u_gal' * Anbc * u_gal);
%fprintf(' solution energy %9.5f\n',energy);
% fprintf('Solving took %.5f sec\n',toc(solvingTime));

end
