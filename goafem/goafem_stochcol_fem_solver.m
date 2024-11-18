function [u_gal, z_gal, Anbc, Qnbc, qenbc, fnbc] = goafem_stochcol_fem_solver(coord, paras_fem, a_fun, rhs, qoi, varargin)
% GOAFEM_STOCHCOL_FEM_SOLVER assembles and solves linear system of Galerkin
% discretization at collocation point coord
%
% input:
%      coord     collocation point 
%  paras_fem     parameters of current mesh
%      a_fun     handle of diffusion coefficient
%        rhs     handle of right-hand side of primal problem
%        qoi     handle of right-hand side of dual problem (QoI)
% output:
%      u_gal         vector of coefficients for primal problem at point coord
%      z_gal         vector of coefficients for dual problem at point coord
%       Anbc         stiffness matrix
%       Qnbc         mass matrix
%      genbc         rhs vector for dual problem
%       fnbc         rhs vector for primal

%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function imposebc.
%
% TR, AS; 28 June 2024


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


end
