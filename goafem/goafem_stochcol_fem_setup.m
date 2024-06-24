function [A, f, g, Q, qe] = goafem_stochcol_fem_setup(xy, evt, coeff_fun, rhs_fun, varargin)
% COMMENTS NEED FINISHING
%
%   [A, f, g, Q, qe] = goafem_stochcol_fem_setup(xy, evt, diff_fun, rhs_fun, qoi_fun)
%
%   input
%                xy  nodal coordinate vector
%               evt  element mapping matrix
%          p_method  approximation method (1 - linear; 2 - quadratic)
%         coeff_fun  function handle for the diffusion coefficient
%           rhs_fun  RHS function cell
%           qoi_fun  QOI function cell
%   output
%          A         stiffness matrix
%          f         rhs vector
%          g         rhs vector
%          Q         mass matrix
%         qe         mass matrix (local)
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function imposebc.
if numel(varargin) == 0
    qoi_fun = rhs_fun;
    qoi_fun{5} = -1;
elseif numel(varargin) == 1
    qoi_fun = varargin{1};
    qoi_fun{5} = -1;
else
    qoi_fun = varargin{1};
    u_gal = varargin{2};
end

[nel, dim] = size(evt); % % no. of elements and no. of nodes per element
if dim == 3
    p_method = 1;
elseif dim == 6
    p_method = 2;
end
%if p_method == 1
%    fprintf('\nsetting up P1 matrices... \n ')
%elseif p_method == 2
%    fprintf('\nsetting up P2 matrices... \n ')
%end
x = xy(:,1);
y = xy(:,2);
n = length(x); % number of nodes

% Gauss point integration rules (3/7/19/73 Gaussian points)
ngpt = 19;
[s,t,wt] = triangular_gausspoints(ngpt);
% inner loop over elements
% (xl_v(i,j),yl_v(i,j)) the coordinate of vertex j of element i
for ivtx = 1:3 % each triangular element has 3 vertices(P1 or P2)
    xl_v(:,ivtx) = x(evt(:,ivtx));
    yl_v(:,ivtx) = y(evt(:,ivtx));
end
% initialise element stiffness matrices
Ae = zeros(nel,dim,dim);
if ismember(qoi_fun{5},[0, 1])
    Be = zeros(nel,dim,dim);
end
qe = zeros(nel,dim,dim);
fe = zeros(nel,dim);
ge = zeros(nel,dim);
% loop over Gauss points
for igpt = 1:ngpt
    sigpt = s(igpt);
    tigpt = t(igpt);
    wtigpt = wt(igpt);
    %  evaluate derivatives etc, of FE basis
    [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v); % P1
    if p_method == 2
        [phi,dphidx,dphidy] = tqderiv(sigpt,tigpt,xl_v,yl_v); % P2
    end
    [coeff] = tgauss(sigpt, tigpt, xl_v, yl_v, coeff_fun);
    [rhs, qoi] = goafem_stochcol_tgauss(sigpt,tigpt,xl_v,yl_v, rhs_fun, qoi_fun);
    for j = 1:dim
        for i = 1:dim
            Ae(:,i,j) = Ae(:,i,j) + wtigpt*coeff(:).*dphidx(:,i).*dphidx(:,j).*invjac(:);
            Ae(:,i,j) = Ae(:,i,j) + wtigpt*coeff(:).*dphidy(:,i).*dphidy(:,j).*invjac(:);
            if qoi_fun{5} == 1
                Be(:,i,j) = Be(:,i,j) + wtigpt.*qoi(:, 1).*(dphidx(:,i)+ dphidy(:,i)) .* phi(:, j);
            elseif qoi_fun{5} == 0
                Be(:,i,j) = Be(:,i,j) + wtigpt.*qoi(:, 1).*phi(:,i) .* phi(:,j) .* jac(:);
            end
            qe(:,i,j) = qe(:,i,j) + wtigpt.*phi(:,i).*phi(:,j).*jac(:);
        end
    end
    for i = 1:dim
        fe(:,i) = fe(:,i) + wtigpt * rhs(:,1) .* phi(:,i).* jac(:); %L2
        fe(:,i) = fe(:,i) - wtigpt * rhs(:,2) .* dphidx(:,i) - wtigpt * rhs(:,3) .* dphidy(:,i); %H1
        ge(:,i) = ge(:,i) + wtigpt * qoi(:,1) .* phi(:,i).*jac(:); %L2
        ge(:,i) = ge(:,i) - wtigpt * qoi(:,2) .* dphidx(:,i) - wtigpt * qoi(:,3) .* dphidy(:,i); %H1
    end
end
% assembly element contributions into global matrix/vector
A = assembly(Ae,evt,n,evt,n);
Q = assembly(qe,evt,n,evt,n);
f = assembly(fe,evt,n);
if qoi_fun{5} == 0
    B = assembly(Be,evt,n,evt,n);
    g = 2*B*u_gal;
elseif qoi_fun{5} == 1
    B = assembly(Be,evt,n,evt,n);
    g = (B+B.')*u_gal;
elseif qoi_fun{5} == 2
    g = assembly(ge,evt,n);
    int_w_u = g.'*u_gal;
    g = 200*int_w_u*g;
elseif qoi_fun{5} == 3
    E = qoi_fun{6};
    g = assembly(ge,evt,n);
    int_w_u = g.'*u_gal;
    g = 200*(int_w_u-E)*g;
else
    g = assembly(ge,evt,n);
end
% A = sparse(n,n);
% f = zeros(n,1);
% g = zeros(n,1);
% for krow=1:dim
%     nrow=evt(:,krow);
%     for kcol=1:dim
%         ncol=evt(:,kcol);
%         A = A + sparse(nrow,ncol,Ae(:,krow,kcol),n,n);
%     end
%     for els=1:nel
%         f(nrow(els),1) = f(nrow(els),1) + fe(els,krow);
%         g(nrow(els),1) = g(nrow(els),1) + ge(els,krow);
%     end
% end
end
