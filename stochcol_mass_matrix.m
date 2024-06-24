function [Q, P] = stochcol_mass_matrix(xy,mv,varargin)
%STOCHCOL_MASS_MATRIX finite element matrices generator
%
%   [Q] = stochcol_mass_matrix(xy, mv)
%
%   input
%                xy  nodal coordinate vector
%                mv  element mapping matrix
%   output
%          Q         mass matrix
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function imposebc.
%    TIFISS function: AB; 14 December 2021.
% Copyright (c) 2017 A. Bespalov, D.J. Silvester

[nel, dim] = size(mv); % % no. of elements and no. of nodes per element
if dim == 3
    p_method = 1;
elseif dim == 6
    p_method = 2;
end

x = xy(:,1);
y = xy(:,2);
n = length(x); % number of nodes
% Gauss point integration rules (3/7/19/73 Gaussian points)
ngpt = 7;
[s,t,wt] = triangular_gausspoints(ngpt);
% inner loop over elements
for ivtx = 1:3 % each triangular element has 3 vertices(P1 or P2)
    xl_v(:,ivtx) = x(mv(:,ivtx));
    yl_v(:,ivtx) = y(mv(:,ivtx));
end
% initialise element mass matrices
Qe = zeros(nel,dim,dim);
Pe = zeros(nel,dim,dim);
% loop over Gauss points
for igpt = 1:ngpt
    sigpt = s(igpt);
    tigpt = t(igpt);
    wtigpt = wt(igpt);
    if ~isempty(varargin)
        [coeff] = tgauss(sigpt, tigpt, xl_v, yl_v, varargin{1});
        wtigpt = wtigpt*coeff(:);
    end
    %  evaluate derivatives etc, of FE basis
    [jac,~,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v); % P1
    if p_method == 2
        [phi,~,~] = tqderiv(sigpt,tigpt,xl_v,yl_v); % P2
    end
    for j = 1:dim
        for i = 1:dim
            Qe(:,i,j) = Qe(:,i,j) + ...
                wtigpt .* phi(:,i) .* phi(:,j) .* jac(:);
            Pe(:,i,j) = Pe(:,i,j) + ...
                wtigpt .* phi(:,i) .* (dphidx(:,j) +  dphidy(:,j));
        end
    end
end
% assembly element contributions into global matrix/vector
Q = assembly(Qe,mv,n,mv,n);
P = assembly(Pe,mv,n,mv,n);

end
