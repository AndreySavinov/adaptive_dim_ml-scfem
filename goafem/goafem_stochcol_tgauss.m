function [rhs, qoi] = goafem_stochcol_tgauss(s, t, xl_v, yl_v, rhs_fun, qoi_fun)
% GOAFEM_STOCHCOL_TGAUSS calculates right-hand side and QoI at physical coordinates of points in all physical elements
%
%       input
%           s           x coordinate of Gaussian point in refence element
%           t           y coordinate of Gaussian point in refence element
%           xl_v        x coordinate of vertex in physical element
%           yl_v        y coordinate of vertex in physical element
%           rhs_fun     RHS function cell
%           qoi_fun     QOI function cell
%
%       output
%           rhs         Evaluated RHS cell
%           qoi         Evaluated QOI cell
%
% TR, AS; 28 June 2024

nel = length(xl_v(:,1)); % number of elements
zero_v = zeros(nel,1); 
x1 = zero_v; 
x2 = x1;
% vectorized linear shape functions
[phi_e,~,~] = tshape(s, t);
% calculate physical coordinates of points in all physical elements which
% are corresponding to the point (s, t) in the reference element using
% isoparametric relation
for ivtx = 1:3 % size(phi_e, 2) == 3 
    x1 = x1 + phi_e(:, ivtx) .* xl_v(:, ivtx);
    x2 = x2 + phi_e(:, ivtx) .* yl_v(:, ivtx);
end

% Evaluated function handles are split up into the following components:
% L2 | H1(1) | H1(2) | DIV-H1
rhs = [rhs_fun{1}(x1,x2),rhs_fun{2}(x1,x2),rhs_fun{3}(x1,x2),rhs_fun{4}(x1,x2)];
qoi = [qoi_fun{1}(x1,x2),qoi_fun{2}(x1,x2),qoi_fun{3}(x1,x2),qoi_fun{4}(x1,x2)];

end