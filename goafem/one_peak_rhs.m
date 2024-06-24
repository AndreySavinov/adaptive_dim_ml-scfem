function f = one_peak_rhs(x1,x2,y)
%RHS_SC   nonzero RHS function for SC test problem
%   f = rhs_sc(x,y);
%   input
%          x1         first coordinate vector
%          x2         second coordinate vector
%          y          y vectore of parameters of f(\bm{x},\bm{y})
%   T-IFISS function: DJS; 14 December 2021
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi
% reference point
y1 = y(1); y2 = y(2);

%whos

% parameters
beta = 50/16;
alpha = (9*y1 + 11)/2;

% RHS (scaled by factor of 1/16)
%d = 32*beta*(1+alpha) - 64*beta^2*( alpha^2*(x-y1).^2 + (y-y2).^2 );
d = 2*beta*(1+alpha) - 4*beta^2*( alpha^2*(x1-y1).^2 + (x2-y2).^2 );
f = d .* exp(-beta*( alpha*(x1-y1).^2 + (x2-y2).^2 ));
end % end function
