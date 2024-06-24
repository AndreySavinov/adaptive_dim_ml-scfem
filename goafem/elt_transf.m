function [sigptloc,tigptloc] = elt_transf(s,t)
%SUBELT_TRANSF locate 2D Gaussian points on the right reference element

% Top sub-element
x = [  0, 0, 1];    % x-coordinates sub-element
y = [  1, 0, 0];    % y-coordinates sub-element  
  
% Affine mapping from the reference element K^ to the right sub-element
  B = [x(2) - x(1),  y(2) - y(1); ...
       x(3) - x(1),  y(3) - y(1)];
  c = [x(1); y(1)];
  F = c + B'*[s;t];

% Assign new points
  sigptloc = F(1);
  tigptloc = F(2);

end  % end function
