%LARGECRACK outlines large crack domain (-1,1)^2 \ (-1,0)x{0} (crack on the left)
%   TIFISS scriptfile: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi
 
  global slope1crack
  global slope2crack
 
  hold on  
  
% Lower diagonal edge
  xP1 = -1; 
  yP1 = -slope1crack;
  plot([0,xP1],[0,yP1],'-k');
 
% Upper diagonal edge
  xP2 = -1; 
  yP2 = -slope2crack;
  plot([0,xP2],[0,yP2],'-k');
     
  plot([xP1,-1], [yP1,-1], '-k');
  plot([-1,1],   [-1,-1],  '-k');
  plot([1,1],    [-1,1],   '-k');
  plot([1,-1],   [1,1],    '-k');
  plot([-1,xP2], [1,yP2],  '-k');
  axis square;

  hold off

% end scriptfile