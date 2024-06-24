%UNITSQUARE outlines square domain (-4,4)^2
%   TIFISS scriptfile: AS; 17 November 2022
% Copyright (c) 2018 A. Bespalov, L. Rocchi, A. Savinov

  hold on
  plot([-4 4], [-4,-4],'-k');
  plot([4, 4], [-4,4],'-k');
  plot([4,-4], [4,4],'-k');
  plot([-4,-4], [4,-4],'-k');
  grid on
  hold off

% end scriptfile