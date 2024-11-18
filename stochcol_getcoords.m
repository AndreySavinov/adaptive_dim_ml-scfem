function coords = stochcol_getcoords(grid, rule_id)
%STOCHCOL_GETCOORDS get coordinates of collocation nodes in parametric space
%from their grids
%
%  coords = stochcol_getcoords(grid,rule_id)
%
%   TIFISS function: AS 28 June 2024
% Copyright (c) 2024 A. Savinov, F. Xu


% coords is one-to-one with grid
if rule_id == 1
    levels = max(max(grid));
    nodes = stochcol_nodes_leja(levels);
elseif rule_id == 2
    if max(max(grid)) == 1
        levels = max(max(grid));
    else
        levels = log2(max(max(grid)-1))+1;
    end
    nodes = stochcol_nodes_cc(levels);
elseif rule_id == 3
    if max(max(grid)) == 1
        levels = max(max(grid));
    else
        levels = (max(max(grid))+1)/2;
    end
    nodes = stochcol_nodes_leja_gauss(levels);
end
[rg, cg]=size(grid);
coords = zeros(rg,cg);
for i = 1:rg
    for j = 1:cg
        coords(i,j) = nodes{1}(grid(i,j));
    end
end
end

