function [polys, nodes] = stochcol_onedlagpolys(max_level,rule_id, varargin)
%STOCHCOL_ONEDLAGPOLYS sets up 1D Lagrange polynomials
polys = cell(max_level, 1);
for i = 1:length(polys)
    polys{i} = cell(i, 1);
    if rule_id == 1
        nodes = stochcol_nodes_leja(i);
    elseif rule_id == 2
        nodes = stochcol_nodes_cc(i);
    elseif rule_id == 3
        nodes = stochcol_nodes_leja_gauss(i, varargin{1});
    end
    
    for j = 1:length(nodes{1})
        [polys{i}{j}, ~]= stochcol_1Dlagpoly(nodes{1}(j),nodes{1});
    end
end
