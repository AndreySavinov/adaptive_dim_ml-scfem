function paras_sg = stochcol_sg(X, rule_id, varargin)
%STOCHCOL_SG generate attributes of the general sparse grid interpolation
% given an index set and a sequence of 1D collocation nodes
%
% paras_sg = stochcol_sg(X, rule_id)
%
%   Latest update: AS; 24 November 2022
% Copyright (c) 2019 F. Xu

% indset is a subset of X, accounts for nonzero c_i
% sum of tensor products
[cset, indset, mindset] = stochcol_gsgterms(X, rule_id); 
gridd = stochcol_getgrid(mindset);

if numel(varargin) == 0    
    % sum of single-term polynomials
    [gridd, clincombiset, indlincombiset, mindlincombiset, multiweights]...
        = stochcol_multilag(gridd, cset, indset, mindset, rule_id); 
else
    PDF = varargin{1};
    % sum of single-term polynomials, calculation of weights with respect
    % to propability measure
    [gridd, clincombiset, indlincombiset, mindlincombiset, multiweights]...
        = stochcol_multilag(gridd, cset, indset, mindset, rule_id, PDF); 
end

 % coordinates of underlying grid
 coords = stochcol_getcoords(gridd,rule_id);
 % outputs
 paras_sg = {cset, indset, mindset, gridd, clincombiset, ...
             indlincombiset, mindlincombiset, multiweights, coords};
    