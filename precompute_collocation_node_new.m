function [polys, list] = precompute_collocation_node_new(rule_id, rv_id, L, varargin)

% 1D Lagrange polynomials
if numel(varargin) == 0
    if rule_id == 1 || rule_id == 3
        max_level = 9;
    elseif rule_id == 2
        max_level = 7;
    end
else
    max_level = varargin{1};
end

if ~ismember(rv_id,[1,2,3])
    error('Wrong selection for type of random variable!')
elseif rv_id == 1
    pdf = @(x) 1/(2*L);
    str_L = num2str(round(L, 2));
    str_L(str_L == '.') = '_';
    surname = ['uniform', str_L];
elseif rv_id == 2
    sigma = varargin{2};
    pdf = @(x) ...
        exp(-x.^2./sigma./sigma./2)...
        /sigma/sqrt(2*pi)/erf(L/sqrt(2)/sigma);

    str_L = num2str(round(L, 2));
    str_L(str_L == '.') = '_';
    
    str_sigma = num2str(round(sigma, 2));
    str_sigma(str_sigma == '.') = '_';
    surname = ['gaussian', str_L, '_', str_sigma];
elseif rv_id == 3
    pdf = @(x) ...
        exp(-x.^2./2)/sqrt(2*pi)/erf(L/sqrt(2));
    
    surname = ['gaussian_full', num2str(L), '_', num2str(1)];
end

[polys, nodes] = stochcol_onedlagpolys(max_level, rule_id, 1);
list = stochcol_uni_int(pdf, polys, L);

if rule_id == 1 || rule_id == 3
    filename = ['precomputation_','leja_', num2str(max_level),'_'];
elseif rule_id == 2
    filename = ['precomputation_','cc_', num2str(max_level),'_'];
else
    error('Not supported type of collocation nodes')
end

filename = [filename, surname];

save(filename, 'polys', 'list', 'nodes');