%STOCHCOL_ADAPTIVE_GLOBAL_SETTINGSP1 sets parameters for test problems and approximations
%
% Latest update: AB; 12 January 2023
% Copyright (c) 2022 A. Bespalov, D. Silvester, F. Xu

% parameter space [-L, L]^M
L = 1; % If L is changed, modifications will be needed to e.g. stochcol_1Dlagpoly.m
if sn ~= 3 
    if sn == 8
        fprintf('\nFor this test problem, the dimension of parametric space is fixed as 9');
        M = 9;
    elseif sn == 9
        fprintf('\nFor this test problem, the dimension of parametric space is fixed as 4');
        M = 4;
    else
        M = default('Dimension of parametric space (default is 1)', 1);
        Q = default('Number of considered additional dimensions of parametric space (default is 1)', 1);
    end

    % function handle for the right-hand side w.r.t. spatial variables
    rhs = @(x1, x2) right_hand_side(x1, x2);

    % selection of distribution
    fprintf('\nchoose type of random variable ');
    fprintf('\n     1.  Uniform ')
    fprintf('\n     2.  Truncated Gaussian');
    fprintf('\n     3.  Gaussian\n');
    rv_id = default('Default is uniform', 1);
else % sn = 3
    fprintf('\nFor this test problem, the dimension of parametric space is fixed as 2');
    M = 2;
    
    % function handle for the right-hand side w.r.t. spatial and parametric variables
    rhs = @(x1, x2, y) right_hand_side(x1, x2, y);

    fprintf('\nDistribution of random variable is uniform\n');
    rv_id = 1;
end


if ~ismember(rv_id,[1,2,3])
    error('Wrong selection for type of random variable!')
elseif rv_id == 1
    fun_p = @(x) 1/(2*L);
elseif rv_id == 2
    sigma_g = default('Sigma for 1-D Gaussian distribution (default is 1)', 1);
    fun_p = @(x) ...
        exp(-x.^2./sigma_g./sigma_g./2)...
        /sigma_g/sqrt(2*pi)/erf(L/sqrt(2)/sigma_g);
elseif rv_id == 3
    fun_p = @(x) ...
        exp(-x.^2./2)/sqrt(2*pi)/erf(10/sqrt(2));
end

% KL expansion (coefficient model)
if expansion_type == 1 % exponent field of seperable exponents
    sigma = default('SE standard deviation (default is 0.5)', 0.5);
    ax = 1;
    ay = 1;
    correl_x = default('correlation length in x (default is 1)', 1);
    correl_y = default('correlation length in y (default is 1)', 1);
    simulationParams = [M, ax, ay, correl_x, correl_y, sigma];

elseif expansion_type == 2 % Eigel
    sigma = default('Eigel standard deviation (default is 0.547)', 0.547);
    simulationParams = [M, sigma, 2];

elseif expansion_type == 3 % Bachmair
    c = default('Coefficient c (default is 0.1)',0.1);
    alpha = default('Coefficient alpha (default is 1)',1);
    ell = 0;
    k1 = 0;
    k2 = 0;
    i = 0;
    while length(ell)<M
       i=i+1;
       ell = [ell,i*ones(1,(2^(i+1)-1)^2)];
       centers_hat = 0:0.5:(2^i-1);
       [k2new,k1new] = meshgrid(centers_hat,centers_hat);
       k1 = [k1,(k1new(:))'];
       k2 = [k2,(k2new(:))'];
        % remove triples with k1 and k2 both non-integers
       good_triples = logical((floor(k1)==k1)+(floor(k2)==k2));
       ell = ell(good_triples);
       k1 = k1(good_triples);
       k2 = k2(good_triples);
    end
    simulationParams = {M, c, alpha, ell, k1, k2};

elseif expansion_type == 0 % unit coefficient a(x,y) = 1 or cookies problems coeffs
    simulationParams = M;
end

fprintf('\nPiecewise linear (P1) finite element approximation\n');
pmethod = 1;

% fprintf('\nchoose type of finite element approximation');
% fprintf('\n     1.  P1 ')
% fprintf('\n     2.  P2\n');
% pmethod = default('default is P1', 1);
% if ~ismember(pmethod,[1,2])
%     error('Wrong selection for approximation method!')
% end

% Red/Bisec3 for spatial error estimation 1/2? See:
% Ainsworth, Oden, A posteriori error estimation in finite element analysis,
% Wiley, 2000 - Figure 5.2 (p. 87) for the basis functions in both cases.
subdivPar = 2;

% Error estimation and estimators type
paras_fem_errest = subdivPar;
if pmethod == 1
    pestim = default('\nError estimation: linear/quadratic bubble functions 1/2? (default 1)',1);
    paras_fem_errest = [paras_fem_errest, pestim];
    if pestim == 1
        % Error estimation type (for P1 approximations only)
        % 1 - eY hierarchical estimator (elementwise residual problems)
        % 2 - eY hierarchical estimator (assembled system for the residual problems)
        % 3 - 2-level error estimator
        fprintf('Estimator type:\n');
        fprintf('   1. hierarchical estimator (elementwise residual problems)\n');
        fprintf('   2. hierarchical estimator (fully assembled system for residual problem)\n');
        fprintf('   3. 2-level estimator\n');
        estimtype = default('(default 1)',1);
        paras_fem_errest = [paras_fem_errest, estimtype];
        if ~ismember(estimtype,[1,2,3]), error('Estimator type not allowed!'); end
        %
        % Marking elements or edges 1/2?
        % This depends on the estimator type:
        % - ypestim=1   -> only elements can be marked, i.e., markedgelem=1;
        % - ypestim=2/3 -> both elements and edges can be marked, i.e.,markedgelem=1/2.
        if estimtype == 1
            % Marking elements only
            markedgelem = 1;
        elseif estimtype == 2
            markedgelem = default('Marking elements/edges 1/2 (default 1)',1);
        else%estimtype = 3
            markedgelem = default('Marking elements/edges 1/2 (default 2)',2);
        end
        paras_fem_errest = [paras_fem_errest, markedgelem];
        if ~ismember(markedgelem,[1,2])
            error('Marking type not allowed!');
        end
    elseif pestim == 2
        % Marking elements
        fprintf('Using hierarchical estimator (elementwise residual problems)\n');
        markedgelem = 1;
        paras_fem_errest = [paras_fem_errest, markedgelem];
    else
        error('Estimation type not allowed!');
    end
elseif pmethod == 2
    % Marking elements
    fprintf('Using hierarchical estimator (elementwise residual problems)\n');
    markedgelem = 1;
    paras_fem_errest = [paras_fem_errest, markedgelem];
end

% Marking threshold parameters for both elements/edges and indices:
% 1 - maximum strategy:
%     large threshold -> small set of marked elements/edges
%     small threshold -> large set of marked elements/edges
% 2 - Doerfler (equilibration) strategy:
%     large threshold -> large set of marked elements/edges
%     small threshold -> small set of marked elements/edges
markstrat   = default('Marking strategy: maximum or equilibration 1/2? (default 2)',2);
smthreshold = default('Threshold parameter (default 0.3)',0.3);
paras_fem_errest = [paras_fem_errest, markstrat, smthreshold];

% domain type
if dom_type == 1
    dom_paras = dom_type;
elseif dom_type == 2
    dom_paras = dom_type;
    mesh_type = 1;%default('\nStructured/unstructured mesh 1/2 (default 1)',1);
    dom_paras = [dom_paras, mesh_type];
elseif dom_type == 3
    dom_paras = dom_type;
end

if rv_id ~= 3
    fprintf('\nchoose type of collocation nodes');
    fprintf('\n     1.  Leja ')
    fprintf('\n     2.  CC\n');
    rule_id = default('default is Leja nodes',2);
    if ~ismember(rule_id,[1,2]) 
        error('Wrong type of collocation nodes!') 
    end 
else
    fprintf('\nLeja collocation nodes for Gaussian weight.\n');
    rule_id = 3;
end

% 1D Lagrange polynomials
if rule_id == 1 || rule_id == 3
    max_level = 9;
    max_level_p = max_level;
elseif rule_id == 2
    max_level = 7;
    max_level_p = max_level;
end

try 
    if rule_id == 1 && rv_id == 1 && L == 1 && max_level <= 9
        load('precomputation_leja_9_uniform1.mat')
    elseif rule_id == 1 && rv_id == 2 && L == 1 && sigma_g == 1 && max_level <= 9
        load('precomputation_leja_9_gaussian1_1.mat')
    elseif rv_id == 3 && max_level <= 9 
        load('precomputation_leja_9_gaussian_full10_1.mat')
    elseif rule_id == 2 && rv_id == 1 && L == 1 && max_level <= 7
        load('precomputation_cc_7_uniform1.mat')
    elseif rule_id == 2 && rv_id == 2 && L == 1 && sigma_g == 1 && max_level <= 7
        load('precomputation_cc_7_gaussian1_1.mat')
    else
        error('m')
    end
catch 
    decider = default('Would you like to wait for additional computation about an hour (default NO)? YES/NO (1/2) ', 2);
    if decider == 1
        if rv_id == 1
            [polys, list] = precompute_collocation_node(rule_id, 1, L, max_level);
        elseif rv_id == 2
            [polys, list] = precompute_collocation_node(rule_id, 2, L, max_level, sigma_g);
        elseif rv_id == 3
            [polys, list] = precompute_collocation_node_new(3, 3, 10, max_level);
        end
    else
        error('Please, rerun SinglelevelSC with a set of parameters for which points are precomputed');
    end
end

% function handle for diffusion coefficient w.r.t. spatial(x) and
% parametric(y) variabels [h(x,y) = h_0 + \sum_{i=1}^{M} h_i(x) * w(y)]
a = @(x1, x2, yy) stochcol_diffusion_coeff_spatial_expansion(x1, x2, yy, simulationParams); 

% function handle for gradients of diffusion coefficient w.r.t. spatial and parametric
% variabels (\frac{\partial h(x, y)}/{\partial x_1}, \frac{\partial h(x, y)}/{\partial x_2})
ax1 = @(x1, x2, yy) stochcol_diffusion_grad_x1_spatial_expansion(x1, x2, yy, simulationParams);
ax2 = @(x1, x2, yy) stochcol_diffusion_grad_x2_spatial_expansion(x1, x2, yy, simulationParams);

% a = g(h(x,y)), \frac{\partial a(x, y)}/{\partial x_1}, \frac{\partial a(x, y)}/{\partial x_2}
[aa, aax1, aax2] = stochcol_diffusion_grad_and_coeff(a, ax1, ax2);

% Marking threshold parameters for indices
pmthreshold = default('Threshold parameter for marking indices (default 0.3)', 0.3);
