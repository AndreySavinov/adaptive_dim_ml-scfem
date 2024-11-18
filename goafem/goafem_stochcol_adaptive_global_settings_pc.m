% Script specifiyng the main parameters of goafem algotrihm for windows
% machine
% AS, TR; 28 June 2024
clear
close ALL
clc % AB
% set parameters
fprintf('\n Numerical solution of goal-oriented stochastic diffusion problem.');
fprintf('\n Choose specific example or specify your choice:');
fprintf('\n     1.  L-shape domain: [-1,1]^2 - [-1,0)^2, standard truncated affine Eigel expansion, unit RHS, integral over a square subdomain')
fprintf('\n     2.  Left-Crack domain [-1,1]^2 - [-1,0)x{0}, exponent of truncated Eigel expansion, unit RHS, estimation of pointwise values');
fprintf('\n     3.  Unit square domain: [0,1]^2, standard truncated affine Eigel expansion, Mommer-Stevenson RHS, second moment of a linear goal functional')
fprintf('\n     4.  L-shape domain: [-1,1]^2 - [-1,0)^2, exponent of truncated Eigel expansion, unit RHS, non-linear convection')
fprintf('\n     5.  One-Peak Problem: square domain: [-4,4]^2, unit coefficient, random one-peak source, intergral of u^2');
fprintf('\n     0.  Your Choice\n')
sn = default('default = 1', 1);

if ismember(sn,[1, 2, 3, 4, 5])
   QOI_type = sn;
   if ismember(sn,[1, 3])
       expansion_type = 2;
       rf_type = 3;
   elseif ismember(sn,[2, 4])
       expansion_type = 2;
       rf_type = 1;
   else
       rf_type = 3;
       expansion_type = 3;
       QOI_type = 7;
       M = 2;
       Q = 0;
   end
   if ismember(sn,[1, 4])
       dom_type = 2;
       RHS_type = 1;
   elseif sn == 2
       dom_type = 3;
       RHS_type = 1;
   elseif sn == 3
       dom_type = 1;
       RHS_type = 2;
   else
       dom_type = 5;
       RHS_type = 4;
   end
else
    % domain type
    fprintf('\n Choose the spatial domain of the problem:\n')
    fprintf('1. Unit square domain: [0,1]^2 \n')
    fprintf('2. L-shape domain: [-1,1]^2 - [-1,0)^2 \n')
    fprintf('3. Left-Crack domain [-1,1]^2 - [-1,0)x{0} \n')
    fprintf('4. Large square domain: [-1,1]^2 \n')
    fprintf('5. Enormous square domain: [-4,4]^2 \n')
    dom_type = default('(default is 2)', 2);
    
    
    fprintf('choose type of diffusion coefficient ');
    fprintf('\n     1.  exponential of truncated affine expansion')
    fprintf('\n     2.  square of truncated affine expansion');
    fprintf('\n     3.  standard truncated affine (K-L) expansion\n')
    rf_type = default('default = 3',3);
    if ~ismember(rf_type,[1, 2, 3])
        error('Wrong selection for type of diffusion coefficient!')
    end
    
    fprintf('\nchoose type of spatial expansion coefficient ');
    fprintf('\n     1.  separable exponential ')
    fprintf('\n     2.  synthetic (Eigel slow) expansion');
    fprintf('\n     3.  unit\n');
    expansion_type  = default('default = 2',2);
    if ~ismember(expansion_type,[1,2,3])
        error('Wrong selection for expansion type!')
    end
end    
    
if rf_type == 1
    system('copy .\test_problems\exponent_trunc_affine.m               .\stochcol_diffusion_grad_and_coeff.m');
elseif rf_type == 2
    system('copy .\test_problems\square_trunc_affine.m                 .\stochcol_diffusion_grad_and_coeff.m');
elseif rf_type == 3
   system('copy .\test_problems\standard_trunc_affine.m               .\stochcol_diffusion_grad_and_coeff.m');
end


if expansion_type == 1 % seperable exponential
    sigma = default('SE standard deviation (default is 0.5)', 0.5);
    ax = 1;
    ay = 1;
    correl_x = default('correlation length in x (default is 1)', 1);
    correl_y = default('correlation length in y (default is 1)', 1);
    system('copy .\test_problems\separable_exponential_spatial_expansion_coeff.m       .\stochcol_diffusion_coeff_spatial_expansion.m');
    system('copy .\test_problems\separable_exponential_spatial_expansion_grad_x1.m     .\stochcol_diffusion_grad_x1_spatial_expansion.m');
    system('copy .\test_problems\separable_exponential_spatial_expansion_grad_x2.m     .\stochcol_diffusion_grad_x2_spatial_expansion.m');
elseif expansion_type == 2 % Eigel
    sigma = default('Eigel standard deviation (default is 0.547)', 0.547);
    system('copy .\test_problems\eigel_spatial_expansion_coeff.m       .\stochcol_diffusion_coeff_spatial_expansion.m');
    system('copy .\test_problems\eigel_spatial_expansion_grad_x1.m     .\stochcol_diffusion_grad_x1_spatial_expansion.m');
    system('copy .\test_problems\eigel_spatial_expansion_grad_x2.m     .\stochcol_diffusion_grad_x2_spatial_expansion.m');
elseif expansion_type == 3 % unit
    rf_type = 3;
    system('copy .\test_problems\unit_coeff.m                          .\stochcol_diffusion_coeff_spatial_expansion.m');
    system('copy .\test_problems\zero_grad_x1.m                        .\stochcol_diffusion_grad_x1_spatial_expansion.m');
    system('copy .\test_problems\zero_grad_x2.m                        .\stochcol_diffusion_grad_x2_spatial_expansion.m');
end

fprintf('\nchoose type of random variable ');
fprintf('\n     1.  Uniform ')
fprintf('\n     2.  Truncated Gaussian');
fprintf('\n     3.  Gaussian\n');
rv_id = default('default = 1', 1);
if ~ismember(rv_id,[1,2,3])
    error('Wrong selection for type of random variable!')
end
if rv_id == 2
    sigma = 0.5;
end

fprintf('\n Proceeding with P1 finite element approximations. \n')
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
    fprintf('\n Using linear bubble functions for error estimation. \n')
    pestim = 1;
%     pestim = default('\nError estimation: linear/quadratic bubble functions 1/2? (default 1)',1);
    paras_fem_errest = [paras_fem_errest, pestim];
    if pestim == 1
        fprintf('\n 2-level error estimators. \n')
        estimtype = 3;        
        % Error estimation type (for P1 approximations only)
        % 1 - eY hierarchical estimator (elementwise residual problems)
        % 2 - eY hierarchical estimator (assembled system for the residual problems)
        % 3 - 2-level error estimator
%         fprintf('Estimator type:\n');
%         fprintf('   1. hierarchical estimator (elementwise residual problems)\n');
%         fprintf('   2. hierarchical estimator (fully assembled system for residual problem)\n');
%         fprintf('   3. 2-level estimator\n');
%         estimtype = default('(default 1)',1);
        paras_fem_errest = [paras_fem_errest, estimtype];
        if ~ismember(estimtype,[1,2,3]), error('Estimator type not allowed!'); end
        %
        % Marking elements or edges 1/2?
        % This depends on the estimator type:
        % - ypestim=1   -> only elements can be marked, i.e., markedgelem=1;
        % - ypestim=2/3 -> both elements and edges can be marked, i.e.,markedgelem=1/2
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


if dom_type > 2
    if dom_type == 5
        dom_type = 3;
    else
        dom_type = dom_type + 1;
    end
end

dom_paras = dom_type;
if dom_type == 2
    dom_paras = dom_type;
    mesh_type = default('\nStructured/unstructured mesh 1/2 (default 1)',1);
    dom_paras = [dom_paras, mesh_type];
end

% parameter space [-L, L]^M
L = 1;
if sn~=5
    M = default('dimension of parametric space (default is 4)',4);
    Q = default('Number of considered additional dimensions of parametric space (default is 0)', 0);
end


if rv_id ~= 3
    fprintf('\nchoose type of collocation nodes');
    fprintf('\n     1.  Leja ')
    fprintf('\n     2.  CC\n');
    rule_id = default('default is CC nodes',2);
    if ~ismember(rule_id,[1,2]) 
        error('Wrong type of collocation nodes!') 
    end 
else
    fprintf('\nLeja collocation nodes for Gaussian weight.\n');
    rule_id = 3;
end
if ~ismember(rule_id,[1,2,3]) % AB
    error('Wrong type of collocation nodes!') % AB
end % AB

% KL expansion
if expansion_type == 1 % seperable exponential
    input = [M, ax, ay, correl_x, correl_y, sigma];
elseif expansion_type == 2 % Eigel
    input = [M, sigma, 2];
elseif expansion_type == 3 % Eigel
    input = M;
end
%KL_DATA = stoch_kl(input,expansion_type);
% function handle for diffusion coefficient w.r.t. spatial and parametric
% variabels
%[aa, a] = stochcol_diffusion_coeff_fun(KL_DATA, rf_type, M);
% function handle for gradients of diffusion coefficient
%[aax1, aax2] = stochcol_diffusion_gradcoeff_fun(KL_DATA, rf_type, M, aa, a);

%%%%%%%%%%%%%%%%%%%%%%
% function handle for diffusion coefficient w.r.t. spatial(x) and
% parametric(y) variabels [h(x,y) = h_0 + \sum_{i=1}^{M} h_i(x) * w(y)]
a = @(x1, x2, yy) stochcol_diffusion_coeff_spatial_expansion(x1, x2, yy, input); 

% function handle for gradients of diffusion coefficient w.r.t. spatial and parametric
% variabels (\frac{\partial h(x, y)}/{\partial x_1}, \frac{\partial h(x, y)}/{\partial x_2})
ax1 = @(x1, x2, yy) stochcol_diffusion_grad_x1_spatial_expansion(x1, x2, yy, input);
ax2 = @(x1, x2, yy) stochcol_diffusion_grad_x2_spatial_expansion(x1, x2, yy, input);

% a = g(h(x,y)), \frac{\partial a(x, y)}/{\partial x_1}, \frac{\partial a(x, y)}/{\partial x_2}
[aa, aax1, aax2] = stochcol_diffusion_grad_and_coeff(a, ax1, ax2);
%%%%%%%%%%%%%%%%%%%%%%%

goafem_stochcol_qtychoice;

% probability density function
if rv_id == 1 % uniform
    fun_p = @(x) 0.5;
    rv_name = ['uniform'];
elseif rv_id == 2 % truncated Gaussian
    sigma_g = 1;
    fun_p = @(x) ...
        exp(-x.^2./sigma_g./sigma_g./2)...
        /sigma_g/sqrt(2*pi)/erf(L/sqrt(2)/sigma_g);
    rv_name = ['gaussian_','sig_', num2str(sigma_g)];
elseif rv_id == 3 %Gaussian
    fun_p = @(x) ...
        exp(-x.^2./2)/sqrt(2*pi)/erf(10/sqrt(2));
    rv_name = ['gaussian_full','sig_', num2str(1)];
end


% 1D Lagrange polynomials
if rule_id == 1 || rule_id == 3
    max_level = 9;
elseif rule_id == 2
    max_level = 7;
end

try 
    if rule_id == 1 && rv_id == 1 && L == 1 && max_level <= 9
        load('precomputation_leja9_uniform.mat')
    elseif rule_id == 1 && rv_id == 2 && L == 1 && sigma_g == 1 && max_level <= 9
        load('precomputation_leja9_gaussian_sig_1.mat')
    elseif rule_id == 2 && rv_id == 1 && L == 1 && max_level <= 7
        load('precomputation_cc7_uniform.mat')
    elseif rule_id == 2 && rv_id == 2 && L == 1 && sigma_g == 1 && max_level <= 7
        load('precomputation_cc7_gaussian_sig_1.mat')
    elseif rv_id == 3 && max_level <= 9 
        load('precomputation_leja9_gaussian_full_sig_1.mat')
    else
        error('Precomputation data does not exist.')
    end
catch 
    fprintf('\nWarning: Precomputation data does not exist. \n');
    pause(1);
    fprintf('This data can be computed, but might take a long time depending on parameter choice (possibly ~6 hours). \n');
    pause(1);
    run = default('Do you still want to continue (default YES)? YES/NO (1/0)', 1);
    if run == 1
        polys = stochcol_onedlagpolys(max_level, rule_id);
        [list, listy, listy2, listy3, listy4] = goafem_stochcol_uni_int(fun_p, polys, L);
        if rule_id == 1
            filename = ['precomputation_','leja', num2str(max_level),'_',rv_name];
        elseif rule_id == 2
            filename = ['precomputation_','cc', num2str(max_level),'_',rv_name];
        else
            error('\nPrecomputation data not generated.')
        end
    	save(filename, 'max_level', 'rule_id', 'rv_id', 'fun_p', 'polys', 'L', 'list', 'listy', 'listy2', 'listy3', 'listy4')
    else
    	error('\nPrecomputation data does not exist.');
    end
end

% Marking threshold parameters for indices
pmthreshold = default('\nThreshold parameter for marking indices (default 0.3)',0.3);