%STOCH_TESTPROBLEMS setting up stochastic diffusion test problems (Unix version)
%
% NOTE:
% sn: example problem
%
% expansion_type: random-field (coefficient) type
% - 0 a(x,y) = 1, everywhere or cookies coeff
% - 1 Karhunen-Loeve expansion 
% - 2 Synthetic Eigel expansion
% - 3 Bachmayr's Example
%
% dom_type: domain type
% - 1 square domain       (0,1)^2
% - 2 L-shaped domain     (-1,1)^2 \ (-1,0]^2
% - 3 Large square domain (-4,4)^2
%
% right_hand_side.m: define f(x) = 1 or f(x,y) for one-peak-problem
%
% boundary conditions: homogeneous Dirichlet
%
% Latest update: AB; 30 November 2022
% Copyright (c) 2022 A. Bespalov, D. Silvester, F. Xu


fprintf('\n Numerical solution of reference stochastic diffusion problem.');
fprintf('\n Choose specific example:');
fprintf('\n   1. Square domain (0,1)^2, affine random coefficient (Eigel expansion), constant source');
fprintf('\n   2. L-shaped domain, exponential random coefficient (with analytic KL-expansion), constant source');
fprintf('\n   3. Square domain (-4,4)^2, constant coefficient, random one-peak source');
fprintf('\n   4. Square domain (0,1)^2, quadratic random coefficient (with Eigel expansion), constant source');
fprintf('\n   5. Square domain (0,1)^2, exponential random coefficient (with Eigel expansion), constant source');
fprintf('\n   6. L-shaped domain, quadratic random coefficient (with analytic KL-expansion), constant source');
fprintf('\n   7. L-shaped domain, affine random coefficient (Eigel expansion), constant source');
fprintf('\n   8. Square domain (0,1)^2, cookie problem (9 round inclusions), constant source');
fprintf('\n   9. Square domain (0,1)^2, cookie problem (4 square inclusions), discontinuous source');
fprintf('\n  10. Square domain (0,1)^2, affine random coefficient (Bachmayr expansion), constant source\n');

sn = default('',1);

if sn == 1
    % Part I 5.1
    % Square domain (0,1)^2; Eigel affine expansion; unit f(x,y)
    !/bin/cp ./test_problems/standard_trunc_affine.m               ./stochcol_diffusion_grad_and_coeff.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_coeff.m       ./stochcol_diffusion_coeff_spatial_expansion.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_grad_x1.m     ./stochcol_diffusion_grad_x1_spatial_expansion.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_grad_x2.m     ./stochcol_diffusion_grad_x2_spatial_expansion.m
    !/bin/cp ./test_problems/unit_rhs.m                            ./right_hand_side.m
    
    expansion_type  = 2;% Eigel
    dom_type = 1;  % 'domain type: Square (0 , 1)^2'
    
elseif sn == 2
    % Part I 5.2
    % L-shaped domain (-1,1)^2\(-1,0]^2; Exponent of Separable exponents expansion; unit f(x,y)
    !/bin/cp ./test_problems/exponent_trunc_affine.m                               ./stochcol_diffusion_grad_and_coeff.m
    !/bin/cp ./test_problems/separable_exponential_spatial_expansion_coeff.m       ./stochcol_diffusion_coeff_spatial_expansion.m
    !/bin/cp ./test_problems/separable_exponential_spatial_expansion_grad_x1.m     ./stochcol_diffusion_grad_x1_spatial_expansion.m
    !/bin/cp ./test_problems/separable_exponential_spatial_expansion_grad_x2.m     ./stochcol_diffusion_grad_x2_spatial_expansion.m
    !/bin/cp ./test_problems/unit_rhs.m                                            ./right_hand_side.m
    
    expansion_type  = 1;% Separable exponent KL-expansion
    dom_type = 2;  % 'domain type: L-shaped '
    
elseif sn == 3
    % one-peak problem Part II 4.3
    !/bin/cp ./test_problems/standard_trunc_affine.m               ./stochcol_diffusion_grad_and_coeff.m
    !/bin/cp ./test_problems/unit_coeff.m                          ./stochcol_diffusion_coeff_spatial_expansion.m
    !/bin/cp ./test_problems/zero_grad_x1.m                        ./stochcol_diffusion_grad_x1_spatial_expansion.m
    !/bin/cp ./test_problems/zero_grad_x2.m                        ./stochcol_diffusion_grad_x2_spatial_expansion.m
    !/bin/cp ./test_problems/one_peak_rhs.m                        ./right_hand_side.m
    
    expansion_type  = 0;% a(x,y) = 1 everywhere
    dom_type = 3;  % 'domain type: Square (-4, 4)^2
    
elseif sn == 4
    % Square domain (0,1)^2; Squared Eigel affine expansion; unit f(x,y)
    !/bin/cp ./test_problems/square_trunc_affine.m                 ./stochcol_diffusion_grad_and_coeff.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_coeff.m       ./stochcol_diffusion_coeff_spatial_expansion.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_grad_x1.m     ./stochcol_diffusion_grad_x1_spatial_expansion.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_grad_x2.m     ./stochcol_diffusion_grad_x2_spatial_expansion.m
    !/bin/cp ./test_problems/unit_rhs.m                            ./right_hand_side.m
    
    expansion_type  = 2;% Eigel
    dom_type = 1;  % 'domain type: Square (0 , 1)^2'
    
elseif sn == 5
    % Square domain (0,1)^2; Exponent of Eigel affine expansion; unit f(x,y)
    !/bin/cp ./test_problems/exponent_trunc_affine.m               ./stochcol_diffusion_grad_and_coeff.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_coeff.m       ./stochcol_diffusion_coeff_spatial_expansion.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_grad_x1.m     ./stochcol_diffusion_grad_x1_spatial_expansion.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_grad_x2.m     ./stochcol_diffusion_grad_x2_spatial_expansion.m
    !/bin/cp ./test_problems/unit_rhs.m                            ./right_hand_side.m
    
    expansion_type  = 2;% Eigel
    dom_type = 1;  % 'domain type: Square (0 , 1)^2'
    
elseif sn == 6
    % L-shaped domain (-1,1)^2\(-1,0]^2; Square of Separable exponents expansion; unit f(x,y)
    !/bin/cp ./test_problems/square_trunc_affine.m                                 ./stochcol_diffusion_grad_and_coeff.m
    !/bin/cp ./test_problems/separable_exponential_spatial_expansion_coeff.m       ./stochcol_diffusion_coeff_spatial_expansion.m
    !/bin/cp ./test_problems/separable_exponential_spatial_expansion_grad_x1.m     ./stochcol_diffusion_grad_x1_spatial_expansion.m
    !/bin/cp ./test_problems/separable_exponential_spatial_expansion_grad_x2.m     ./stochcol_diffusion_grad_x2_spatial_expansion.m
    !/bin/cp ./test_problems/unit_rhs.m                                            ./right_hand_side.m
    
    expansion_type  = 1;% Separable exponent KL-expansion
    dom_type = 2;  % 'domain type: L-shaped '
    
elseif sn == 7
    % L-shaped domain (-1,1)^2\(-1,0]^2; Eigel affine expansion; unit f(x,y)
    !/bin/cp ./test_problems/standard_trunc_affine.m               ./stochcol_diffusion_grad_and_coeff.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_coeff.m       ./stochcol_diffusion_coeff_spatial_expansion.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_grad_x1.m     ./stochcol_diffusion_grad_x1_spatial_expansion.m
    !/bin/cp ./test_problems/eigel_spatial_expansion_grad_x2.m     ./stochcol_diffusion_grad_x2_spatial_expansion.m
    !/bin/cp ./test_problems/unit_rhs.m                            ./right_hand_side.m
    
    expansion_type  = 2;% Eigel
    dom_type = 2;  % 'domain type: L-shaped '
    
elseif sn == 8
    % Square domain (0,1)^2; Cookies coeffs (9 round subdomains), unit f(x,y)
    !/bin/cp ./test_problems/standard_trunc_affine.m              ./stochcol_diffusion_grad_and_coeff.m
    !/bin/cp ./test_problems/round_cookies9_expansion_coeff.m     ./stochcol_diffusion_coeff_spatial_expansion.m
    !/bin/cp ./test_problems/zero_grad_x1.m                       ./stochcol_diffusion_grad_x1_spatial_expansion.m
    !/bin/cp ./test_problems/zero_grad_x2.m                       ./stochcol_diffusion_grad_x2_spatial_expansion.m
    !/bin/cp ./test_problems/unit_rhs.m                           ./right_hand_side.m
    
    expansion_type  = 0;% Cookies 9 round subdomains
    dom_type = 1;  % 'domain type: Square (0 , 1)^2'
    
elseif sn == 9
    % Square domain (0,1)^2; Square cookies coeffs (4 square subdomains), f(x,y) = 100*F_{\chi_F}
    !/bin/cp ./test_problems/standard_trunc_affine.m              ./stochcol_diffusion_grad_and_coeff.m
    !/bin/cp ./test_problems/square_cookies4_expansion_coeff.m    ./stochcol_diffusion_coeff_spatial_expansion.m
    !/bin/cp ./test_problems/zero_grad_x1.m                       ./stochcol_diffusion_grad_x1_spatial_expansion.m
    !/bin/cp ./test_problems/zero_grad_x2.m                       ./stochcol_diffusion_grad_x2_spatial_expansion.m
    !/bin/cp ./test_problems/central_burst_rhs.m                  ./right_hand_side.m
    
    expansion_type  = 0;% Cookies 4 cookies subdomains
    dom_type = 1;  % 'domain type: Square (0 , 1)^2'
    
elseif sn == 10
    % Square domain (0,1)^2; Bachmayr's coefficients, unit f(x,y)
    !/bin/cp ./test_problems/standard_trunc_affine.m              ./stochcol_diffusion_grad_and_coeff.m
    !/bin/cp ./test_problems/bachmayr_coeff.m                     ./stochcol_diffusion_coeff_spatial_expansion.m
    !/bin/cp ./test_problems/bachmayr_grad_x1.m                   ./stochcol_diffusion_grad_x1_spatial_expansion.
    !/bin/cp ./test_problems/bachmayr_grad_x2.m                   ./stochcol_diffusion_grad_x2_spatial_expansion.m
    !/bin/cp ./test_problems/unit_rhs.m                           ./right_hand_side.m
    
    expansion_type  = 3;% Bachmayr's coefficients
    dom_type = 1;  % 'domain type: Square (0 , 1)^2'
    
else
    error('...Reference problem datafile not found!');
end