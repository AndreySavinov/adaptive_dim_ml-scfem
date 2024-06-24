function dadx2 = stochcol_diffusion_coeff_spatial_expansion(x1, x2, yy, input)
    M = size(yy, 2);;
    c = input{2};
    alpha = input{3};
    ell = input{4};
    k1 = input{5};
    k2 = input{6};
    % Set da_m/dx1
    hat_function_1D = @(x) max(0,1-abs(2*x-1));
    dadx2 = zeros(size(x1));
    for m = 1:M
      dadx2 = dadx2 + c*(2^(-(alpha-1)*ell(m)))...
                    *(-2*((2^ell(m))*x2-k2(m)>0.5).*((2^ell(m))*x2-k2(m)<=1) + 2*((2^ell(m))*x2-k2(m)>=0)...
                    .*((2^ell(m))*x2-k2(m)<0.5))...
                    .*hat_function_1D((2^ell(m))*x1-k1(m))*yy(m);
    end
end