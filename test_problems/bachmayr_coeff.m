function a = stochcol_diffusion_coeff_spatial_expansion(x1, x2, yy, input)
    M = size(yy, 2);
    c = input{2};
    alpha = input{3};
    ell = input{4};
    k1 = input{5};
    k2 = input{6};
    
    a = ones(size(x1));%a_0

    % Set a_m * y_m
    hat_function_1D = @(x) max(0,1-abs(2*x-1));
    for m = 1:M
        a = a + c*(2^(-alpha*ell(m)))*hat_function_1D((2^ell(m))*x1-k1(m)).*hat_function_1D((2^ell(m))*x2-k2(m))*yy(m);
    end
  
end 