function a = stochcol_diffusion_coeff_spatial_expansion(x1, x2, yy, input)
    M = input(1);
    a = 1.1*ones(size(x1)); % a_0
    for m = 1:M
        % Set a_m
        if m == 1
            a_m = 0.9*(x1<=0.3).*(x1>=0.1).*(x2<=0.3).*(x2>=0.1);
        elseif m == 2
            a_m = 0.6*(x1<=0.9).*(x1>=0.7).*(x2<=0.3).*(x2>=0.1);
        elseif m == 3
            a_m = 0.3*(x1<=0.3).*(x1>=0.1).*(x2<=0.9).*(x2>=0.7);
        elseif m == 4
            a_m = 0.1*(x1<=0.9).*(x1>=0.7).*(x2<=0.9).*(x2>=0.7);
        else
            a_m = 0;
        end
        a = a + a_m * yy(m);
    end
end