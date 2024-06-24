function a = stochcol_diffusion_coeff_spatial_expansion(x1, x2, yy, input)
    M = input(1);
    a = ones(size(x1)); % a_0
    % Prepare cookies
    i = [1 3 5 1 3 5 1 3 5];
    j = [1 1 1 3 3 3 5 5 5];
    for m = 1:M
        % Set a_m
        if m == 5
            a_m = 0.9*((x1 - i(m)/6).^2 + (x2 - j(m)/6).^2 <= 1/(8*8));
        elseif ismember(m,[2,4,6,8]) 
            a_m = 0.7*((x1 - i(m)/6).^2 + (x2 - j(m)/6).^2 <= 1/(8*8));
        elseif ismember(m,[1,3,7,9])
            a_m = 0.5*((x1 - i(m)/6).^2 + (x2 - j(m)/6).^2 <= 1/(8*8));
        else
            a_m = 0;
        end
        a = a + a_m * yy(m);
    end
end