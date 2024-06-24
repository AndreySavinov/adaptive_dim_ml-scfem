function KL_coeff = stochcol_diffusion_coeff_spatial_expansion(x, y, yy, input)
    norv2d = size(yy, 2);
    ax = input(2);
    ay = input(3);
    correl_x = input(4);
    correl_y = input(5);
    sigma = input(6);
    % set anonymous functions
    oddfun  = @(z,a,c) c-z*tan(a*z);
    evenfun = @(z,a,c) z+c*tan(a*z);
    % set the number of r.v. in each 1d direction
    norv1d = norv2d;
    % invert correlation lengths
    cx = 1/correl_x;
    cy = 1/correl_y;
    % initialization
    omega_x = zeros(1, norv1d);
    omega_y = zeros(1, norv1d);
    lambda_x = zeros(1, norv1d);
    lambda_y = zeros(1, norv1d);
    alpha_x = zeros(1, norv1d);
    alpha_y = zeros(1, norv1d);
    ef_x = cell(1, norv1d);
    ef_y = cell(1, norv1d);
    for n = 1:norv1d
        if n==1
            x0 = [0, pi/(2*ax) - 1.0e-08];
            y0 = [0, pi/(2*ay) - 1.0e-08];
        else
            k = floor(n/2);
            x0 = [(2*k-1)*pi/(2*ax) + 1.0e-08, (2*k+1)*pi/(2*ax) - 1.0e-08];
            y0 = [(2*k-1)*pi/(2*ay) + 1.0e-08, (2*k+1)*pi/(2*ay) - 1.0e-08];
        end
        if mod(n,2) == 1
            omega_x(n) = fzero(@(z,a,c) oddfun(z,ax,cx), x0);
            omega_y(n) = fzero(@(z,a,c) oddfun(z,ay,cy), y0);
        else
            omega_x(n) = fzero(@(z,a,c) evenfun(z,ax,cx), x0);
            omega_y(n) = fzero(@(z,a,c) evenfun(z,ay,cy), y0);
        end
        lambda_x(n) = 2*cx/(omega_x(n)^2 + cx^2);
        lambda_y(n) = 2*cy/(omega_y(n)^2 + cy^2);
        alpha_x(n) = 1/sqrt(ax + (-1)^(n-1)*sin(2*ax*omega_x(n))/(2*omega_x(n)));
        alpha_y(n) = 1/sqrt(ay + (-1)^(n-1)*sin(2*ay*omega_y(n))/(2*omega_y(n)));
        if mod(n,2) == 1
            ef_x{n} = @(x) alpha_x(n)*(cos(omega_x(n)*x));
            ef_y{n} = @(y) alpha_y(n)*(cos(omega_y(n)*y));
        else
            ef_x{n} = @(x) alpha_x(n)*(sin(omega_x(n)*x));
            ef_y{n} = @(y) alpha_y(n)*(sin(omega_y(n)*y));
        end
    end
    % find 2d eigenvalues, sort them out keeping track of their 'directional' indices
    eigaux_2d = lambda_x' * lambda_y;
    [eigaux_2d_sorted, ind_sorted] = sort(eigaux_2d(:)','descend');
    [indaux_x, indaux_y] = ind2sub([norv1d norv1d],ind_sorted);
    % assign values to the fields of the structure KL_DATA
    KL_coeff= ones(size(x)); % a_0
    for m = 1:norv2d
        ef_x_temp = ef_x{indaux_x(m)};
        ef_y_temp = ef_y{indaux_y(m)};
        ev2d_temp = eigaux_2d_sorted(m);
        KL_coeff = KL_coeff + sigma*sqrt(ev2d_temp)*ef_x_temp(x).*ef_y_temp(y) * yy(m); % a_m
    end
end