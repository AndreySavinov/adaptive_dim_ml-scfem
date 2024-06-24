function KL_gradcoeff = stochcol_diffusion_grad_x2_spatial_expansion(x1, x2, yy, input)
M = size(yy, 2);;
alpha_bar = input(2);
sigma_tilde = input(3);
KL_gradcoeff = zeros(size(x2)); % a_0_dy
for m = 1:M
    km = floor(-0.5e0+sqrt(0.25e0+2*m));
    beta_x1 = m -km*(km + 1)/2;
    beta_x2 = km -beta_x1;
    KL_gradcoeff= KL_gradcoeff -2*pi*beta_x2*alpha_bar/(m^(sigma_tilde))*cos(2*pi*beta_x1*x1).*sin(2*pi*beta_x2*x2)*yy(m);
end
end
