function [aa, aax1, aax2] = stochcol_diffusion_grad_and_coeff(a, ax1, ax2)
aa = @(x1, x2, yy) a(x1, x2, yy);
aax1 = @(x1, x2, yy) ax1(x1, x2, yy);
aax2 = @(x1, x2, yy) ax2(x1, x2, yy);
end