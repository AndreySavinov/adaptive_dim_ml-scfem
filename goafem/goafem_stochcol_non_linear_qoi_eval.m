function Q = goafem_stochcol_non_linear_qoi_eval(u_gal, paras_fem, qoi)
% works only in spatial domain
% return Q as a square matrix of size number of collocation nodes
M = stochcol_mass_matrix(paras_fem{1}, paras_fem{2});
Q = u_gal.'*M*u_gal;
